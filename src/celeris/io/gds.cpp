#include "celeris/io/gds.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <vector>

namespace celeris {
namespace {

// GDSII record types we emit.
enum Rec : uint8_t {
    HEADER = 0x00, BGNLIB = 0x01, LIBNAME = 0x02, UNITS = 0x03, ENDLIB = 0x04,
    BGNSTR = 0x05, STRNAME = 0x06, ENDSTR = 0x07, BOUNDARY = 0x08,
    LAYER = 0x0D, DATATYPE = 0x0E, XY = 0x10, ENDEL = 0x11,
};
// GDSII data types.
enum Dt : uint8_t {
    NODATA = 0x00, BITARRAY = 0x01, INT16 = 0x02, INT32 = 0x03, REAL8 = 0x05,
    ASCII = 0x06,
};

void put16(std::vector<uint8_t>& b, uint16_t v) {  // big-endian
    b.push_back(static_cast<uint8_t>(v >> 8));
    b.push_back(static_cast<uint8_t>(v & 0xFF));
}
void put32(std::vector<uint8_t>& b, int32_t v) {  // big-endian, two's complement
    uint32_t u = static_cast<uint32_t>(v);
    b.push_back(static_cast<uint8_t>(u >> 24));
    b.push_back(static_cast<uint8_t>(u >> 16));
    b.push_back(static_cast<uint8_t>(u >> 8));
    b.push_back(static_cast<uint8_t>(u));
}

// Encode an IEEE double as a GDSII 8-byte real (sign, 7-bit excess-64 base-16
// exponent, 56-bit mantissa).
void put_real8(std::vector<uint8_t>& b, double value) {
    uint8_t out[8] = {0};
    if (value != 0.0) {
        uint8_t sign = 0;
        if (value < 0) { sign = 0x80; value = -value; }
        int exp = 64;  // excess-64
        while (value >= 1.0) { value /= 16.0; ++exp; }
        while (value < 1.0 / 16.0) { value *= 16.0; --exp; }
        // value now in [1/16, 1); mantissa = value * 2^56.
        long double mant = static_cast<long double>(value);
        out[0] = static_cast<uint8_t>(sign | (exp & 0x7F));
        for (int i = 1; i < 8; ++i) {
            mant *= 256.0L;
            uint8_t byte = static_cast<uint8_t>(mant);
            out[i] = byte;
            mant -= byte;
        }
    }
    for (uint8_t v : out) b.push_back(v);
}

// Write one record: 2-byte total length, type, datatype, then payload.
void record(std::ofstream& f, Rec type, Dt dt, const std::vector<uint8_t>& data) {
    uint16_t len = static_cast<uint16_t>(4 + data.size());
    std::vector<uint8_t> hdr;
    put16(hdr, len);
    hdr.push_back(static_cast<uint8_t>(type));
    hdr.push_back(static_cast<uint8_t>(dt));
    f.write(reinterpret_cast<const char*>(hdr.data()),
            static_cast<std::streamsize>(hdr.size()));
    if (!data.empty())
        f.write(reinterpret_cast<const char*>(data.data()),
                static_cast<std::streamsize>(data.size()));
}

std::vector<uint8_t> i16(int16_t v) { std::vector<uint8_t> d; put16(d, static_cast<uint16_t>(v)); return d; }
std::vector<uint8_t> ascii(std::string s) {
    if (s.size() & 1) s.push_back('\0');  // GDSII strings are even-length
    return std::vector<uint8_t>(s.begin(), s.end());
}
std::vector<uint8_t> timestamps() { return std::vector<uint8_t>(24, 0); }  // 12 INT16 zeros

} // namespace

int write_metalens_gds(const MetalensDesign& lens, const std::string& path,
                       int layer, double min_fill) {
    std::ofstream f(path, std::ios::binary);
    if (!f) return -1;

    record(f, HEADER, INT16, i16(600));            // GDSII v6
    record(f, BGNLIB, INT16, timestamps());
    record(f, LIBNAME, ASCII, ascii("CELERIS"));
    { std::vector<uint8_t> u; put_real8(u, 1e-3); put_real8(u, 1e-9); // db unit = 1 nm, in meters
      record(f, UNITS, REAL8, u); }
    record(f, BGNSTR, INT16, timestamps());
    record(f, STRNAME, ASCII, ascii("METALENS"));

    const double p = lens.period_um;
    const double center = (lens.n_cells - 1) / 2.0;
    const double to_db = 1000.0;  // microns -> db units (nm)
    int count = 0;

    for (int iy = 0; iy < lens.n_cells; ++iy) {
        for (int ix = 0; ix < lens.n_cells; ++ix) {
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * lens.n_cells + ix];
            if (fill < min_fill) continue;
            const double cx = (ix - center) * p;
            const double cy = (iy - center) * p;
            const double half = 0.5 * fill * p;
            const int32_t x0 = static_cast<int32_t>(std::lround((cx - half) * to_db));
            const int32_t x1 = static_cast<int32_t>(std::lround((cx + half) * to_db));
            const int32_t y0 = static_cast<int32_t>(std::lround((cy - half) * to_db));
            const int32_t y1 = static_cast<int32_t>(std::lround((cy + half) * to_db));

            record(f, BOUNDARY, NODATA, {});
            record(f, LAYER, INT16, i16(static_cast<int16_t>(layer)));
            record(f, DATATYPE, INT16, i16(0));
            std::vector<uint8_t> xy;  // closed polygon: 5 points
            put32(xy, x0); put32(xy, y0);
            put32(xy, x1); put32(xy, y0);
            put32(xy, x1); put32(xy, y1);
            put32(xy, x0); put32(xy, y1);
            put32(xy, x0); put32(xy, y0);
            record(f, XY, INT32, xy);
            record(f, ENDEL, NODATA, {});
            ++count;
        }
    }

    record(f, ENDSTR, NODATA, {});
    record(f, ENDLIB, NODATA, {});
    return f.good() ? count : -1;
}

int gds_count_boundaries(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return -1;
    int count = 0;
    bool saw_endlib = false;
    while (f) {
        uint8_t h[4];
        if (!f.read(reinterpret_cast<char*>(h), 4)) break;
        uint16_t len = static_cast<uint16_t>((h[0] << 8) | h[1]);
        if (len < 4) return -1;  // malformed
        uint8_t type = h[2];
        if (type == BOUNDARY) ++count;
        if (type == ENDLIB) saw_endlib = true;
        f.seekg(len - 4, std::ios::cur);  // skip payload
    }
    return saw_endlib ? count : -1;
}

} // namespace celeris
