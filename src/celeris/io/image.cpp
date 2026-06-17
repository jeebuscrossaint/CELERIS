#include "celeris/io/image.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>

namespace celeris {

bool write_pgm(const std::string& path, int nx, int ny,
               const std::vector<double>& values, double gamma) {
    if (nx <= 0 || ny <= 0 ||
        values.size() != static_cast<std::size_t>(nx) * ny)
        return false;

    double maxv = 0.0;
    for (double v : values) maxv = std::max(maxv, v);
    if (maxv <= 0.0) maxv = 1.0;

    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    f << "P5\n" << nx << " " << ny << "\n255\n";

    const double inv_g = (gamma != 0.0) ? 1.0 / gamma : 1.0;
    std::vector<uint8_t> row(static_cast<std::size_t>(nx));
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            double norm = values[static_cast<std::size_t>(j) * nx + i] / maxv;
            norm = std::clamp(norm, 0.0, 1.0);
            if (gamma != 1.0) norm = std::pow(norm, inv_g);
            row[static_cast<std::size_t>(i)] =
                static_cast<uint8_t>(std::lround(norm * 255.0));
        }
        f.write(reinterpret_cast<const char*>(row.data()), nx);
    }
    return f.good();
}

} // namespace celeris
