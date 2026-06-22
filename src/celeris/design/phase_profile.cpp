#include "celeris/design/phase_profile.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace celeris {
namespace {

// Bilinear sample of the freeform grid at physical (x,y); edges clamped.
double sample_freeform(const PhaseProfile& p, double x, double y) {
    const int n = p.freeform_n;
    if (n <= 0 || p.freeform_extent_um <= 0.0 ||
        (int)p.freeform_phase_rad.size() < n * n)
        return 0.0;
    if (n == 1) return p.freeform_phase_rad[0];
    const double W = p.freeform_extent_um;
    // Physical [-W/2, W/2] -> grid coordinate [0, n-1], clamped at the edges.
    double gx = std::clamp((x / W + 0.5) * (n - 1), 0.0, (double)(n - 1));
    double gy = std::clamp((y / W + 0.5) * (n - 1), 0.0, (double)(n - 1));
    const int x0 = (int)std::floor(gx), y0 = (int)std::floor(gy);
    const int x1 = std::min(x0 + 1, n - 1), y1 = std::min(y0 + 1, n - 1);
    const double fx = gx - x0, fy = gy - y0;
    auto at = [&](int ix, int iy) {
        return p.freeform_phase_rad[(std::size_t)iy * n + ix];
    };
    const double top = at(x0, y0) * (1 - fx) + at(x1, y0) * fx;
    const double bot = at(x0, y1) * (1 - fx) + at(x1, y1) * fx;
    return top * (1 - fy) + bot * fy;
}

} // namespace

double phase_profile_value(const PhaseProfile& p, double x, double y,
                           double wavelength_um) {
    const double k = 2.0 * pi / wavelength_um;
    const double r = std::sqrt(x * x + y * y);
    switch (p.kind) {
        case PhaseProfileKind::Focusing:
            // Converging hyperbolic wavefront -> on-axis focus at z = f.
            return -k * (std::sqrt(r * r + p.focal_length_um * p.focal_length_um) -
                         p.focal_length_um);
        case PhaseProfileKind::Quadratic:
            // Parabolic (paraxial) lens. KEY wide-FOV property: an obliquely
            // incident plane wave adds a linear phase k*sin(theta)*x, and
            //   -k*(x^2+y^2)/(2f) + k*sin(theta)*x
            //     = -k/(2f)*[(x - f*sin(theta))^2 + y^2] + const,
            // i.e. the SAME parabola simply recentered at x0 = f*sin(theta).
            // The focal spot shifts laterally with the field angle but the
            // wavefront aberration is angle-independent -> no coma, a wide field
            // of view. The price (vs the hyperbolic lens, which is perfect only
            // on-axis): on-axis spherical aberration that grows with NA and a
            // curved (Petzval) focal surface, so a real wide-FOV design pairs
            // this with an aperture stop offset from the metasurface.
            return -k * r * r / (2.0 * p.focal_length_um);
        case PhaseProfileKind::Vortex: {
            // Azimuthal phase ramp l*atan2(y,x) carries l*hbar of OAM per photon
            // (winds 2*pi*l around the axis). Add the focusing term for a focused
            // vortex (a donut focal spot); f <= 0 leaves a pure (collimated) OAM beam.
            double phi = p.topological_charge * std::atan2(y, x);
            if (p.focal_length_um > 0.0)
                phi += -k * (std::sqrt(r * r + p.focal_length_um * p.focal_length_um) -
                             p.focal_length_um);
            return phi;
        }
        case PhaseProfileKind::Deflector: {
            // Linear ramp: a transverse momentum k*sin(a) along the azimuth direction
            // tilts a normally-incident beam to angle `a` (a blazed grating / prism).
            const double a = p.deflect_deg * pi / 180.0;
            const double az = p.deflect_azimuth_deg * pi / 180.0;
            return k * std::sin(a) * (x * std::cos(az) + y * std::sin(az));
        }
        case PhaseProfileKind::Axicon: {
            // Conical phase: every ray bends toward the axis by `b`, forming a
            // non-diffracting Bessel beam over an extended on-axis line focus.
            const double b = p.axicon_deg * pi / 180.0;
            return -k * std::sin(b) * r;
        }
        case PhaseProfileKind::Freeform:
            return sample_freeform(p, x, y);
    }
    return 0.0;
}

PhaseProfile load_freeform_phase(const std::string& path, double extent_um) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open freeform phase file: " + path);
    std::vector<double> vals;
    std::string line;
    while (std::getline(f, line)) {
        if (std::size_t h = line.find('#'); h != std::string::npos)
            line.resize(h);  // strip an inline comment
        std::istringstream ss(line);
        for (double v; ss >> v;) vals.push_back(v);
    }
    if (vals.empty())
        throw std::runtime_error("freeform phase file is empty: " + path);
    const int n = (int)std::llround(std::sqrt((double)vals.size()));
    if ((std::size_t)n * n != vals.size())
        throw std::runtime_error(
            "freeform phase file is not a square grid (" +
            std::to_string(vals.size()) + " values, not a perfect square): " + path);
    if (extent_um <= 0.0)
        throw std::runtime_error("freeform phase extent must be positive");

    PhaseProfile p;
    p.kind = PhaseProfileKind::Freeform;
    p.freeform_phase_rad = std::move(vals);
    p.freeform_n = n;
    p.freeform_extent_um = extent_um;
    return p;
}

} // namespace celeris
