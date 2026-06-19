#include "textures.hpp"

#include <GLFW/glfw3.h>  // OpenGL (glGenTextures/glTexImage2D, via GLFW's GL include)

#include <algorithm>
#include <cmath>
#include <vector>

using namespace celeris;

namespace celeris::gui {

void upload_psf_texture(const PsfMap& psf, unsigned int& tex) {
    if (psf.n <= 0 || psf.intensity.empty()) return;
    double mx = 0.0;
    for (double v : psf.intensity) mx = std::max(mx, v);
    if (mx <= 0.0) mx = 1.0;
    std::vector<std::uint8_t> rgba(static_cast<std::size_t>(psf.n) * psf.n * 4);
    for (std::size_t i = 0; i < psf.intensity.size(); ++i) {
        double v = std::pow(psf.intensity[i] / mx, 1.0 / 2.2);  // gamma
        auto b = static_cast<std::uint8_t>(std::lround(std::clamp(v, 0.0, 1.0) * 255.0));
        rgba[i * 4 + 0] = b;
        rgba[i * 4 + 1] = b;
        rgba[i * 4 + 2] = b;
        rgba[i * 4 + 3] = 255;
    }
    if (tex == 0) glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, psf.n, psf.n, 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, rgba.data());
}

void upload_wavefront_texture(const WavefrontAnalysis& wf, unsigned int& tex) {
    if (wf.n <= 0 || wf.opd.empty()) return;
    double amp = 1e-9;
    for (std::size_t i = 0; i < wf.opd.size(); ++i)
        if (wf.mask[i]) amp = std::max(amp, std::abs(wf.opd[i]));
    std::vector<std::uint8_t> rgba(static_cast<std::size_t>(wf.n) * wf.n * 4);
    for (std::size_t i = 0; i < wf.opd.size(); ++i) {
        std::uint8_t r, g, b;
        if (!wf.mask[i]) { r = g = b = 70; }  // outside aperture
        else {
            double v = std::clamp(wf.opd[i] / amp, -1.0, 1.0);  // -1..1
            if (v >= 0) { r = 255; g = b = static_cast<std::uint8_t>(255 * (1 - v)); }
            else        { b = 255; r = g = static_cast<std::uint8_t>(255 * (1 + v)); }
        }
        rgba[i * 4 + 0] = r; rgba[i * 4 + 1] = g; rgba[i * 4 + 2] = b; rgba[i * 4 + 3] = 255;
    }
    if (tex == 0) glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, wf.n, wf.n, 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, rgba.data());
}

void hue_to_rgb(double h, std::uint8_t& r, std::uint8_t& g, std::uint8_t& b) {
    double x = (1.0 - std::abs(std::fmod(h * 6.0, 2.0) - 1.0));
    double rr = 0, gg = 0, bb = 0;
    int seg = static_cast<int>(h * 6.0) % 6;
    switch (seg) {
        case 0: rr = 1; gg = x; break;
        case 1: rr = x; gg = 1; break;
        case 2: gg = 1; bb = x; break;
        case 3: gg = x; bb = 1; break;
        case 4: rr = x; bb = 1; break;
        default: rr = 1; bb = x; break;
    }
    r = static_cast<std::uint8_t>(rr * 255);
    g = static_cast<std::uint8_t>(gg * 255);
    b = static_cast<std::uint8_t>(bb * 255);
}

void upload_layout_texture(const MetalensDesign& d, const UnitCellLibrary& lib,
                           int mode, unsigned int& tex) {
    const int n = d.n_cells;
    if (n <= 0 || d.fill_map.empty()) return;
    std::vector<std::uint8_t> rgba(static_cast<std::size_t>(n) * n * 4);
    for (std::size_t i = 0; i < d.fill_map.size(); ++i) {
        double fill = d.fill_map[i];
        std::uint8_t r, g, b;
        if (mode == 0) {
            cdouble t = lib.transmission_for_fill(fill);
            double h = (std::arg(t) + pi) / (2.0 * pi);  // -pi..pi -> 0..1
            hue_to_rgb(std::clamp(h, 0.0, 1.0), r, g, b);
        } else {
            auto v = static_cast<std::uint8_t>(std::clamp(fill, 0.0, 1.0) * 255.0);
            r = g = b = v;
        }
        rgba[i * 4 + 0] = r; rgba[i * 4 + 1] = g; rgba[i * 4 + 2] = b;
        rgba[i * 4 + 3] = 255;
    }
    if (tex == 0) glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, n, n, 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, rgba.data());
}

void upload_caustic_texture(const ThroughFocus& tf, unsigned int& tex) {
    if (tf.caustic_nx <= 0 || tf.caustic.empty()) return;
    const int w = tf.caustic_nx, h = tf.caustic_nz;
    std::vector<std::uint8_t> rgba(static_cast<std::size_t>(w) * h * 4);
    for (std::size_t i = 0; i < tf.caustic.size(); ++i) {
        double v = std::pow(std::clamp((double)tf.caustic[i], 0.0, 1.0), 1.0 / 1.6);
        auto ch = [](double t) {
            return static_cast<std::uint8_t>(std::clamp(t, 0.0, 1.0) * 255.0);
        };
        rgba[i * 4 + 0] = ch(v * 3.0);          // red first
        rgba[i * 4 + 1] = ch(v * 3.0 - 1.0);    // then green
        rgba[i * 4 + 2] = ch(v * 3.0 - 2.0);    // then blue -> white
        rgba[i * 4 + 3] = 255;
    }
    if (tex == 0) glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE,
                 rgba.data());
}

} // namespace celeris::gui
