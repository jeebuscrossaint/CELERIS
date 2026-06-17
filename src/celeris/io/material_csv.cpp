#include "celeris/io/material_csv.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace celeris {

Material load_material_csv(const std::string& path, const std::string& name) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("load_material_csv: cannot open " + path);

    std::vector<double> wl, n, k;
    std::string line;
    while (std::getline(f, line)) {
        // Treat commas and semicolons as whitespace, then parse numbers.
        for (char& c : line)
            if (c == ',' || c == ';' || c == '\t') c = ' ';
        std::istringstream ss(line);
        double a, b, c = 0.0;
        if (!(ss >> a >> b)) continue;  // header/comment/blank -> skip
        ss >> c;                        // optional k
        wl.push_back(a);
        n.push_back(b);
        k.push_back(c);
    }
    if (wl.size() < 2)
        throw std::runtime_error("load_material_csv: <2 valid rows in " + path);

    return Material::tabulated(std::move(wl), std::move(n), std::move(k), name);
}

} // namespace celeris
