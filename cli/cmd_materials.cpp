#include "cli.hpp"

// celeris materials [--wavelength <µm>]: list the built-in material library
// (analytic dielectrics + real tabulated n,k), with the complex index sampled
// at a reference wavelength so a user can sanity-check each entry at a glance.
int cmd_materials(int argc, char** argv) {
    const double lambda = std::atof(arg_value(argc, argv, "--wavelength", "0.532"));
    std::println("CELERIS material library  (n + ik sampled at λ = {} µm)", lambda);
    std::println("  any of these names work with --pillar / --substrate; tabulated\n"
                 "  entries carry real loss (k), analytic dielectrics are lossless.\n");
    std::println("  {:<14} {:<6} {:>16} {:>12}   {}", "name", "type", "n+ik @λ",
                 "valid (µm)", "source");
    std::println("  {}", std::string(78, '-'));
    for (const auto& m : materials::catalog()) {
        std::string idx = "  (file missing)";
        std::string range = "broadband";
        if (m.lambda_min_um > 0.0 || m.lambda_max_um > 0.0)
            range = std::format("{:.2f}-{:.2f}", m.lambda_min_um, m.lambda_max_um);
        if (!m.tabulated || m.available) {
            try {
                cdouble n = materials::by_name(m.name).index(lambda);
                idx = std::format("{:.3f}{:+.3f}i", n.real(), n.imag());
            } catch (const std::exception&) { idx = "  (error)"; }
        }
        std::println("  {:<14} {:<6} {:>16} {:>12}   {}", m.name,
                     m.tabulated ? "data" : "model", idx, range, m.description);
    }
    std::println("\n  load any other refractiveindex.info CSV with --pillar-csv <file>.");
    return 0;
}
