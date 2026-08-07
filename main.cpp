// CELERIS command-line front-end: a thin dispatcher over the subcommands in
// cli/*.cpp. Shared helpers and all subcommand declarations live in cli/cli.hpp.
// Run `celeris help` for usage.

#include <string>

#include "cli/cli.hpp"

int main(int argc, char** argv) {
    const std::string cmd = argc > 1 ? argv[1] : "help";
    if (cmd == "selftest") {
        bool quick = false;
        for (int i = 2; i < argc; ++i)
            if (std::string(argv[i]) == "--quick") quick = true;
        return run_selftest(quick);
    }
    if (cmd == "shapeconv") return cmd_shapeconv(argc, argv);
    if (cmd == "validate") return cmd_validate(argc, argv);
    if (cmd == "reproduce") return cmd_reproduce(argc, argv);
    if (cmd == "materials") return cmd_materials(argc, argv);
    if (cmd == "design") return cmd_design(argc, argv);
    if (cmd == "efficiency") return cmd_efficiency(argc, argv);
    if (cmd == "fieldmap") return cmd_fieldmap(argc, argv);
    if (cmd == "widefov") return cmd_widefov(argc, argv);
    if (cmd == "birefringence") return cmd_birefringence(argc, argv);
    if (cmd == "polardesign") return cmd_polardesign(argc, argv);
    if (cmd == "pbdesign") return cmd_pbdesign(argc, argv);
    if (cmd == "achromatic") return cmd_achromatic(argc, argv);
    if (cmd == "pbachromatic") return cmd_pb_achromatic(argc, argv);
#ifdef CELERIS_USE_CUDA
    if (cmd == "gpubench") return run_gpubench(argc, argv);
#endif
#ifdef CELERIS_USE_CUDA_KERNELS
    if (cmd == "psfbench") return run_psfbench(argc, argv);
#endif
    print_help();
    return (cmd == "help" || cmd == "--help" || cmd == "-h") ? 0 : 1;
}
