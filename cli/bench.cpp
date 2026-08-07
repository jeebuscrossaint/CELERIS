#include "cli.hpp"

#ifdef CELERIS_USE_CUDA
// Honest head-to-head for the batched GPU eigensolve: the metalens library
// sweep is a batch of independent same-size general-complex eigenproblems, so
// this times CPU-sequential, CPU-parallel (the path the real builder uses), and
// GPU-batched over a representative batch. Usage: celeris gpubench [--n N]
// [--batch B] [--streams S].
int run_gpubench(int argc, char** argv) {
    int n = 242, batch = 32, streams = 4;  // n ~ 2N at M=5 (a real 2D RCWA size)
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&] { return (i + 1 < argc) ? std::atoi(argv[++i]) : 0; };
        if (a == "--n") n = next();
        else if (a == "--batch") batch = next();
        else if (a == "--streams") streams = next();
    }
    if (!cuda::available()) { std::println("gpubench: no CUDA device available"); return 1; }

    const std::size_t nn = static_cast<std::size_t>(n) * n;
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<cdouble> As(static_cast<std::size_t>(batch) * nn);
    for (auto& v : As) v = cdouble{dist(rng), dist(rng)};

    std::println("GPU batched eigensolve benchmark");
    std::println("  batch = {} matrices, {}x{} general complex, streams = {}",
                 batch, n, n, streams);

    auto eig_one = [&](int b) {
        Eigen::Map<const Eigen::MatrixXcd> M(As.data() + static_cast<std::size_t>(b) * nn, n, n);
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(M);
        return ces.eigenvalues()(0);  // touch a result so it isn't optimized away
    };

    // CPU sequential.
    auto t0 = std::chrono::steady_clock::now();
    cdouble sink{0, 0};
    for (int b = 0; b < batch; ++b) sink += eig_one(b);
    auto t1 = std::chrono::steady_clock::now();
    double cpu_seq = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // CPU parallel (same std::async fan-out the unit-cell library builder uses).
    auto t2 = std::chrono::steady_clock::now();
    {
        unsigned hw = std::max(1u, std::thread::hardware_concurrency());
        int workers = std::min<int>(static_cast<int>(hw), batch);
        std::vector<std::future<void>> jobs;
        for (int w = 0; w < workers; ++w)
            jobs.push_back(std::async(std::launch::async, [&, w] {
                for (int b = w; b < batch; b += workers) (void)eig_one(b);
            }));
        for (auto& j : jobs) j.get();
    }
    auto t3 = std::chrono::steady_clock::now();
    double cpu_par = std::chrono::duration<double, std::milli>(t3 - t2).count();

    // GPU batched (warm up device first so init isn't charged to the timing).
    std::vector<cdouble> ws(static_cast<std::size_t>(batch) * n);
    std::vector<cdouble> vrs(static_cast<std::size_t>(batch) * nn);
    cuda::geev_batched(As.data(), n, 1, ws.data(), vrs.data(), 1);
    auto t4 = std::chrono::steady_clock::now();
    bool ok = cuda::geev_batched(As.data(), n, batch, ws.data(), vrs.data(), streams);
    auto t5 = std::chrono::steady_clock::now();
    double gpu = std::chrono::duration<double, std::milli>(t5 - t4).count();

    // Correctness: eigenvalues of matrix 0 must match Eigen (order-independent).
    Eigen::Map<const Eigen::MatrixXcd> M0(As.data(), n, n);
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces0(M0);
    std::vector<cdouble> ec(ces0.eigenvalues().data(), ces0.eigenvalues().data() + n);
    std::vector<cdouble> eg(ws.begin(), ws.begin() + n);
    auto cmp = [](cdouble a, cdouble b) {
        return a.real() != b.real() ? a.real() < b.real() : a.imag() < b.imag();
    };
    std::sort(ec.begin(), ec.end(), cmp);
    std::sort(eg.begin(), eg.end(), cmp);
    double max_diff = 0.0;
    for (int i = 0; i < n; ++i) max_diff = std::max(max_diff, std::abs(ec[i] - eg[i]));

    std::println("");
    std::println("  CPU sequential : {:8.1f} ms   ({:.2f} ms / solve)", cpu_seq, cpu_seq / batch);
    std::println("  CPU parallel   : {:8.1f} ms   ({:.1f}x vs seq)", cpu_par, cpu_seq / cpu_par);
    std::println("  GPU batched    : {:8.1f} ms   ({:.2f}x vs CPU parallel, {:.2f}x vs seq)",
                 gpu, cpu_par / gpu, cpu_seq / gpu);
    std::println("  correctness    : max|d eigenvalue| = {:.2e}  ok={}", max_diff, ok);
    (void)sink;
    return ok ? 0 : 1;
}
#endif

#ifdef CELERIS_USE_CUDA_KERNELS
// Honest head-to-head for the GPU far-field propagation kernel: build a real
// metalens, then compute its focal PSF on the CPU (analyze_focus path, double,
// parallel over cores) vs the GPU kernel (float), comparing the maps and timing.
// Usage: celeris psfbench [--diameter D] [--focal F] [--grid N].
int run_psfbench(int argc, char** argv) {
    double diameter = 120.0, focal = 50.0, wavelength = 0.532, period = 0.35;
    int grid = 161, M = 5, fill_samples = 16;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        auto nd = [&] { return (i + 1 < argc) ? std::atof(argv[++i]) : 0.0; };
        auto ni = [&] { return (i + 1 < argc) ? std::atoi(argv[++i]) : 0; };
        if (a == "--diameter") diameter = nd();
        else if (a == "--focal") focal = nd();
        else if (a == "--grid") grid = ni();
    }
    // cuda::available() lives in the (opt-in) cuSOLVER header; without it the
    // propagation kernel just falls back to CPU internally if there's no device.
#ifdef CELERIS_USE_CUDA
    if (!cuda::available()) { std::println("psfbench: no CUDA device available"); return 1; }
#endif

    std::println("PSF propagation benchmark (building lens...)");
    auto pillar = Material::constant(cdouble{2.40, 0.0}, "TiO2~");
    auto lib = build_unit_cell_library(pillar, materials::air(), materials::air(),
                                       materials::bk7(), period, wavelength, 0.6,
                                       0.08, 0.92, fill_samples, M);
    auto lens = design_metalens(lib, focal, diameter);
    const double dl = wavelength * focal / diameter;
    const double W = std::max(6.0 * dl, 4.0);

    // Aperture pillar list (same construction analyze_focus/compute_psf use).
    const double p = lens.period_um;
    const double center = (lens.n_cells - 1) / 2.0;
    const double R_ap = diameter / 2.0;
    std::vector<double> px, py;
    std::vector<cdouble> pt;
    for (int iy = 0; iy < lens.n_cells; ++iy)
        for (int ix = 0; ix < lens.n_cells; ++ix) {
            double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            double fill = lens.fill_map[(std::size_t)iy * lens.n_cells + ix];
            px.push_back(x); py.push_back(y);
            pt.push_back(lib.transmission_for_fill(fill));
        }
    std::println("  {} x {} cells, {} pillars in aperture, {}x{} focal grid",
                 lens.n_cells, lens.n_cells, px.size(), grid, grid);

    // CPU (double, parallel across cores).
    auto t0 = std::chrono::steady_clock::now();
    auto cpu = compute_psf(lens, lib, focal, wavelength, diameter, grid, W);
    auto t1 = std::chrono::steady_clock::now();
    double cpu_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // GPU (float kernel). Warm up first so device init isn't timed.
    std::vector<double> gpu((std::size_t)grid * grid);
    double k = 2.0 * pi / wavelength;
    cuda::propagate_psf(px.data(), py.data(), pt.data(), (int)px.size(), 0.0, 0.0,
                        focal, k, grid, W, gpu.data());  // warm-up
    auto t2 = std::chrono::steady_clock::now();
    bool ok = cuda::propagate_psf(px.data(), py.data(), pt.data(), (int)px.size(),
                                  0.0, 0.0, focal, k, grid, W, gpu.data());
    auto t3 = std::chrono::steady_clock::now();
    double gpu_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

    // Compare peak-normalized maps.
    double cmax = 0, gmax = 0;
    for (double v : cpu.intensity) cmax = std::max(cmax, v);
    for (double v : gpu) gmax = std::max(gmax, v);
    double maxrel = 0.0;
    if (cmax > 0 && gmax > 0)
        for (std::size_t i = 0; i < gpu.size(); ++i)
            maxrel = std::max(maxrel, std::abs(cpu.intensity[i] / cmax - gpu[i] / gmax));

    std::println("");
    std::println("  CPU (double, parallel) : {:8.1f} ms", cpu_ms);
    std::println("  GPU (float kernel)     : {:8.1f} ms   ({:.1f}x faster)", gpu_ms, cpu_ms / gpu_ms);
    std::println("  agreement              : max|d normalized PSF| = {:.2e}  ok={}", maxrel, ok);
    return ok ? 0 : 1;
}
#endif

// Fast shape-convergence diagnostic: sweep M for one shape and print each row as
// soon as it is computed (flushed), so a heavy eigensolve sweep shows progress
// instead of buffering to exit. Usage: celeris shapeconv [circle|cross|square|ring]
int cmd_shapeconv(int argc, char** argv) {
    const std::string shape = argc > 2 ? argv[2] : "cross";
    const auto tio2 = Material::constant(cdouble{2.45, 0.0}, "TiO2~");
    MetaShape ms = shape == "circle" ? MetaShape::Ellipse
                 : shape == "ring"   ? MetaShape::Ring
                 : shape == "square" ? MetaShape::Rectangle
                                     : MetaShape::Cross;
    double fill = shape == "cross" ? 0.8 : 0.7;
    double param = shape == "cross" ? 0.4 : 0.5;
    std::printf("shape=%s  (Λ=0.35 λ=0.532 TiO2/SiO2)\n  %3s  %10s  %8s  %10s  %9s\n",
                shape.c_str(), "M", "T0", "phase°", "ΣDE", "solve(ms)");
    std::fflush(stdout);
    for (int m : {6, 8, 10, 12}) {
        Rcwa2DStack s{0.35, 0.35, {RectCell2D{tio2, materials::air(), fill, fill, 0.6, ms, param}}};
        auto t0 = std::chrono::steady_clock::now();
        auto r = solve_rcwa_2d(materials::air(), s, materials::fused_silica(), 0.532,
                               0.0, 0.0, 1.0, 0.0, m, m);
        auto t1 = std::chrono::steady_clock::now();
        double ms_ = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::printf("  %3d  %10.4f  %8.1f  %10.6f  %9.0f\n", m, r.de_t0,
                    std::arg(r.tx0) * 180.0 / pi, r.sum_de, ms_);
        std::fflush(stdout);
    }
    return 0;
}
