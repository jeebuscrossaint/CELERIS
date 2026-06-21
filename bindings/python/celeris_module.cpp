// Python bindings for the CELERIS engine (celeris_core), via pybind11.
//
// Goal: let the metasurface/CREOL community drive the SAME validated C++ engine
// the GUI and CLI use, from a notebook -- `import celeris`, build a material,
// sweep a pillar library, design a focusing metalens, analyze the focus, and get
// PSF/efficiency maps back as numpy arrays for matplotlib. No physics lives here:
// this file is a thin, faithful surface over the public headers in src/celeris.
//
// The module is named `_celeris` (the compiled extension); the `celeris` Python
// package (python/celeris/__init__.py) re-exports it with a friendlier surface.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>       // std::vector <-> list
#include <pybind11/complex.h>   // std::complex<double> <-> complex
#include <pybind11/numpy.h>     // py::array_t for the 2D maps

#include "celeris/core.hpp"
#include "celeris/materials/material.hpp"
#include "celeris/materials/database.hpp"
#include "celeris/rcwa/grating1d.hpp"
#include "celeris/rcwa/rcwa1d.hpp"
#include "celeris/rcwa/rcwa2d.hpp"
#include "celeris/design/metalens.hpp"
#include "celeris/design/phase_profile.hpp"
#include "celeris/design/pb_metalens.hpp"
#include "celeris/analysis/focal.hpp"

#include <stdexcept>

namespace py = pybind11;
using namespace celeris;

// Copy a row-major n_rows x n_cols vector into a freshly-owned numpy 2D array.
static py::array_t<double> to_2d_array(const std::vector<double>& v, int n_rows,
                                       int n_cols) {
    py::array_t<double> a({n_rows, n_cols});
    std::copy(v.begin(), v.end(), a.mutable_data());
    return a;
}

// Flatten a square 2D numpy array (row-major) into a vector<double>, returning the
// side length n. Used to load a freeform phase map from numpy.
static std::vector<double> from_square_2d(const py::array_t<double>& arr, int& n) {
    auto buf = py::array_t<double, py::array::c_style | py::array::forcecast>(arr);
    if (buf.ndim() != 2 || buf.shape(0) != buf.shape(1))
        throw std::runtime_error("freeform phase map must be a square 2D array");
    n = static_cast<int>(buf.shape(0));
    const double* p = buf.data();
    return std::vector<double>(p, p + static_cast<size_t>(n) * n);
}

PYBIND11_MODULE(_celeris, m) {
    m.doc() = "CELERIS — GPU-ready metalens / metasurface design via rigorous "
              "coupled-wave analysis. Bindings over the validated C++ engine.";

    // ---- enums -------------------------------------------------------------
    py::enum_<Pol>(m, "Pol", "Polarization relative to the plane of incidence.")
        .value("TE", Pol::TE, "s-polarization (E perpendicular to plane of incidence)")
        .value("TM", Pol::TM, "p-polarization (E in the plane of incidence)");

    py::enum_<MetaShape>(m, "MetaShape", "Meta-atom cross-section shape.")
        .value("Rectangle", MetaShape::Rectangle)
        .value("Ellipse", MetaShape::Ellipse)
        .value("Cross", MetaShape::Cross)
        .value("Ring", MetaShape::Ring);

    py::enum_<PhaseProfileKind>(m, "PhaseProfileKind",
        "Target wavefront phi(x,y) a metasurface imprints.")
        .value("Focusing", PhaseProfileKind::Focusing, "hyperbolic lens -> focal spot")
        .value("Vortex", PhaseProfileKind::Vortex, "OAM plate -> donut / OAM beam")
        .value("Deflector", PhaseProfileKind::Deflector, "blazed grating -> tilted beam")
        .value("Axicon", PhaseProfileKind::Axicon, "conical phase -> Bessel / line focus")
        .value("Freeform", PhaseProfileKind::Freeform,
               "arbitrary loaded phi(x,y) map (hologram / CGH), bilinearly sampled");

    // ---- Material ----------------------------------------------------------
    py::class_<Material>(m, "Material",
        "Optical material: complex refractive index n+ik as a function of vacuum "
        "wavelength (microns). Construct via the constant/sellmeier/tabulated "
        "factories.")
        .def_static("constant", &Material::constant, py::arg("n"),
                    py::arg("name") = "constant",
                    "Non-dispersive material with a fixed complex index (use a "
                    "Python complex like 2.4+0.01j to add loss).")
        .def_static("sellmeier", &Material::sellmeier, py::arg("bc_terms"),
                    py::arg("name") = "sellmeier",
                    "Sellmeier dispersion (lossless): n^2 = 1 + sum_i B_i*l^2/(l^2-C_i), "
                    "l in microns. bc_terms is a list of [B_i, C_i] pairs.")
        .def_static("tabulated", &Material::tabulated, py::arg("wavelength_um"),
                    py::arg("n"), py::arg("k"), py::arg("name") = "tabulated",
                    "Tabulated n,k vs wavelength (microns, strictly ascending), "
                    "linearly interpolated; out-of-range queries clamp.")
        .def("index", &Material::index, py::arg("wavelength_um"),
             "Complex refractive index n+ik at the given wavelength (microns).")
        .def("permittivity", &Material::permittivity, py::arg("wavelength_um"),
             "Relative permittivity eps = n^2 at the given wavelength (microns).")
        .def_property_readonly("name", &Material::name)
        .def("__repr__", [](const Material& mat) {
            return "<celeris.Material '" + mat.name() + "'>";
        });

    // Built-in materials (exact, published dispersion models). Returned by value
    // (Material is value-semantic / copyable), so callers own a stable copy.
    auto mats = m.def_submodule("materials", "Built-in published materials.");
    mats.def("air", [] { return materials::air(); }, "n = 1 exactly.");
    mats.def("bk7", [] { return materials::bk7(); }, "Schott N-BK7 (Sellmeier).");
    mats.def("fused_silica", [] { return materials::fused_silica(); },
             "Fused silica SiO2 (Malitson 1965).");
    mats.def("silicon_nitride", [] { return materials::silicon_nitride(); },
             "Si3N4 LPCVD (Luke 2015).");

    // ---- 1D RCWA -----------------------------------------------------------
    py::class_<BinaryGrating1D>(m, "BinaryGrating1D",
        "A 1D binary (lamellar) grating layer: a ridge of width fill*period in a "
        "groove, periodic in x.")
        .def(py::init([](Material ridge, Material groove, double period_um,
                         double fill, double thickness_um) {
                 return BinaryGrating1D{std::move(ridge), std::move(groove),
                                        period_um, fill, thickness_um};
             }),
             py::arg("ridge"), py::arg("groove"), py::arg("period_um"),
             py::arg("fill"), py::arg("thickness_um"))
        .def_readwrite("ridge", &BinaryGrating1D::ridge)
        .def_readwrite("groove", &BinaryGrating1D::groove)
        .def_readwrite("period_um", &BinaryGrating1D::period_um)
        .def_readwrite("fill", &BinaryGrating1D::fill)
        .def_readwrite("thickness_um", &BinaryGrating1D::thickness_um);

    py::class_<GratingLayer1D>(m, "GratingLayer1D",
        "One layer of a 1D multilayer stack (shares the stack's period).")
        .def(py::init([](Material ridge, Material groove, double fill,
                         double thickness_um) {
                 return GratingLayer1D{std::move(ridge), std::move(groove), fill,
                                       thickness_um};
             }),
             py::arg("ridge"), py::arg("groove"), py::arg("fill"),
             py::arg("thickness_um"))
        .def_static("homogeneous", &GratingLayer1D::homogeneous, py::arg("material"),
                    py::arg("thickness_um"), "A uniform (unpatterned) layer.")
        .def_readwrite("ridge", &GratingLayer1D::ridge)
        .def_readwrite("groove", &GratingLayer1D::groove)
        .def_readwrite("fill", &GratingLayer1D::fill)
        .def_readwrite("thickness_um", &GratingLayer1D::thickness_um);

    py::class_<Rcwa1DStack>(m, "Rcwa1DStack",
        "A 1D multilayer stack: a shared period and an ordered list of layers "
        "(top/incident side first).")
        .def(py::init([](double period_um, std::vector<GratingLayer1D> layers) {
                 return Rcwa1DStack{period_um, std::move(layers)};
             }),
             py::arg("period_um"), py::arg("layers"))
        .def_readwrite("period_um", &Rcwa1DStack::period_um)
        .def_readwrite("layers", &Rcwa1DStack::layers);

    py::class_<Rcwa1DResult>(m, "Rcwa1DResult",
        "Per-order 1D diffraction efficiencies.")
        .def_readonly("orders", &Rcwa1DResult::orders)
        .def_readonly("de_r", &Rcwa1DResult::de_r, "reflected DE per order")
        .def_readonly("de_t", &Rcwa1DResult::de_t, "transmitted DE per order")
        .def_readonly("sum_de", &Rcwa1DResult::sum_de, "R+T (==1 if lossless)");

    m.def("solve_rcwa_1d",
          py::overload_cast<const Material&, const BinaryGrating1D&, const Material&,
                            double, double, int, Pol>(&solve_rcwa_1d),
          py::arg("incident"), py::arg("grating"), py::arg("substrate"),
          py::arg("wavelength_um"), py::arg("theta0_rad"), py::arg("n_harmonics"),
          py::arg("pol"), "Solve a single 1D binary grating layer.");
    m.def("solve_rcwa_1d",
          py::overload_cast<const Material&, const Rcwa1DStack&, const Material&,
                            double, double, int, Pol>(&solve_rcwa_1d),
          py::arg("incident"), py::arg("stack"), py::arg("substrate"),
          py::arg("wavelength_um"), py::arg("theta0_rad"), py::arg("n_harmonics"),
          py::arg("pol"), "Solve a 1D multilayer stack (S-matrix recursion).");

    // ---- 2D RCWA -----------------------------------------------------------
    py::class_<RectCell2D>(m, "RectCell2D",
        "One biperiodic layer: a pillar (shape inside the fill_x x fill_y box) of "
        "`pillar` material embedded in `background`.")
        .def(py::init([](Material pillar, Material background, double fill_x,
                         double fill_y, double thickness_um, MetaShape shape,
                         double shape_param) {
                 return RectCell2D{std::move(pillar), std::move(background), fill_x,
                                   fill_y, thickness_um, shape, shape_param};
             }),
             py::arg("pillar"), py::arg("background"), py::arg("fill_x"),
             py::arg("fill_y"), py::arg("thickness_um"),
             py::arg("shape") = MetaShape::Rectangle, py::arg("shape_param") = 0.5)
        .def_static("homogeneous", &RectCell2D::homogeneous, py::arg("material"),
                    py::arg("thickness_um"))
        .def_readwrite("pillar", &RectCell2D::pillar)
        .def_readwrite("background", &RectCell2D::background)
        .def_readwrite("fill_x", &RectCell2D::fill_x)
        .def_readwrite("fill_y", &RectCell2D::fill_y)
        .def_readwrite("thickness_um", &RectCell2D::thickness_um)
        .def_readwrite("shape", &RectCell2D::shape)
        .def_readwrite("shape_param", &RectCell2D::shape_param);

    py::class_<Rcwa2DStack>(m, "Rcwa2DStack",
        "A biperiodic stack: periods in x and y, plus an ordered list of layers.")
        .def(py::init([](double period_x_um, double period_y_um,
                         std::vector<RectCell2D> layers) {
                 return Rcwa2DStack{period_x_um, period_y_um, std::move(layers)};
             }),
             py::arg("period_x_um"), py::arg("period_y_um"), py::arg("layers"))
        .def_readwrite("period_x_um", &Rcwa2DStack::period_x_um)
        .def_readwrite("period_y_um", &Rcwa2DStack::period_y_um)
        .def_readwrite("layers", &Rcwa2DStack::layers);

    py::class_<Rcwa2DResult>(m, "Rcwa2DResult", "2D RCWA efficiencies + order-0 amplitudes.")
        .def_readonly("R", &Rcwa2DResult::R)
        .def_readonly("T", &Rcwa2DResult::T)
        .def_readonly("sum_de", &Rcwa2DResult::sum_de)
        .def_readonly("de_t0", &Rcwa2DResult::de_t0, "zeroth-order transmitted DE")
        .def_readonly("de_r0", &Rcwa2DResult::de_r0)
        .def_readonly("tx0", &Rcwa2DResult::tx0, "complex order-0 transmitted Ex")
        .def_readonly("ty0", &Rcwa2DResult::ty0, "complex order-0 transmitted Ey");

    m.def("solve_rcwa_2d", &solve_rcwa_2d, py::arg("incident"), py::arg("stack"),
          py::arg("substrate"), py::arg("wavelength_um"), py::arg("theta_rad"),
          py::arg("phi_rad"), py::arg("Ex0"), py::arg("Ey0"), py::arg("Mx"),
          py::arg("My"),
          "Solve a biperiodic stack. Incident E at order 0 is (Ex0, Ey0); "
          "Mx/My are the harmonic half-counts per direction.");

    // ---- design ------------------------------------------------------------
    py::class_<UnitCellLibrary>(m, "UnitCellLibrary",
        "Phase library: the complex transmission (phase + |t|) vs pillar fill for "
        "one fab process (period/wavelength/height).")
        .def_readonly("period_um", &UnitCellLibrary::period_um)
        .def_readonly("wavelength_um", &UnitCellLibrary::wavelength_um)
        .def_readonly("thickness_um", &UnitCellLibrary::thickness_um)
        .def_readonly("fill", &UnitCellLibrary::fill)
        .def_readonly("phase", &UnitCellLibrary::phase, "transmission phase (rad)")
        .def_readonly("amplitude", &UnitCellLibrary::amplitude, "|t| of order 0")
        .def("phase_span", &UnitCellLibrary::phase_span)
        .def("coverage", &UnitCellLibrary::coverage,
             "Effective phase coverage (rad) = 2*pi minus the largest gap on the circle.")
        .def("lookup", &UnitCellLibrary::lookup, py::arg("target_phase_rad"))
        .def("lookup_weighted", &UnitCellLibrary::lookup_weighted,
             py::arg("target_phase_rad"), py::arg("amplitude_weight"))
        .def("transmission_for_fill", &UnitCellLibrary::transmission_for_fill,
             py::arg("fill"));

    m.def("build_unit_cell_library", &build_unit_cell_library, py::arg("pillar"),
          py::arg("background"), py::arg("incident"), py::arg("substrate"),
          py::arg("period_um"), py::arg("wavelength_um"), py::arg("thickness_um"),
          py::arg("fill_min"), py::arg("fill_max"), py::arg("n_samples"),
          py::arg("M"), py::arg("shape") = MetaShape::Rectangle,
          py::arg("shape_param") = 0.5,
          "Sweep a pillar's fill and tabulate complex transmission (one 2D RCWA "
          "solve per sample, normal incidence, x-polarized).");

    py::class_<HeightSweepEntry>(m, "HeightSweepEntry")
        .def_readonly("thickness_um", &HeightSweepEntry::thickness_um)
        .def_readonly("coverage_deg", &HeightSweepEntry::coverage_deg)
        .def_readonly("mean_transmittance", &HeightSweepEntry::mean_transmittance);

    py::class_<HeightOptResult>(m, "HeightOptResult",
        "Result of the full-2pi pillar-height search.")
        .def_readonly("best_thickness_um", &HeightOptResult::best_thickness_um)
        .def_readonly("coverage_deg", &HeightOptResult::coverage_deg)
        .def_readonly("mean_transmittance", &HeightOptResult::mean_transmittance)
        .def_readonly("reached_target", &HeightOptResult::reached_target)
        .def_readonly("coverage_target_deg", &HeightOptResult::coverage_target_deg)
        .def_readonly("sweep", &HeightOptResult::sweep)
        .def_readonly("best_library", &HeightOptResult::best_library);

    m.def("optimize_height_for_2pi", &optimize_height_for_2pi, py::arg("pillar"),
          py::arg("background"), py::arg("incident"), py::arg("substrate"),
          py::arg("period_um"), py::arg("wavelength_um"), py::arg("thick_lo"),
          py::arg("thick_hi"), py::arg("n_heights"), py::arg("fill_min"),
          py::arg("fill_max"), py::arg("fill_samples"), py::arg("M"),
          py::arg("coverage_target_deg") = 330.0,
          py::arg("shape") = MetaShape::Rectangle, py::arg("shape_param") = 0.5,
          "Sweep pillar height for full-2pi coverage; pick the height clearing the "
          "coverage target at the highest transmittance. Returns the winning library too.");

    py::class_<MetalensDesign>(m, "MetalensDesign", "A fabricable pillar map + fidelity.")
        .def_readonly("n_cells", &MetalensDesign::n_cells)
        .def_readonly("period_um", &MetalensDesign::period_um)
        .def_readonly("rms_phase_error_deg", &MetalensDesign::rms_phase_error_deg)
        .def_readonly("mean_amplitude", &MetalensDesign::mean_amplitude)
        .def_property_readonly("fill_map",
            [](const MetalensDesign& d) {
                return to_2d_array(d.fill_map, d.n_cells, d.n_cells);
            },
            "n_cells x n_cells numpy array of pillar fills (row-major).");

    m.def("design_metalens",
          py::overload_cast<const UnitCellLibrary&, double, double, double>(
              &design_metalens),
          py::arg("lib"), py::arg("focal_length_um"), py::arg("diameter_um"),
          py::arg("amplitude_weight") = 0.25,
          "Design a hyperbolic focusing metalens from a library: at each lattice "
          "site pick the pillar whose phase best matches the ideal profile.");

    // ---- phase profiles (shared by both design paths) ----------------------
    py::class_<PhaseProfile>(m, "PhaseProfile",
        "Target phase profile phi(x,y) a metasurface imprints. The SAME profile "
        "drives both the propagation-phase (design_metalens) and geometric-phase "
        "(design_pb_metalens) paths. Use the static factories or set fields directly.")
        .def(py::init<>())
        .def_readwrite("kind", &PhaseProfile::kind)
        .def_readwrite("focal_length_um", &PhaseProfile::focal_length_um,
                       "Focusing/Vortex: focal length (<=0 => pure OAM for Vortex).")
        .def_readwrite("topological_charge", &PhaseProfile::topological_charge,
                       "Vortex: OAM charge l (winds 2*pi*l around the axis).")
        .def_readwrite("deflect_deg", &PhaseProfile::deflect_deg,
                       "Deflector: beam deflection angle from normal.")
        .def_readwrite("deflect_azimuth_deg", &PhaseProfile::deflect_azimuth_deg,
                       "Deflector: in-plane direction of the deflection.")
        .def_readwrite("axicon_deg", &PhaseProfile::axicon_deg,
                       "Axicon: cone half-angle.")
        .def_readwrite("freeform_n", &PhaseProfile::freeform_n)
        .def_readwrite("freeform_extent_um", &PhaseProfile::freeform_extent_um,
                       "Freeform: full physical width the map spans (centered).")
        .def_property("freeform_phase",
            [](const PhaseProfile& p) {
                return to_2d_array(p.freeform_phase_rad, p.freeform_n, p.freeform_n);
            },
            [](PhaseProfile& p, const py::array_t<double>& arr) {
                int n = 0;
                p.freeform_phase_rad = from_square_2d(arr, n);
                p.freeform_n = n;
            },
            "Freeform: the n x n target phase grid (radians), as a numpy 2D array.")
        // -- factories mirroring the CLI --
        .def_static("focusing",
            [](double focal_length_um) {
                PhaseProfile p;
                p.kind = PhaseProfileKind::Focusing;
                p.focal_length_um = focal_length_um;
                return p;
            }, py::arg("focal_length_um"), "Hyperbolic focusing lens.")
        .def_static("vortex",
            [](int charge, double focal_length_um) {
                PhaseProfile p;
                p.kind = PhaseProfileKind::Vortex;
                p.topological_charge = charge;
                p.focal_length_um = focal_length_um;
                return p;
            }, py::arg("charge"), py::arg("focal_length_um") = 0.0,
            "OAM vortex of charge l (focal_length_um>0 => focused vortex / donut).")
        .def_static("deflector",
            [](double deflect_deg, double azimuth_deg) {
                PhaseProfile p;
                p.kind = PhaseProfileKind::Deflector;
                p.deflect_deg = deflect_deg;
                p.deflect_azimuth_deg = azimuth_deg;
                return p;
            }, py::arg("deflect_deg"), py::arg("azimuth_deg") = 0.0,
            "Blazed-grating beam deflector.")
        .def_static("axicon",
            [](double axicon_deg) {
                PhaseProfile p;
                p.kind = PhaseProfileKind::Axicon;
                p.axicon_deg = axicon_deg;
                return p;
            }, py::arg("axicon_deg"), "Conical (axicon) phase -> Bessel/line focus.")
        .def_static("freeform",
            [](const py::array_t<double>& phase_rad, double extent_um) {
                PhaseProfile p;
                p.kind = PhaseProfileKind::Freeform;
                int n = 0;
                p.freeform_phase_rad = from_square_2d(phase_rad, n);
                p.freeform_n = n;
                p.freeform_extent_um = extent_um;
                return p;
            }, py::arg("phase_rad"), py::arg("extent_um"),
            "Arbitrary loaded phi(x,y) (CGH/hologram) from a square numpy array "
            "(radians) spanning extent_um in x and y, bilinearly sampled.");

    m.def("phase_profile_value", &phase_profile_value, py::arg("profile"),
          py::arg("x"), py::arg("y"), py::arg("wavelength_um"),
          "Target phase phi(x,y) [rad] at a lattice point for the given profile.");

    m.def("load_freeform_phase", &load_freeform_phase, py::arg("path"),
          py::arg("extent_um"),
          "Load a freeform phase map from a whitespace text grid (radians, "
          "perfect-square count -> n x n) spanning extent_um, centered on origin.");

    m.def("design_metalens",
          py::overload_cast<const UnitCellLibrary&, const PhaseProfile&, double,
                            double>(&design_metalens),
          py::arg("lib"), py::arg("profile"), py::arg("diameter_um"),
          py::arg("amplitude_weight") = 0.25,
          "Design a metalens for an ARBITRARY phase profile on the propagation-phase "
          "path: at each site pick the library pillar whose phase best matches "
          "phi(x,y). Carries a finite residual phase error (discrete library).");

    // ---- focal analysis ----------------------------------------------------
    py::class_<FocalAnalysis>(m, "FocalAnalysis", "Standard focal-plane metrics.")
        .def_readonly("strehl", &FocalAnalysis::strehl,
                      "transmission-weighted Strehl (conservative; 1.0=perfect)")
        .def_readonly("phase_strehl", &FocalAnalysis::phase_strehl,
                      "phase-only Strehl (1.0=diffraction-limited wavefront)")
        .def_readonly("fwhm_um", &FocalAnalysis::fwhm_um)
        .def_readonly("diffraction_limit_um", &FocalAnalysis::diffraction_limit_um)
        .def_readonly("encircled_energy", &FocalAnalysis::encircled_energy);

    m.def("analyze_focus", &analyze_focus, py::arg("lens"), py::arg("lib"),
          py::arg("focal_length_um"), py::arg("wavelength_um"), py::arg("diameter_um"),
          "Propagate the design to focus and report Strehl/FWHM/encircled-energy.");

    py::class_<PsfMap>(m, "PsfMap",
        "Focal-plane intensity map (the PSF), n x n over +/- half_window_um.")
        .def_readonly("n", &PsfMap::n)
        .def_readonly("half_window_um", &PsfMap::half_window_um)
        .def_property_readonly("intensity",
            [](const PsfMap& p) { return to_2d_array(p.intensity, p.n, p.n); },
            "n x n numpy array of |E|^2 (row-major).");

    m.def("compute_psf", &compute_psf, py::arg("lens"), py::arg("lib"),
          py::arg("focal_length_um"), py::arg("wavelength_um"), py::arg("diameter_um"),
          py::arg("n"), py::arg("half_window_um"),
          "Full 2D PSF on an n x n window (GPU kernel when built with CUDA, else CPU).");

    // ---- Pancharatnam-Berry (geometric-phase) path -------------------------
    py::class_<JonesMatrix>(m, "JonesMatrix",
        "Zeroth-order 2x2 linear-basis Jones transmission matrix.")
        .def_readonly("xx", &JonesMatrix::xx)
        .def_readonly("xy", &JonesMatrix::xy)
        .def_readonly("yx", &JonesMatrix::yx)
        .def_readonly("yy", &JonesMatrix::yy);

    m.def("solve_jones", &solve_jones, py::arg("incident"), py::arg("stack"),
          py::arg("substrate"), py::arg("wavelength_um"), py::arg("M"),
          "Zeroth-order Jones matrix of a 2D cell at normal incidence (two RCWA solves).");

    py::class_<HwpAtom>(m, "HwpAtom",
        "The half-wave-plate meta-atom PB optics rotates per site (max spin-flip "
        "conversion). Ideal HWP: |t_x|=|t_y|=1, retardance=180 deg.")
        .def_readonly("fill_x", &HwpAtom::fill_x)
        .def_readonly("fill_y", &HwpAtom::fill_y)
        .def_readonly("thickness_um", &HwpAtom::thickness_um)
        .def_readonly("t_x", &HwpAtom::t_x)
        .def_readonly("t_y", &HwpAtom::t_y)
        .def_readonly("retardance_deg", &HwpAtom::retardance_deg)
        .def_readonly("conversion_efficiency", &HwpAtom::conversion_efficiency,
                      "|t_x - t_y|^2 / 4 (spin-flip, in [0,1])");

    m.def("find_hwp_atom", &find_hwp_atom, py::arg("pillar"), py::arg("background"),
          py::arg("incident"), py::arg("substrate"), py::arg("period_um"),
          py::arg("wavelength_um"), py::arg("thickness_um"), py::arg("fill_min"),
          py::arg("fill_max"), py::arg("n_samples"), py::arg("M"),
          "Search a (fill_x, fill_y) grid for the best HWP atom (max spin-flip "
          "conversion efficiency).");

    py::class_<PbVerifyPoint>(m, "PbVerifyPoint",
        "One RCWA-measured sample of the geometric-phase relation (RCP input).")
        .def_readonly("rotation_deg", &PbVerifyPoint::rotation_deg)
        .def_readonly("cross_phase_deg", &PbVerifyPoint::cross_phase_deg,
                      "phase of the spin-flipped output (tracks -2*theta)")
        .def_readonly("conversion_eff", &PbVerifyPoint::conversion_eff)
        .def_readonly("copol_leakage", &PbVerifyPoint::copol_leakage);

    m.def("verify_pb_phase", &verify_pb_phase, py::arg("pillar"),
          py::arg("background"), py::arg("incident"), py::arg("substrate"),
          py::arg("period_um"), py::arg("wavelength_um"), py::arg("atom"),
          py::arg("rotations_rad"), py::arg("M"),
          "RCWA-verify the PB relation: solve the ROTATED atom at each angle and "
          "report the spin-flip phase/efficiency (the proof of the 2*theta map).");

    py::class_<PbMetalensDesign>(m, "PbMetalensDesign",
        "A PB metasurface: one HWP atom rotated per site (geometric phase is exact).")
        .def_readonly("n_cells", &PbMetalensDesign::n_cells)
        .def_readonly("period_um", &PbMetalensDesign::period_um)
        .def_readonly("atom", &PbMetalensDesign::atom)
        .def_readonly("rms_phase_error_deg", &PbMetalensDesign::rms_phase_error_deg,
                      "~0: geometric phase is exact")
        .def_readonly("conversion_efficiency", &PbMetalensDesign::conversion_efficiency,
                      "spin-flip efficiency (focusing-efficiency cap)")
        .def_property_readonly("rotation_deg",
            [](const PbMetalensDesign& d) {
                std::vector<double> deg(d.rotation_rad.size());
                for (size_t i = 0; i < d.rotation_rad.size(); ++i)
                    deg[i] = d.rotation_rad[i] * 180.0 / 3.14159265358979323846;
                return to_2d_array(deg, d.n_cells, d.n_cells);
            },
            "n_cells x n_cells numpy array of per-site atom rotation (degrees).");

    m.def("design_pb_metalens",
          py::overload_cast<const HwpAtom&, double, double, const PhaseProfile&,
                            double, int>(&design_pb_metalens),
          py::arg("atom"), py::arg("period_um"), py::arg("wavelength_um"),
          py::arg("profile"), py::arg("diameter_um"), py::arg("handedness") = 1,
          "Design a PB metasurface for an ARBITRARY phase profile: a fixed HWP atom "
          "rotated per site so the cross-circular output carries phi(x,y). "
          "handedness=+1 for RCP illumination (cross phase = -2*theta), -1 for LCP.");

    m.def("design_pb_metalens",
          py::overload_cast<const HwpAtom&, double, double, double, double, int>(
              &design_pb_metalens),
          py::arg("atom"), py::arg("period_um"), py::arg("wavelength_um"),
          py::arg("focal_length_um"), py::arg("diameter_um"),
          py::arg("handedness") = 1,
          "Convenience overload: a PB focusing lens of the given focal length.");
}
