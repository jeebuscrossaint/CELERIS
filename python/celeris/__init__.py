"""CELERIS — GPU-ready metalens / metasurface design via RCWA.

Thin Python surface over the validated C++ engine (the same one the desktop GUI
and CLI use). The compiled extension `_celeris` is imported and re-exported here.

Typical workflow (see examples/python/)::

    import celeris as cel

    tio2 = cel.Material.constant(2.40 + 0.0j, "TiO2")
    lib = cel.build_unit_cell_library(
        pillar=tio2, background=cel.materials.air(),
        incident=cel.materials.air(), substrate=cel.materials.fused_silica(),
        period_um=0.35, wavelength_um=0.532, thickness_um=0.6,
        fill_min=0.1, fill_max=0.9, n_samples=24, M=8)
    lens = cel.design_metalens(lib, focal_length_um=50.0, diameter_um=20.0)
    foc = cel.analyze_focus(lens, lib, 50.0, 0.532, 20.0)
    print("Strehl", foc.strehl, "FWHM", foc.fwhm_um, "um")
"""

from ._celeris import (  # noqa: F401  (re-export)
    Pol,
    MetaShape,
    PhaseProfileKind,
    Material,
    materials,
    BinaryGrating1D,
    GratingLayer1D,
    Rcwa1DStack,
    Rcwa1DResult,
    solve_rcwa_1d,
    RectCell2D,
    Rcwa2DStack,
    Rcwa2DResult,
    solve_rcwa_2d,
    UnitCellLibrary,
    build_unit_cell_library,
    HeightSweepEntry,
    HeightOptResult,
    optimize_height_for_2pi,
    MetalensDesign,
    design_metalens,
    PhaseProfile,
    phase_profile_value,
    load_freeform_phase,
    FocalAnalysis,
    analyze_focus,
    PsfMap,
    compute_psf,
    JonesMatrix,
    solve_jones,
    HwpAtom,
    find_hwp_atom,
    PbVerifyPoint,
    verify_pb_phase,
    PbMetalensDesign,
    design_pb_metalens,
    DispersiveAtom,
    MetaAtomSpec,
    DispersiveLibrary,
    build_dispersive_library,
    build_dispersive_library_from_specs,
    build_single_etch_library,
    AchromaticDesign,
    design_achromatic_metalens,
    AchromaticFocalPoint,
    verify_achromatic_focus,
    to_metalens_design,
)

__version__ = "0.1.0"

__all__ = [
    "Pol",
    "MetaShape",
    "PhaseProfileKind",
    "Material",
    "materials",
    "BinaryGrating1D",
    "GratingLayer1D",
    "Rcwa1DStack",
    "Rcwa1DResult",
    "solve_rcwa_1d",
    "RectCell2D",
    "Rcwa2DStack",
    "Rcwa2DResult",
    "solve_rcwa_2d",
    "UnitCellLibrary",
    "build_unit_cell_library",
    "HeightSweepEntry",
    "HeightOptResult",
    "optimize_height_for_2pi",
    "MetalensDesign",
    "design_metalens",
    "PhaseProfile",
    "phase_profile_value",
    "load_freeform_phase",
    "FocalAnalysis",
    "analyze_focus",
    "PsfMap",
    "compute_psf",
    "JonesMatrix",
    "solve_jones",
    "HwpAtom",
    "find_hwp_atom",
    "PbVerifyPoint",
    "verify_pb_phase",
    "PbMetalensDesign",
    "design_pb_metalens",
    "DispersiveAtom",
    "MetaAtomSpec",
    "DispersiveLibrary",
    "build_dispersive_library",
    "build_dispersive_library_from_specs",
    "build_single_etch_library",
    "AchromaticDesign",
    "design_achromatic_metalens",
    "AchromaticFocalPoint",
    "verify_achromatic_focus",
    "to_metalens_design",
]
