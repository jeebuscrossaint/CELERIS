# Cover letter — CELERIS submission to Computer Physics Communications

Dear Editors,

Please consider the enclosed manuscript, *"CELERIS: an integrated, validated RCWA
pipeline for metalens design, analysis, and fabrication layout,"* for publication in
*Computer Physics Communications* as a Computer Programs in Physics (CPiP) paper.

Rigorous coupled-wave analysis (RCWA / the Fourier modal method) is the standard rigorous
method for metasurface unit cells, and several open-source solvers exist. CELERIS is not
another solver kernel: it is the first from-scratch, validated, **end-to-end pipeline** that
takes a metalens specification through solver, meta-atom library, phase mapping and design
(including achromatic and Pancharatnam–Berry recipes), a full optical analysis battery, and
fabrication-ready GDSII layout — a workflow that existing solvers leave the user to
hand-assemble. Its correctness is enforced by a locked self-test suite cross-checked against
closed-form physics, energy conservation, and two independent RCWA codes (grcwa and Stanford
S⁴, which it matches to 7×10⁻⁸ in zeroth-order transmittance), and the complete pipeline is
demonstrated by reproducing two canonical published TiO₂ metalenses end to end. We believe
this integration is of direct, broad utility to the metasurface and computational-optics
community that CPC serves.

The software is open source (Apache-2.0), archived with a permanent DOI
(10.5281/zenodo.21891345), and continuously tested on Linux. All numerical results, figures,
and tables in the manuscript are reproducible from the provided scripts.

This manuscript is original, has not been published previously, and is not under
consideration for publication elsewhere. The sole author is the corresponding author and has
approved the submission. We declare no competing interests.

Thank you for your consideration.

Sincerely,
Amarnath Patel
University of Central Florida — am397421@ucf.edu
ORCID 0009-0008-9460-082X
