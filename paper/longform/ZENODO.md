# Zenodo upload: long-form report

Metadata for uploading the long-form report (`paper.pdf`, built from `paper.tex`) to
Zenodo as a citable report with its own DOI (issue #42). This is a separate record from
the software archive (v1.0.0, 10.5281/zenodo.21891346). Paste the fields below into the
Zenodo "New upload" form.

## Build the PDF

```sh
cd paper/longform
latexmk -pdf paper.tex      # needs elsarticle.cls and elsarticle-num.bst on the TeX path
```

Upload `paper/longform/paper.pdf` only. Do not upload the build files.

## Fields

| Zenodo field | Value |
|---|---|
| Resource type | Publication → Report |
| Title | CELERIS: an integrated, validated RCWA pipeline for metalens design, analysis, and fabrication layout |
| Creator | Patel, Amarnath. Affiliation: University of Central Florida. ORCID 0009-0008-9460-082X |
| Publication date | the upload date |
| Version | v1 |
| Language | English |
| License | Creative Commons Attribution 4.0 International (CC BY 4.0) |
| Keywords | metasurface; metalens; rigorous coupled-wave analysis; Fourier modal method; achromatic metalens; group delay; inverse design; GDSII |

### Description

> Technical report for CELERIS, an open-source (Apache-2.0) C++/Python pipeline for
> metalens design by rigorous coupled-wave analysis (RCWA). It describes the solver
> (1D and 2D vectorial Fourier modal method with Li/Liu–Fan factorization and a
> Redheffer S-matrix stack), the meta-atom library builder, the phase-mapping and
> inverse-design layer (including achromatic and Pancharatnam–Berry designs), the
> optical analysis battery, GPU far-field propagation, and GDSII export. Validation
> covers closed-form physics, energy conservation, and the independent solvers grcwa and
> Stanford S4. The report reproduces two published TiO2 metalenses and computes the
> group-delay limit of single-etch broadband achromats. A sweep of 36 single-etch
> libraries (TiO2, Si3N4 and GaN; propagation and geometric phase; six etch depths) maps
> the fraction of the Presutti–Monticone delay–bandwidth bound they reach (13–68 %,
> about 30–55 % for etch depths of 0.8 µm and more) and shows that span-based aperture
> estimates overstate propagation-phase achromats by about a factor of two.
>
> Source code and the scripts for every number: https://github.com/jeebuscrossaint/CELERIS
> The author used Anthropic Claude to assist with software, the study scripts, and
> drafting; see the declaration in the report.

### Related identifiers

| Relation | Identifier | Resource type |
|---|---|---|
| Documents | 10.5281/zenodo.21891346 | Software |
| Is supplemented by | https://github.com/jeebuscrossaint/CELERIS | Software |

## After the upload

Record the new DOI in `paper/paper.md` (research impact) and in
`memory/celeris-resume.md`. If the report changes later, upload a new version of the same
Zenodo record ("New version"), so the concept DOI stays the same.
