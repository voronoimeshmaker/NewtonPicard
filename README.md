# NewtonPicard

Research source code, case configurations and numerical results for Picard,
Newton and spectral-radius-based solver selection in semilinear finite-volume systems.

## Repository and archived release

- Source repository: https://github.com/voronoimeshmaker/NewtonPicard
- Archived release v1.0.0: https://doi.org/10.5281/zenodo.23306169
- Version-independent Zenodo DOI: https://doi.org/10.5281/zenodo.23306168
- Git tag: `v1.0.0`
- Commit referenced by that tag: `cc2ae381c79f10cf8f131d4b4906a6abc202881b`
- Immutable source reference: https://github.com/voronoimeshmaker/NewtonPicard/commit/cc2ae381c79f10cf8f131d4b4906a6abc202881b

The v1.0.0 tag identifies the original archived release. Subsequent local changes
must not be described as part of that release without a corresponding new archive.

## Computational framework

The article campaigns use ADOL-C for source-term differentiation, PETSc for
numerical linear algebra and SLEPc for spectral analysis. The framework also uses
VoronoiMeshMaker and CGAL for mesh generation and VoroGradToolkit for gradient
reconstruction. See [README-ADOLC.md](README-ADOLC.md) for differentiation and build notes.

## Citation

Cite the version-specific DOI above for the archived release and include the
GitHub repository URL for access to the maintained source. Machine-readable
citation metadata are in [CITATION.cff](CITATION.cff).

## Associated manuscript

A spectral criterion for the solution branch reached by the Picard iteration in semilinear unstructured finite-volume systems

## Licence

Code: GPL-3.0-or-later. Author-owned data, case configurations and documentation:
CC-BY-4.0. Existing third-party notices remain applicable. See [LICENSING.md](LICENSING.md)
and [LICENSE](LICENSE).
