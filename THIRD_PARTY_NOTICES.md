# Attribution and software boundaries

The MIT license applies to original source code, configuration examples and authored documentation in this repository. It does not grant rights to COMSOL, DWSIM, their libraries or the publications cited below.

- **COMSOL Multiphysics** is separately licensed software. Users supply a compatible installation and licenses for the physics interfaces in their models. No COMSOL binaries, official model files, manual text or license material are bundled. The Java model builders here were authored for this workflow. Official documentation is linked as API and modeling reference.
- **DWSIM** is a separate project with its own license: [official repository](https://github.com/DanWBR/dwsim). The optional adapter loads a user-installed runtime; no DWSIM assemblies or databases are bundled or relicensed.
- **Research publications** remain subject to their publisher licenses. This repository includes bibliographic metadata, bounded factual parameter summaries and separately labeled assumptions. It does not redistribute full-text papers, screenshots or supplementary datasets. Source-specific access and evidence limits are recorded in the JSON registers.
- **Python** and the Windows/.NET toolchain are supplied by the user. The reaction runners use Python's standard library; no Python packages or compiled dependencies are vendored.

Chemistry references: [Han et al. 2025](https://www.nature.com/articles/s41467-025-58813-3), [Shi et al. 2024](https://www.nature.com/articles/s41467-024-47179-7), and additional sources in [the catalyst register](skills/ai-comsol-skill/references/catalyst-reactor-data-v2.json). A source about a different catalyst, electrolyte or flow cell is a bounded prior, not validation of this implementation.

The optional numerical checks establish agreement with declared equations. They do not establish experimental calibration, industrial suitability or engineering acceptance.
