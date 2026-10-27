# Repository instructions

- The functional skill is `skills/ai-comsol-skill`. Keep installation location independent of the user's selected calculation workspace.
- Preserve evidence boundaries: literature, synthetic examples, explicit assumptions, numerical verification and experimental validation are separate. Do not claim generic CFD or calibrated catalyst performance from the planar surrogate.
- Keep AI process assets, downloaded sources, generated models, temporary files and logs under the active project's `AGENTCODING`. Do not commit them. Functional source, assets, references and repository checks remain in their normal repository locations.
- Do not overwrite an existing output or source model. Native runners require an explicit installation or `COMSOL_ROOT`/`DWSIM_ROOT`; no simulator binaries are bundled.
- On Windows use `py -3.12 -B tools/validate_repo.py` before committing. This checks text/JSON integrity and runs behavioral tests without requiring COMSOL. Physics changes additionally require native solver, conservation and readback checks in a licensed environment.
- Keep UTF-8 text and relative documentation links. Add third-party source metadata with access, license, evidence grade, adoption decision and applicability boundary.
