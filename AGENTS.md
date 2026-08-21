# AGENTS.md

## Repository purpose

This repository is a mirror of the actual MTPLite codebase. The source of truth
is the external/original repository, not this mirror.

## Working rules

- Treat changes here as mirror maintenance, documentation, inspection, or
  verification work unless the user explicitly says otherwise.
- Do not assume that a local code change should be treated as the canonical
  implementation. If code is changed here, clearly report that it must be
  propagated back to the actual codebase.
- Preserve the existing layout and filenames so the mirror remains easy to
  compare with the source repository.
- Avoid broad refactors, dependency upgrades, generated-file churn, or cleanup
  unrelated to the requested task.
- Keep user-created or untracked files intact unless the user explicitly asks
  for them to be changed or removed.

## Repository layout

- `mtp_lite_v1.1/` — mirrored implementation and scripts.
- `docs/` — supporting documentation and run-order notes.
- `report.pdf` — project report/reference material.

## v1.1 configuration conventions

- The v1.1 entry points use visible, user-editable placeholders; do not
  reintroduce a run-directory resolver, environment-variable path model, or
  positional path arguments unless explicitly requested.
- One project artifact root must be used consistently: `DIR` in
  `extractUnikmers.sh`, Python-stage `base_dir` values, and
  `barcode_assembler.py`'s `DEFAULT_BASE_DIR` should identify the same root.
  `barcode_assembler.py` retains `--base-dir` as an optional override.
- `final_read_selection.py` must use `<project-root>/output` as `output_dir`,
  because it reads `anchors.pkl` from there. Its selected-read FASTA is
  `<project-root>/output/<prefix>.fasta`; configure `assembly.sh`'s
  `READ_FILE` to that exact path, retain the same `PREFIX`, and choose its
  `OUT_DIR` independently for assembly results.
- When changing path configuration, update the root README and relevant files
  under `docs/` in the same maintenance change.

## Verification

When making a change, run the narrowest relevant checks available and state
clearly what was verified. If no automated test suite exists for the affected
area, say so rather than inventing one.
