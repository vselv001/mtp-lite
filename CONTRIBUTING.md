# Contributing to MTPLite

Thank you for helping improve the MTPLite research record. This repository is a
mirror of the original implementation, so changes here are not canonical unless
they are also propagated to the source repository.

## Appropriate contributions

Useful changes include:

- corrections or clarifications to documentation and reproducibility records;
- narrowly scoped fixes with a clear scientific or operational rationale;
- additional validation, benchmark, or provenance material that identifies its
  inputs, parameters, environment, and limitations;
- tests or checks that do not require committing sequencing data or large
  generated artifacts.

Avoid broad refactors, dependency upgrades, formatting-only churn, or generated
file updates unrelated to a specific research or maintenance goal.

## Before opening a change

1. Describe the scientific, reproducibility, or maintenance motivation.
2. Keep raw reads, references, binary indices, logs, and generated assemblies
   out of version control; use checksums and provenance instead.
3. Preserve the current layout and filenames so the mirror remains comparable
   with its source.
4. State whether the change must also be applied to the canonical repository.
5. Run the narrowest relevant validation and report the result.

For documentation changes, check internal links, command examples, and claims
against the checked-in implementation. For code changes, distinguish clearly
between a reproduced result and a new experiment with modified inputs or
parameters.

## Reporting results

When sharing a new MTPLite result, include the information in
[docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md): repository commit, input
and reference provenance/checksums, exact configuration, resolved environment,
hardware context, logs, and evaluation command. Do not generalize the v1.1
reference result to another dataset without supporting evidence.

## License

By contributing, you agree that your contribution may be distributed under the
[MIT License](LICENSE).
