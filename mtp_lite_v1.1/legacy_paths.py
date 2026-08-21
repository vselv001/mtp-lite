"""Path helpers for isolated MTPLite v1.1 legacy runs.

The original v1.1 scripts wrote their artifacts into the repository root.
Keeping a run root separate lets legacy results coexist with the optimized
pipeline, while reads and references remain shared, read-only inputs.
"""

import os
from pathlib import Path


PIPELINE_DIR = Path(__file__).resolve().parent
REPOSITORY_DIR = PIPELINE_DIR.parent
DEFAULT_RUN_DIR = PIPELINE_DIR / "runs" / "active"
DEFAULT_READS_PATH = REPOSITORY_DIR / "input" / "hifi_1000x_chr1.fastq"


def _resolve_path(path):
    return os.path.abspath(os.path.expanduser(os.fspath(path)))


def resolve_run_dir(run_dir=None):
    """Return the artifact root for a legacy run.

    ``MTP_LITE_V1_RUN_DIR`` makes the existing stage commands reusable for a
    specific run without placing generated files in the repository root.
    """
    return _resolve_path(
        run_dir or os.environ.get("MTP_LITE_V1_RUN_DIR") or DEFAULT_RUN_DIR
    )


def resolve_reads_path(reads_path=None):
    """Return the shared reads input, optionally overridden for a rerun."""
    return _resolve_path(
        reads_path or os.environ.get("MTP_LITE_V1_READS") or DEFAULT_READS_PATH
    )
