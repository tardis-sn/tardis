#!/bin/bash
set -euo pipefail

# Set a few environment variables that are not set due to
# https://github.com/conda/conda-build/issues/3993
PIP_NO_BUILD_ISOLATION=False \
    PIP_NO_DEPENDENCIES=True \
    PIP_IGNORE_INSTALLED=True \
    PIP_NO_INDEX=True \
    PYTHONDONTWRITEBYTECODE=True \
    ${PYTHON} -m pip install . --no-deps -vv
