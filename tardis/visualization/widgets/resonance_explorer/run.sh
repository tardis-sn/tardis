#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "$script_dir/../../../.." && pwd)"
env_prefix="${TARDIS_MERGER_WORKFLOW_ENV:-${HOME}/micromamba/envs/tardis-merger-workflow}"
python_bin="$env_prefix/bin/python"
port="${RESONANCE_EXPLORER_PORT:-5118}"

if [[ ! -x "$python_bin" ]]; then
    printf 'Python not found in tardis-merger-workflow environment: %s\n' "$python_bin" >&2
    exit 1
fi

cd "$repo_root"
export PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/resonance-explorer-matplotlib-${UID}}"
mkdir -p "$MPLCONFIGDIR"

exec "$python_bin" -m tardis.visualization.widgets.resonance_explorer \
    --no-browser --port "$port" "$@"
