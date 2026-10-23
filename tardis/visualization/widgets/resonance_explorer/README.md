# Single-shell resonance explorer

This Panel widget explores the resonance equation for one straight packet
flight through one shell with `v(r) = m r + q`. It lives alongside TARDIS's
other visualization widgets and requires no atomic data or network service.
The scientific assumptions and lessons are described below. A future
nonhomology documentation section can import this widget without hosting a
live application in the published docs.

## Use in a notebook

Like `ShellInfoWidget` and `LineInfoWidget`, construct the widget and display
its Panel layout in a local Jupyter notebook:

```python
from tardis.visualization import ResonanceExplorerWidget

explorer = ResonanceExplorerWidget()
explorer.display()
```

The published documentation cannot run the Panel controls. The notebook must
use a Python environment with this TARDIS checkout and the visualization
dependencies installed.

## Run as a local server

From the checkout root, run the included script:

```shell
./tardis/visualization/widgets/resonance_explorer/run.sh
```

Open `http://127.0.0.1:5118/` in Firefox. Keep the terminal open while using
the explorer; press **Ctrl+C** there to stop it. The script uses the
`tardis-merger-workflow` environment and imports this TARDIS checkout,
regardless of the directory from which it is launched. Set
`RESONANCE_EXPLORER_PORT=5006` before the command to use another port. If the
environment is installed elsewhere, set `TARDIS_MERGER_WORKFLOW_ENV` to its
absolute prefix. The script calls that environment's Python directly, so it
does not need the `conda` command on `PATH`. The server binds only to
`127.0.0.1` and serves its JavaScript and CSS inline for offline use.

## Guided lessons

The initial screen starts with a controlled comparison of homologous expansion
and linked nonhomology. Switching lessons 1 and 2 changes only the radial-law
intercept; the shell, packet ray, line target, plot limits, and shell color
scale stay fixed. The same scales remain in place when the comparison toggle
shows or hides the other velocity law. The main plot shows the
active forward segment of projected velocity, the line target, numbered
verified crossings, and gradient reversals. Smaller views show the shell ray,
the radial law and its edge velocities, and the projected gradient. Between the
shell view and radial-law plot, a compact shell-free construction uses the same
spatial limits as the shell view. It shows each distinct verified resonance on
the full packet trajectory. Radial arrows use the local shell-field color and
a thin black outline; red arrows show their signed projection along packet travel,
`v_parallel = v_r x/r = w`. All arrows share one display scale, and red
projection arrows lie directly on the thin gray packet trajectory. On a radial
ray, radial and projected vectors coincide; the radial stroke is slightly
wider so it remains visible beneath the red one. Gray start and end ticks sit
behind the vectors; the shell view carries the prominent trajectory markers.
The shell view has no coordinate axes, and the vector construction has no
vertical axis; the data curves retain their labeled axes without box frames.

The seven lesson cards cover homology, three linked crossings, launch clipping,
target miss and extremum contact, whole-shell collapse, homologous collapse,
and special solver or event branches. Paired lesson states expose their one
changed input. Below the plot, four production-ordered slot chips describe
each algebraic result and highlight its coordinate. The decision path keeps
production acceptance distinct from independent verification of the original
unsquared equation. A displayed potential resonance is not an absorbed
packet.

## Inputs

All lengths and velocities use consistent **illustrative units**. Each lesson
state sets the shell edges and their signed radial velocities, impact parameter `b`, signed
ray launch coordinate `x_0`, and projected line target `w`. The app derives
`r_0 = hypot(x_0,b)`, `mu_0 = x_0/r_0`, `m`, and `q`. Distance `s` grows
forward with `x = x_0 + s`, including when the ray moves inward in radius.
The displayed intercept uses the packet transport arithmetic
`q = v_out - m*r_out`, with `m = (v_out-v_in)/(r_out-r_in)`. Algebraically the
inner-edge expression is the same, but its floating-point rounding can differ
and change an exact `Q=0` solver route.
The target corresponds to `w = c(1 - nu_line/nu_packet)`; frequencies are not
input controls. An optional assumed electron scatter truncates the active
flight at `s_stop`.

## Interpretation

Displayed resonance crossings are **potential** line encounters. The app
applies shell, active-flight, earlier-stop, and boundary ownership rules but
does not calculate opacity, optical depth, or a stochastic packet absorption.
Packet transport treats a contact as singular for Sobolev depth only when the
projected gradient is exactly `G == 0.0`; a small finite `G` remains finite.

## Check

From the checkout root, run focused checks in the same environment:

```shell
PYTHONPATH=. /home/connor/micromamba/envs/tardis-merger-workflow/bin/python -m pytest -q \
  tardis/visualization/widgets/resonance_explorer/tests
/home/connor/micromamba/envs/tardis-merger-workflow/bin/ruff check \
  tardis/visualization/widgets/resonance_explorer
```

The `PYTHONPATH` and `python -m pytest` form selects this TARDIS checkout even
if the environment has another editable checkout installed.
