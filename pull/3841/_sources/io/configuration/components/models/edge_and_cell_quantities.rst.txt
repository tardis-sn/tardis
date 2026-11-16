.. _model-edge-cell-quantities:

************************************
Edge-Defined and Cell-Defined Values
************************************

TARDIS discretizes the ejecta with a one-dimensional finite-volume grid. The
grid consists of :math:`N` spherical shells (cells), and each model quantity is
either *edge-defined* or *cell-defined*:

- An **edge-defined** quantity is a value at a shell boundary. For :math:`N`
  shells, TARDIS stores :math:`N + 1` edge values, and shell :math:`i` spans
  edges :math:`i` and :math:`i + 1`.
- A **cell-defined** quantity is a single value that represents the whole
  shell. TARDIS stores :math:`N` cell values and treats each value as uniform
  throughout its shell. A cell value is not attached to a point inside the
  shell, including the shell midpoint.

Hydrodynamics and radiative-transfer codes do not all use these conventions.
Some codes tabulate values at grid points (nodes), at cell centres, or on a
staggered grid. Each TARDIS importer therefore maps its source data onto TARDIS
edges and cells. This page documents those mappings.

.. contents::
    :local:

Quantity Conventions
====================

.. list-table::
    :header-rows: 1
    :widths: 30 15 55

    * - Quantity
      - Defined on
      - Notes
    * - Velocity (``velocity``, ``v_inner``, ``v_outer``)
      - Edges
      - ``v_inner`` and ``v_outer`` are views of the same :math:`N + 1` edge
        values.
    * - Radius (``radius``, ``r_inner``, ``r_outer``)
      - Edges
      - Homologous models derive the radius from :math:`r = v t_\mathrm{exp}`.
        Nonhomologous models store radii independently of velocities.
    * - Velocity inside a shell
      - Continuous
      - Homologous models use :math:`v = r / t_\mathrm{exp}`. Nonhomologous
        models interpolate linearly in radius between the two edge velocities.
    * - Mass density
      - Cells
      - The mean density of the shell.
    * - Elemental and isotopic mass fractions
      - Cells
      -
    * - Electron density and temperature read from a file
      - Cells
      -
    * - Radiative temperature :math:`T_\mathrm{rad}` and dilution factor
        :math:`W`
      - Cells
      - Supplied values and Monte Carlo estimator updates are both shell
        values.
    * - Plasma state and line opacities
      - Cells
      - Computed from the cell-defined quantities above.
    * - Volume
      - Cells
      - :math:`\frac{4}{3}\pi(r_\mathrm{outer}^3 - r_\mathrm{inner}^3)`,
        computed from the edges.

The ``v_middle`` and ``r_middle`` properties are the arithmetic means of the
two edge values of each shell. They are not volume-weighted or mass-weighted
centroids.

.. _model-leading-row:

Tabular Inputs and the Leading Row
==================================

Tabular input formats store edge-defined and cell-defined columns in the same
table, so each table has the same number of rows for both. TARDIS reads such a
table with :math:`N + 1` rows as follows:

- Velocity (and radius, if present) in row :math:`j` is an edge. Row 0 holds
  the inner edge of the innermost shell, and row :math:`j \ge 1` holds the
  outer edge of shell :math:`j - 1`.
- Cell-defined values in row :math:`j \ge 1` belong to the shell between the
  edges in rows :math:`j - 1` and :math:`j`.
- Cell-defined values in row 0 nominally describe the region inside the
  innermost edge. TARDIS discards them, so they can contain placeholder
  values.

For example, consider a table with the following rows:

.. code-block:: text

    velocity  density
    9000      5.0e-10   <- inner edge; density discarded
    10500     2.0e-10   <- shell 0 spans 9000-10500 km/s
    12000     9.0e-11   <- shell 1 spans 10500-12000 km/s

TARDIS applies this convention to the
:ref:`simple ASCII density format <densitycust>`, the
:ref:`ASCII and CSV abundance formats <abundancecust>`, the CSV part of a
:ref:`CSVY model <csvy-model>`, and the ARTIS and CMFGEN file types.

.. _model-generated-cell-values:

Generated Cell Values
=====================

When TARDIS generates cell values instead of reading them, it samples a
point in each shell:

- Built-in density profiles (``branch85_w7``, ``power_law``, ``exponential``)
  are evaluated at ``v_middle``. The result is a point sample, not a volume
  average over the shell. Therefore, for a steep profile, the shell mass
  :math:`\rho V` differs from the integral of the profile over the shell.
- When :math:`T_\mathrm{rad}` is not supplied and ``initial_t_rad`` is not
  set, the initial radiative temperature is evaluated at ``v_middle``.
- When :math:`W` is not supplied, the geometric dilution factor is evaluated
  at ``r_middle``.
- Uniform abundances assign the same mass fractions to every shell.

Boundaries Inside a Shell
=========================

``v_inner_boundary`` and ``v_outer_boundary`` can fall inside a shell of the
input grid. In that case, TARDIS moves the affected edge to the requested
boundary and keeps the cell-defined values of that shell unchanged. The
truncated shell keeps the density of the full shell, so its mass decreases in
proportion to its volume.

If a requested boundary lies outside the input grid, TARDIS issues a warning
and extends the innermost or outermost shell to that boundary. The extended
shell keeps its cell-defined values, so its mass increases.

Importer Mappings
=================

The following table summarizes how each importer maps its source data onto
TARDIS edges and cells.

.. list-table::
    :header-rows: 1
    :widths: 22 35 43

    * - Importer
      - Source convention
      - TARDIS mapping
    * - ``model.structure`` ``type: specific``
      - Velocity range and number of shells.
      - TARDIS spaces :math:`N + 1` edges linearly between ``start`` and
        ``stop``. Densities are
        :ref:`generated <model-generated-cell-values>` at ``v_middle``.
    * - ``simple_ascii`` density and ASCII or CSV abundance files
      - Velocity is the outer edge of each cell. Density and mass fractions
        are cell means.
      - TARDIS applies the leading-row convention directly.
    * - CSVY (CSV part)
      - Velocity and radius are outer edges. All other columns are cell
        values.
      - TARDIS applies the leading-row convention directly. This includes
        ``t_rad`` and ``dilution_factor``.
    * - CSVY (YAML ``velocity`` section)
      - Velocity range and number of shells.
      - TARDIS spaces :math:`N + 1` edges linearly, as for
        ``type: specific``.
    * - ARTIS (``artis`` file type, :func:`~tardis.io.model.artis.readers.read_artis_model`)
      - Each row gives the outer velocity of an ARTIS cell and the cell-mean
        density and mass fractions. The innermost ARTIS cell extends to
        :math:`v = 0`.
      - TARDIS applies the leading-row convention directly. The first ARTIS
        velocity becomes the inner edge of the TARDIS grid, and TARDIS drops
        the innermost ARTIS cell. No reinterpretation is needed.
    * - CMFGEN (``cmfgen_model`` file type and ``cmfgen2tardis``)
      - CMFGEN tabulates quantities at radial grid points.
      - Each grid-point velocity becomes the outer edge of a TARDIS shell.
        The density, electron density, temperature, and mass fractions at that
        grid point become the cell values of the shell inside it. The values at
        the innermost grid point are discarded. This is a piecewise-constant
        reconstruction from the outer grid point of each shell, not an average
        over the shell.
    * - Blondin toy model (:func:`~tardis.io.model.readers.blondin_toymodel.convert_blondin_toymodel`)
      - Velocity is given at cell centres.
      - The converter replaces each velocity with the midpoint between that
        cell centre and the next cell centre, and extrapolates the outermost
        edge linearly. The other cell values stay in their rows. Through the
        leading-row convention, the innermost toy-model cell is dropped.
    * - Arepo (:func:`~tardis.io.model.arepo.utils.export_profile_to_csvy`)
      - Arepo values are cell values on a Voronoi mesh. The profile functions
        return one row per Arepo cell along a line of sight, and
        :func:`~tardis.io.model.arepo.utils.rebin_profile` optionally reduces
        them to mass-weighted bin means.
      - The exporter does not rebin. It selects ``nshells`` evenly spaced rows
        and writes the velocity of each row (for rebinned data, the
        mass-weighted mean velocity of the bin) to the CSVY
        ``velocity`` column, which TARDIS reads as an outer edge. As a result,
        the TARDIS shell edges coincide with representative cell velocities,
        not with bin boundaries.
    * - STELLA (:func:`~tardis.io.model.readers.stella.read_stella_model`)
      - STELLA labels its columns as cell-centre values (``cell center m``,
        ``cell center R``, ``cell center v``), outer-edge values
        (``outer edge m``, ``outer edge r``), or cell averages.
      - The reader returns the STELLA table unchanged and does not build a
        TARDIS model. When you convert it, build the edges from the outer-edge
        columns, not from ``cell center v``.
    * - SNEC (:mod:`tardis.io.model.snec`)
      - SNEC is a Lagrangian hydrodynamics code. Its outputs are indexed by
        ``cell_id``.
      - The readers return the SNEC data unchanged and do not build a TARDIS
        model. A converter must decide which SNEC variables are edge-defined.

Writing a New Importer
======================

A new importer or converter should do the following:

1. Document whether the source code defines each variable at edges, at cell
   centres, or as cell averages.
2. Produce :math:`N + 1` velocity edges (and radius edges for nonhomologous
   models) and :math:`N` cell values. When writing a tabular TARDIS format,
   include a leading row whose cell values are placeholders.
3. Convert point values to cell values explicitly, for example by averaging
   over the shell volume or mass, and document the method. If a point value is
   used directly as a cell value, document which point it comes from.
4. Do not write cell-centred velocities into an edge column without
   documenting this reinterpretation.

The finite-volume discussion in the
`CCA summer school notes on advection <https://zingale.github.io/cca-summer-school/advection/advection-finite-volumes.html>`_
explains the difference between cell averages and point values.
