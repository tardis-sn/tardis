Element-selected absorbing chains
=================================

Type IIP transport can skip internal macro-atom jumps while preserving the
probabilities of eventual emission and cooling. In partial mode, let S contain
all represented real atomic levels of the requested elements. Only internal
edges with both endpoints in S are eliminated. For example, this would be
jumping over internal transitions that happen in a specified element.
The shared ionization (i) and thermal (k) pools stay outside S by default.
The independent ``include_ionization_pool`` and ``include_thermal_pool`` options
add their respective pool states to S in partial mode. Either pool or both can
be included without selecting additional atomic elements.

For each cell, Q contains the eliminated transition probabilities, L contains
the original probabilities with just those edges removed, and r is the sum of
each selected state's remaining outgoing channels. The first-exit matrix B and
conditional channel probabilities C satisfy

.. math::

    (I-Q)B = \operatorname{diag}(r), \qquad
    C_{j,c} = L_{j,c}/r_j, \qquad
    (I-Q)^{-1} L = BC.

Rows of B are activation states; columns identify exit-source states.
A draw from B selects the global source of the
next retained channel, followed by an independent draw from C. A retained
internal transition continues the interaction; only a negative transition code
ends it. Each return to S starts a fresh, independent pair of draws. Outside S,
transport samples the original channel probabilities directly.

With k included, selected-atom to k to selected-atom paths are collapsed inside
the matrix. With i included, the analogous ionization/recombination paths are
collapsed. If both are included, internal transfers between the pools also enter
Q. A pool can be the exit-source state selected by B: emission and cooling from
it are retained channels, as are jumps to unselected atoms or an excluded pool.
Those exits are sampled explicitly. Direct activation at an included pool uses
its reduced row; activation at an excluded pool uses its original channels.

Shared pools couple the selected elements into one reduction. The final emission
distribution is preserved, but the explicit history of visits to included pools
is removed. Diagnostics that count individual pool visits would require separate
accounting for the collapsed transitions.

Emission-only states and global level zero are included. Explicit int64 forward
and inverse maps preserve global identities, including gaps in the level
indices. A zero residual row stays zero and cannot supply an exit. Metadata
and numerical probabilities are rebuilt consistently on every opacity solve;
continuum activation indices and estimator dimensions keep their original
ordering and shell count.

Sparse LU solves are checked for finite, normalized rows and a scale-aware
equation residual. Nonnegativity checks apply to networks with nonnegative
original coefficients in every reduction mode. Signed networks retain the
legacy behavior described below. A reachable closed class without a physical
exit raises a diagnostic containing the cell and states. Channels are never
removed by a numerical probability threshold. The sampler allows a CDF endpoint
deficit of at most 1e-10 due to floating-point summation; larger deficits or
activations without outgoing probability raise an error.

Exact small-network tests compare eventual channel distributions with an
independent dense solve of the unreduced network. Compiled sampling tests check
an analytic 50 percent / 18.75 percent / 31.25 percent emission distribution with
direct thermal emission, pool activation, and re-entry through an unselected
element. Every pool-option combination is checked with binomial uncertainty
bounds. Full mode preserves the existing default and its two-draw
sampling behavior. The different modes preserve distributions, rather than
individual random packet histories.

A smaller matrix does not guarantee faster transport: explicit cycling through
unselected states and shared pools can dominate. Compare setup and solve cost,
selected-state count, dense matrix bytes, explicit steps, and total transport
cost with the existing ASV remote benchmarking workflow on the same network and
packet count. No hardware-independent speedup is implied.

Signed-channel limitation
-------------------------

Some accepted plasma states produce negative stimulated-recombination-corrected
photoionization coefficients (for example, the existing 1500 K outer-shell
cutoff scenario). Full, partial, and explicit modes preserve signed coefficients
without clipping or sign-based rejection. Partial reduction applies the same
algebraic elimination to the selected block; negative retained channels or
negative reduced entries do not introduce additional termination checks.
Including the ionization or thermal pool remains optional for signed networks.

For example, a direct selected-state coefficient of -0.1 can combine with a
two-step contribution of +0.3 to give an exit coefficient of +0.2. Other
selections can leave signed entries. The existing cumulative sampler is retained
in all modes for compatibility. Running a signed case successfully does not
establish a valid categorical sampling law or mode-invariant sampled spectra.
This limitation is shared by full and partial reduction; resolving it is a
separate physics change. General physical treatment requires consistently separating
stimulated recombination from absorption, including corresponding changes to
packet opacities and emission, rather than only removing a sign check. See
`Lucy (2002), section 4.1 <https://arxiv.org/html/astro-ph/0107377v3#S4.SS1>`_.

Exceptionally near-trapping partial chains can lose tiny exits through
floating-point cancellation. Their diagonals use actual exits and off-diagonal
jumps; if float64 LU violates probability conservation, a Decimal
elimination with at least 80 digits (40 guard digits beyond the smallest
positive exit) retries the conservation-form system. This preserves tiny physical
channels without inventing probabilities or removing transitions. Full mode
retains the original numerical calculation for regression compatibility,
including its cancellation limits for exceptionally near-trapping chains.
