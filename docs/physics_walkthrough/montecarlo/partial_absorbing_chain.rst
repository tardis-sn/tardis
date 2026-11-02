Element-selected absorbing chains
=================================

Type IIP transport can skip internal macro-atom jumps while preserving the
probabilities of eventual emission and cooling. In partial mode, let S contain
all represented real atomic levels of the requested elements. Only internal
edges with both endpoints in S are eliminated. The shared ionization (i) and
thermal (k) pools stay outside S, even when every atomic element is selected.

For each cell, Q contains the eliminated transition probabilities, L contains
the original probabilities with just those edges removed, and r is the sum of
each selected state's remaining outgoing channels. The first-exit matrix B and
conditional channel probabilities C satisfy

.. math::

    (I-Q)B = \operatorname{diag}(r), \qquad
    C_{j,c} = L_{j,c}/r_j, \qquad
    (I-Q)^{-1} L = BC.

Rows of B are activation states; columns identify exit-source states. This
matrix need not be symmetric. A draw from B selects the global source of the
next retained channel, followed by an independent draw from C. A retained
internal transition continues the interaction; only a negative transition code
ends it. Each return to S starts a fresh, independent pair of draws. Outside S,
transport samples the original channel probabilities directly.

Emission-only states and global level zero are included. Explicit int64 forward
and inverse maps preserve global identities, including gaps in the level
indices. A zero residual row stays zero and cannot supply an exit. Metadata
and numerical probabilities are rebuilt consistently on every opacity solve;
continuum activation indices and estimator dimensions keep their original
ordering and shell count.

Sparse LU solves are checked for finite, nonnegative, normalized probabilities
and a scale-aware equation residual. A reachable closed class without a physical
exit raises a diagnostic containing the cell and states. Channels are never
removed by a numerical probability threshold. The sampler allows a CDF endpoint
deficit of at most 1e-10 due to floating-point summation; larger deficits or
activations without outgoing probability raise an error.

Exact small-network tests compare eventual channel distributions with an
independent dense solve of the unreduced network. Compiled sampling tests check
the analytic 50 percent / 25 percent / 25 percent re-entry example with binomial
uncertainty bounds. Full mode preserves the existing default and its two-draw
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
cutoff scenario). Full mode preserves the legacy aggregate elimination of these
signed internal coefficients and their legacy signed matrix. That legacy
case does not carry the probability nonnegativity guarantee of the new modes.
Partial and explicit modes reject negative original probabilities with a cell
and transition diagnostic: a signed channel cannot be sampled as an explicit
probability. Select full mode for those states. The reduction does not clip
coefficients or change the underlying plasma/opacity physics.

Exceptionally near-trapping partial chains can lose tiny exits through
floating-point cancellation. Their diagonals use actual exits and off-diagonal
jumps; if float64 LU violates probability conservation, a Decimal
elimination with at least 80 digits (40 guard digits beyond the smallest
positive exit) retries the conservation-form system. This preserves tiny physical
channels without inventing probabilities or removing transitions. Full mode
retains the original numerical calculation for regression compatibility,
including its cancellation limits for exceptionally near-trapping chains.
