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
is removed.

Emission-only states and global level zero are included. Explicit int64 forward
and inverse maps preserve global identities, including gaps in the level
indices. A zero residual row stays zero and cannot supply an exit. Metadata
and numerical probabilities are rebuilt consistently on every opacity solve;
continuum activation indices and estimator dimensions keep their original
ordering and shell count.

Sparse LU solves are checked for finite, normalized rows and a scale-aware
equation residual. Input coefficients and reduced entries are not rejected
based on their sign in any reduction mode. A reachable closed class without a
physical exit raises a diagnostic containing the cell and states. Channels are
never removed by a numerical probability threshold. The sampler allows a CDF endpoint
deficit of at most 1e-10 due to floating-point summation; larger deficits or
activations without outgoing probability raise an error.

Exact small-network tests compare eventual channel distributions with an
independent dense solve of the unreduced network. Compiled sampling tests check
an analytic 50 percent / 18.75 percent / 31.25 percent emission distribution with
direct thermal emission, pool activation, and re-entry through an unselected
element. Every pool-option combination is checked with binomial uncertainty
bounds. Full mode preserves the existing default and its two-draw
sampling behavior. For nonnegative kernels, the different modes preserve
distributions, rather than individual random packet histories.