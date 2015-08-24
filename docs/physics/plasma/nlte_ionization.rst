***************
NLTE Ionization
***************

Equation to solve (He used as example - 3 ionisation states):

.. math::
    \left( \begin{array}{cccc}
    -r_{12} & r_{12} & 0 & \\
    r_{12} & -(r_{21}+r_{23}) & r_{32} \\
    0 & r_{23} & -r_{32} \end{array} \right)
    \left( \begin{array}{c}
    N_{1} \\
    N_{2} \\
    N_{3} \end{array} \right) =
    \left( \begin{array}{c}
    0 \\
    0 \\
    0 \end{array} \right)

For now only transitions between adjacent ionisation states are considered, i.e. :math:`r_{13}` etc. = 0.

Three types of transition are currently considered, where `i` is the lower ionization state and `k` the upper.

Spontaneous Recombination:

.. math::
    \alpha_{i}^{sp.} = 4\pi\Phi_{ik}\int_{\nu_{i}}^{\infty}\frac{a_{ik(\nu)}}{h\nu}\frac{2h\nu^{3}}{c^2}\exp(-\frac{h\nu}{kT_{e}})

Stimulated Recombination:

.. math::
    \alpha_{i}^{st.} = 4\pi\Phi_{ik}\int_{\nu_{i}}^{\infty}\frac{a_{ik(\nu)}}{h\nu}J_{\nu}\exp(-\frac{h\nu}{kT_{e}})

Photoionization:

.. math::
    \gamma_{i} = 4\pi\int_{\nu_{i}}^{\infty}\frac{a_{ik(\nu)}}{h\nu}J_{\nu}\exp(-\frac{h\nu}{kT_{e}})

:math:`\Phi_{ik}` is the ion population ratio in LTE:

.. math::
    \Phi_{ik} = \left(\frac{n_{i}}{n_{k}n_{e}}\right)_{LTE}

(This is the inverse of the `PhiSahaLTE` property output `phi` in Tardis.)
