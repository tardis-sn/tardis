.. _nlte:

NLTE treatment
--------------

NLTE treatment of lines is available both in ~LTEPlasma and the ~NebularPlasma class. This can be enabled by specifying
which species should be treated as NLTE with a simple list of tuples (e.g. [(20,1)] for Ca II).

First let's dive into the basics:

There are two rates to consider from a given level.

.. math::

    R_{\textrm{upper}\rightarrow\textrm{lower}} &= \underbrace{A_{ul} n_u}_\textrm{spontaneous emission}
            + \underbrace{B_{ul} n_u \bar{J}_\nu}_\textrm{stimulated emission} +
            \underbrace{C_{ul} n_u n_e}_\textrm{collisional deexcitation}\\
            &= n_u \underbrace{(A_{ul} + B_{ul}\bar{J}_\nu + C_{ul} n_e)}_{r_{ul}} \\

    R_{\textrm{lower}\rightarrow\textrm{upper}} &= \underbrace{B_{lu} n_l \bar{J}_\nu}_\textrm{stimulated absorption} +
                \underbrace{C_lu\,n_l\,n_e}_\textrm{collisional excitation}\\
                &= n_l \underbrace{(B_{lu}\bar{J}_\nu + C_{ul}\ne)}_{r_{lu}},

where :math:`\bar{J}_\nu` (in LTE this is :math:`B(\nu, T)`) denotes the mean intensity at the frequency of the line and
:math:`n_e` the number density of electrons.

Next, we calculate the rate of change of a level by adding up all outgoing and all incoming transitions from level :math:`j`.


.. math::

    \frac{dn_j}{dt} = \underbrace{\sum_{i \ne j} R_{ij}}_\textrm{incoming rate} -
                        \underbrace{\sum_{i \ne j} R_{ji}}_\textrm{outgoing rate}

In a statistical equilibrium all incoming rates and outgoing rates add up to 0 (:math:`\frac{dn_j}{dt}=0`). We use this to
calculate the level populations using the rate coefficients (:math:`r_ij, r_ji`).

For a three level atom we have:

.. math::

    \frac{dn_1}{dt} &= \underbrace{n_2 r_{21} + n_3 r_{31}}_\textrm{incoming rate}
                    - \underbrace{(n_1 r_{12} + n_1 r_{13})}_\textrm{outgoing rate} = 0\\

    \frac{dn_2}{dt} &= \underbrace{n_1 r_{12} + n_3 r_{32}}_\textrm{incoming rate}
                    - \underbrace{(n_2 r_{21} + n_2 r_{23})}_{outgoing rate} = 0\\

    \frac{dn_3}{dt} &= \underbrace{n_1 r_{13} + n_2 r_{23}}_\textrm{incoming rate}
                    - \underbrace{(n_3 r_{32} + n_3 r_{31})}_\textrm{outgoing rate} = 0,


which can be written in matrix from:

.. math::

    \left(\begin{matrix}
        -(r_{12} + r_{13}) & r_{21} & r_{31}\\
        r_{12} & -(r_{21} + r_{23}) & r_{32}\\
        r_{13} & r_{23} & -(r_{31} + r_{32}) \\
    \end{matrix}\right)
    \left(
    \begin{matrix}
        n_1\\
        n_2\\
        n_3\\
    \end{matrix}
    \right)
    =
    \left(
    \begin{matrix}
        0\\
        0\\
        0\\
    \end{matrix}
    \right)

To solve for the level populations we need an additional constraint: :math:`n_1 + n_2 + n_3 = N`. By seting :math:`N = 1`:
we can get the relative rates:

.. math::

    \left(\begin{matrix}
        1 & 1 & 1\\
        r_{12} & -(r_{21} + r_{23}) & r_{32}\\
        r_{13} & r_{23} & -(r_{31} + r_{32}) \\
    \end{matrix}\right)
    \left(
    \begin{matrix}
        n_1\\
        n_2\\
        n_3\\
    \end{matrix}
    \right)
    =
    \left(
    \begin{matrix}
        1\\
        0\\
        0\\
    \end{matrix}
    \right)


Now we go back and look at the rate coefficients used for a level population - as an example :math:`\frac{dn_2}{dt}`:

.. math::

    \frac{dn_2}{dt} &= n_1 r_{12} - n_2 (r_{21} + r_{23}) + n_3 r_{32}\\
                &= n_1 B_{12} \bar{J}_{12} + n_1 C_{12} n_e - n_2 A_{21} - n_2 B_{21} \bar{J}_{21} - n_2 C_{21} n_e
                        - n_2 B_{23} \bar{J}_{23} - n_2 C_{23} n_e + n_3 A_{32} + n_3 B_{32} \bar{J}_{32} + n_3 C_{32} n_e,


&=
                         + n_3 A_{32}  + n_3 C_{32} n_e,

Next we will group the stimulated emission and stimulated absorption terms as we can assume :math:`\bar{J_{12}} = \bar{J_{21}}`:

.. math::

    \frac{dn_2}{dt} &= n_1 (B_{12} \bar{J}_{12} \underbrace{(1 - \frac{n_2}{n_1}\frac{B_{21}}{B_{12}})}_\textrm{stimulated emission term} + C_{12} n_e)
                    - n_2 (A_{21} + C_{23} n_e + n_2 B_{23} \bar{J}_{23} \underbrace{(1 - \frac{n_3}{n_2}\frac{B_{32}}{B_{23}})}_\textrm{stimulated emission term})
                    + n_3 (A_{32} + C_{32} n_e)





