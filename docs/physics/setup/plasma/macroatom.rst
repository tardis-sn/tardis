.. _macroatom:

Macro Atom
----------

The macro atom is described in detail in :cite:`Lucy2002`. The basic principal is that when an energy packet
is absorbed that the macro atom is on a certain level. Three probabilities govern the next step the P\ :sub:`up`,
P\ :sub:`down` and P\ :sub:`down emission` being the probability for going to a higher level, a lower level and a lower
level and emitting a photon while doing this respectively (see Figure 1 in :cite:`Lucy2002` ).

The macro atom is the most complex idea to implement as a data structure. The setup is done in `~tardisatomic`, but
we will nonetheless discuss it here (as `~tardisatomic` is even less documented than this one).

For each level, we look at the line list to see what transitions (upwards or downwards are possible). We create a two arrays,
the first is a long one-dimensional array containing the probabilities. Each level contains a set of probabilities. The first
part of each set contains the upwards probabilities (internal upward), the second part the downwards probabilities
(internal downward), and the last part is the downward and emission probability.


each set is stacked after the other one to make one long one dimensional `~numpy.ndarray`.

The second array is for book-keeping; it has exactly the length as levels (with an example for the Si II level 15):

+--------+------------------+------------+----------------+-----------------+
|Level ID| Probability index|N\ :sub:`up`| N\ :sub:`down` | N\ :sub:`total` |
+========+==================+============+================+=================+
|14001015| ???              |17          | 5              | 17 + 5*2 = 27   |
+--------+------------------+------------+----------------+-----------------+


We now will calculate the transition probabilites, using the radiative rates in Equation 20, 21, and 22
in :cite:`Lucy2002`. Then we calculate the downward emission probability from Equation 5, the downward and
upward internal transition probabilities in :cite:`Lucy2003`.

.. math::
    p_\textrm{emission down}&= {\cal R}_{\textrm{i}\rightarrow\textrm{lower}}\,(\epsilon_\textrm{upper} - \epsilon_\textrm{lower}) / {\cal D}_{i}\\
    p_\textrm{internal down}&= {\cal R}_{\textrm{i}\rightarrow\textrm{lower}}\,\epsilon_\textrm{lower}/{\cal D}_{i}\\,
    p_\textrm{internal up}&={\cal R}_{\textrm{i}\rightarrow\textrm{upper}}\,\epsilon_{i}/{\cal D}_{i}\\,

where :math:`i` is the current level, :math:`\epsilon` is the energy of the level, and :math:`{\cal R}` is the radiative
 rates.


We ignore the probability to emit a k-packet as TARDIS only works with photon packets.
Next we calculate the radiative
rates using equation 10 in :cite:`Lucy2003`.

.. math::
    {\cal R}_{\textrm{upper}\rightarrow\textrm{lower}} &=
    A_{\textrm{upper}\rightarrow\textrm{lower}}\beta_\textrm{Sobolev}n_\textrm{upper}\\
    {\cal R}_{\textrm{lower}\rightarrow\textrm{upper}} &=
    (B_{\textrm{lower}\rightarrow\textrm{upper}}n_\textrm{lower}-
    B_{\textrm{upper}\rightarrow\textrm{lower}}n_\textrm{upper})
    \beta_\textrm{Sobolev} J_{\nu}^{b}\\,

with :math:`\beta_\textrm{Sobolev} = \frac{1}{\tau_\textrm{Sobolev}}(1-e^{-\tau_\textrm{Sobolev}})` .

using the Einstein coefficients

.. math::
    A_{\textrm{upper}\rightarrow\textrm{lower}} &= \frac{8 \nu^2 \pi^2 e^2}{m_e c^3}~
        \frac{g_\textrm{lower}}{g_\textrm{upper}}~f_{\textrm{lower}\rightarrow\textrm{upper}}\\
    A_{\textrm{upper}\rightarrow\textrm{lower}} &= \underbrace{\frac{4 \pi^2 e^2}{m_e c}}_{C_\textrm{Einstein}}~ \frac{2\nu^2}{c^2}
            \frac{g_\textrm{lower}}{g_\textrm{upper}}~f_{\textrm{lower}\rightarrow\textrm{upper}}\\
    B_{\textrm{lower}\rightarrow\textrm{upper}} &= \frac{4\pi^2 e^2}{m_e h\nu c}\,f_{\textrm{lower}\rightarrow\textrm{upper}}\\
    B_{\textrm{lower}\rightarrow\textrm{upper}} &= \underbrace{\frac{4 \pi^2 e^2}{m_e c}}_{C_\textrm{Einstein}}\frac{1}{h\nu} f_{\textrm{lower}\rightarrow\textrm{upper}}\\

    B_{\textrm{upper}\rightarrow\textrm{lower}} &= \frac{4\pi^2 e^2}{m_e h\nu c}\,f_{\textrm{lower}\rightarrow\textrm{upper}}\\
    B_{\textrm{upper}\rightarrow\textrm{lower}} &= \underbrace{\frac{4 \pi^2 e^2}{m_e c}}_{C_\textrm{Einstein}}\frac{1}{h\nu}\frac{g_\textrm{lower}}{g_\textrm{upper}}f_{\textrm{lower}\rightarrow\textrm{upper}}\\

we get

.. math::
    {\cal R}_{\textrm{upper}\rightarrow\textrm{lower}} &=
        C_\textrm{Einstein} \frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{upper}}~f_{\textrm{lower}\rightarrow\textrm{upper}}
        \beta_\textrm{Sobolev}n_\textrm{upper}\\

    {\cal R}_{\textrm{lower}\rightarrow\textrm{upper}} &=
            C_\textrm{Einstein}\frac{1}{h\nu} f_{\textrm{lower}\rightarrow\textrm{upper}}
            (n_\textrm{lower}-\frac{g_\textrm{lower}}{g_\textrm{upper}}n_\textrm{upper})
                        \beta_\textrm{Sobolev} J_{\nu}^{b}\\

This results in the transition probabilities:

.. math::
    p_\textrm{emission down}&= C_\textrm{Einstein} \frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}n_\textrm{i}\,(\epsilon_\textrm{i} - \epsilon_\textrm{lower}) / {\cal D}_{i}\\
    p_\textrm{internal down}&= C_\textrm{Einstein} \frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}n_\textrm{i}\,\epsilon_\textrm{lower}/{\cal D}_{i}\\
    p_\textrm{internal up}&=C_\textrm{Einstein}\frac{1}{h\nu} f_{\textrm{i}\rightarrow\textrm{upper}}
                                        (n_\textrm{i}-\frac{g_\textrm{i}}{g_\textrm{upper}}n_\textrm{upper})
                                                    \beta_\textrm{Sobolev} J_{\nu}^{b}\,\epsilon_{i}/{\cal D}_{i}\\,

and as we will normalise the transition probabilities numerically later,  we can get rid of :math:`C_\textrm{Einstein}`,
 :math:`\frac{1}{{\cal D}_i}` and number densities.

.. math::
    p_\textrm{emission down}&= \frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}\,(\epsilon_\textrm{i} - \epsilon_\textrm{lower})\\
    p_\textrm{internal down}&=  \frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                       \beta_\textrm{Sobolev}\,\epsilon_\textrm{lower}\\
    p_\textrm{internal up}&=\frac{1}{h\nu} f_{\textrm{i}\rightarrow\textrm{upper}}
                                        \underbrace{(1-\frac{g_\textrm{i}}{g_\textrm{upper}}\frac{n_\textrm{upper}}{n_i})}
                                        _\textrm{stimulated emission}
                                                    \beta_\textrm{Sobolev} J_{\nu}^{b}\,\epsilon_{i}\\,




There are two parts for each of the probabilities, one that is pre-computed by `~tardisatomic` and is in the HDF5 File,
and one that is computed during the plasma calculations:

.. math::
        p_\textrm{emission down}&= \underbrace{\frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                           (\epsilon_\textrm{i} - \epsilon_\textrm{lower})}_\textrm{pre-computed}
                                           \,\beta_\textrm{Sobolev}\\
        p_\textrm{internal down} &= \underbrace{\frac{2\nu^2}{c^2} \frac{g_\textrm{lower}}{g_\textrm{i}}~f_{\textrm{lower}\rightarrow\textrm{i}}
                                           \epsilon_\textrm{lower}}_\textrm{pre-computed}\,\beta_\textrm{Sobolev}\\
        p_\textrm{internal up} &= \underbrace{\frac{1}{h\nu} f_{\textrm{i}\rightarrow\textrm{upper}}}_\textrm{pre-computed}
                                                        \beta_\textrm{Sobolev} J_{\nu}^{b}\,
                                                        (1-\frac{g_\textrm{i}}{g_\textrm{upper}}\frac{n_\textrm{upper}}{n_i})
                                                        \,\epsilon_{i}.
