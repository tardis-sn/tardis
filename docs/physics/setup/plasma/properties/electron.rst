*******************
Electron Properties
*******************


Beta Radiation
==============
Calculated From: Radiative Temperature
Equation: :math:`\beta_{rad} = \frac{1}{k_b T_{rad}}`

Within the context of Statistical Mechanics, it is often best to derive a scale factor for certain properties. One of the more well-known of these scale factors is that between energy and temperature

.. math::
    \tau = k_b T
    
where :math:`T` is temperature in Kelvin and :math:`k_b` is the Boltzmann Constant, which is aproximately :math:`1.381 \times 10^{-23}` Joules per Kelvin or :math:`8.167 \times 10^{-5}` electronvolts per Kelvin. This scale factor can help define many different quantities, such as the partition function of a system, the total energy of a system, pressure, and many more. Sometimes, it is best to define a similar quantity

.. math::
    \beta = \dfrac{1}{\tau} = \dfrac{1}{k_b T}
    
which has units :math:`J^{-1}`. TARDIS will calculate this value based on the Radiative Temperature given and implement it whenever necessary. 

G Electron
==========
Calculated From: Beta Radiation
Equation: :math:`g_{electron} = (\frac{2\pi m_e / \beta_{rad}}{h^2})^{3/2}`

Electron Temperature
====================
Calculated From: LinkTRadTElectron, Radiative Temperature
Equation: :math:`T_{electron} = \text{LinkTRadTElectron} \times T_{rad}`

The electron temperature can be defined as the average temperature of a distribution of the velocities of a group of electrons if the velocities follow a Maxwell-Boltzmann distribution or two-thirds of the average energy if the system is not in equilibrium. As described in :cite:`MazzaliLucy93`, in assuming isothermal flow of stellar winds, it is sufficient to simply let 

.. math::
    T_{electron} = 0.9 \times T_{rad}

(see LinkTRadTElectron for more information).
TARDIS will calculate this value and use it if necessary since this temperature can be several orders of magnitude higher than the temperatures of other ions or even neutral atoms as a result of electrons being more susceptible to changes in heat than in ions and energy transfers between electrons being more efficient than atoms or ions due to the mass difference between eletrons and atoms or ions. 
