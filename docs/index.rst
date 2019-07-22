..
  .. image:: graphics/tardis_logo.jpg

.. the "raw" directive below is used to hide the title in favor of just the logo being visible
.. .. raw:: html

..    <style media="screen" type="text/css">
..      h1 { display:none; }
..    </style>

##################################
Tardis Core Package Documentation
##################################

.. image:: graphics/tardis_banner.svg


TARDIS is an open-source Monte Carlo radiative-transfer spectral synthesis code
for 1D models of supernova ejecta. It is designed for rapid spectral modelling
of supernovae. It is developed and maintained by a :ref:`multi-disciplinary team <team>`
including software engineers, computer scientists, statisticians,
and astrophysicists.

If you use this code for any publications or presentations please follow our
citation guidelines in :ref:`tardiscredits`

User modifications and additions that lead to publications **need to be handed
back to the community** by incorporating them into TARDIS.
Please contact the TARDIS team via the `github page
<https://github.com/tardis-sn/tardis>`_ if you have questions or need
assistance.



.. toctree::
    :maxdepth: 2
    :hidden:

    installation
    quickstart/quickstart



.. toctree::
    :maxdepth: 3
    :caption: Using TARDIS
    :hidden:


    running/index


.. toctree::
    :maxdepth: 2
    :caption: The Physics of TARDIS
    :hidden:

    Physics overview <physics/index>
    physics/plasma/index
    physics/montecarlo/index


.. toctree::
    :maxdepth: 2
    :caption: Research with TARDIS
    :hidden:

    Overview <research/index>


.. toctree::
    :maxdepth: 2
    :caption: Team & Credits
    :hidden:

    team
    credits


.. toctree::
    :maxdepth: 2
    :caption: Developers Guide
    :hidden:

    development/index
    CHANGELOG.md

.. toctree::
    :maxdepth: 2
    :caption: API
    :hidden:

    api/modules



.. toctree::
    :caption: References
    :hidden:

    zreferences

====
News
====

.. include:: news.rst

