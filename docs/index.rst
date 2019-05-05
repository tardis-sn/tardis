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
of supernovae. It is developed and maintained by a multi-disciplinary team
including software engineers, computer scientists, statisticians,
and astrophysicists.

If you use this code for any publications or presentations please follow our
citation guidelines in :ref:`tardiscredits`

User modifications and additions that lead to publications need to be handed
back to the community by incorporating them into TARDIS.
Please contact the TARDIS team via the `github page
<https://github.com/tardis-sn/tardis>`_ if you have questions or need
assistance.

============
Using Tardis
============

.. toctree::
    :maxdepth: 2

    installation
    CHANGELOG.md
    quickstart
    running/index
    examples/index
    scripts/index
    credits

======================
Looking under the hood
======================

.. toctree::
    :maxdepth: 2

    atomic/atomic_data
    physics/index
    montecarlo/index


=================
Developing Tardis
=================

.. toctree::
    :maxdepth: 2

    issues
    workflow/development_workflow
    runnints_tests

==========
References
==========

.. toctree::
    zreferences

====
News
====

.. include:: news.rst

