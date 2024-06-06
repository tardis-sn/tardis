..
  .. image:: graphics/tardis_logo.jpg

.. the "raw" directive below is used to hide the title in favor of just the logo being visible
.. .. raw:: html

..    <style media="screen" type="text/css">
..      h1 { display:none; }
..    </style>

##################################
TARDIS Core Package Documentation
##################################

.. image:: graphics/tardis_banner.svg

|

TARDIS is an open-source Monte Carlo radiative-transfer spectral synthesis code
for 1D models of supernova ejecta. It is designed for rapid spectral modelling
of supernovae. It is developed and maintained by a 
`multi-disciplinary team <https://tardis-sn.github.io/team/community_roles/>`_
including software engineers, computer scientists, statisticians,
and astrophysicists.

The tardis package is the core package provided by the TARDIS Project. To learn
more about the complete project, other packages it provides, and the community
behind it, please check the `project website <https://tardis-sn.github.io>`_.

If you use this code for any publications or presentations please follow our
citation guidelines in :ref:`tardiscredits`.

User modifications and additions that lead to publications **need to be handed
back to the community** by incorporating them into TARDIS.
Please contact the TARDIS team via the `GitHub page
<https://github.com/tardis-sn/tardis>`_ if you have questions or need
assistance.

-----------------
Mission Statement
-----------------

    *The mission of the TARDIS community is to develop open-source software 
    instruments to analyze and simulate astronomical transients 
    (supernovae, kilonovae, etc.) for research and education purposes. 
    We aim to build up a diverse group of researchers and developers 
    using an open-community model that emphasizes interdisciplinary 
    research and science reproducibility.*

.. toctree::
    :maxdepth: 2
    :hidden:

    installation
    quickstart
    tutorials
    how_to_guides
    faq
    API <api/modules>


.. toctree::
    :maxdepth: 3
    :caption: Input/Output
    :hidden:
    
    io/hdf/index
    io/configuration/index
    io/model/index
    io/optional/index
    io/visualization/index
    io/output/index
    io/grid/TardisGridTutorial

.. toctree::
    :maxdepth: 2
    :caption: Physics Walkthrough
    :hidden:
    
    physics/intro/index
    physics/setup/index
    physics/montecarlo/index
    physics/update_and_conv/update_and_conv
    physics/spectrum/index


.. toctree::
    :maxdepth: 2
    :caption: Contributing to TARDIS
    :hidden:

    contributing/CONTRIBUTING.md
    contributing/development/index
    contributing/tools/index
    contributing/CHANGELOG.md
    contributing/in_progress/index


.. toctree::
    :caption: Other Resources
    :hidden:
    
    resources/credits
    resources/research_done_using_TARDIS/research_papers
    resources/code_comparison/index
    resources/zreferences
