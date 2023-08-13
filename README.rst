======
TARDIS
======

.. image:: https://img.shields.io/badge/Donate-to%20TARDIS-brightgreen.svg
    :target: https://numfocus.salsalabs.org/donate-to-tardis/index.html

.. image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
    :target: http://numfocus.org

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/tardis-sn/tardis

.. image:: https://img.shields.io/static/v1?logo=visualstudiocode&label=&message=Open%20in%20Visual%20Studio%20Code&labelColor=2c2c32&color=007acc&logoColor=007acc
    :target: https://open.vscode.dev/tardis-sn/tardis
|

TARDIS is a tool that creates synthetic observations (*spectra*) for exploding
stars (*supernovae*).

.. image:: https://codecov.io/gh/tardis-sn/tardis/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/tardis-sn/tardis

.. image:: https://img.shields.io/endpoint?url=https://jsonbin.org/tardis-bot/tardis/badges/docstr-cov
    :target: https://github.com/tardis-sn/tardis/actions/workflows/docstr-cov.yml?query=branch%3Amaster

.. image:: https://github.com/tardis-sn/tardis/actions/workflows/tests.yml/badge.svg
    :target: https://github.com/tardis-sn/tardis/actions/workflows/tests.yml

.. image:: https://github.com/tardis-sn/tardis/actions/workflows/build-docs.yml/badge.svg
    :target: https://tardis-sn.github.io/tardis/index.html

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
|

.. _tardiscredits:

******************************
Credits & Publication Policies
******************************

|DOI_BADGE|

We provide TARDIS as a free, open-source tool. If you are using it, please
adhere to a few policies and acknowledge the TARDIS Team.

Publication Policies
====================

If you use this code for any publications or presentations please acknowledge
it.  Please cite `Kerzendorf & Sim 2014
<http://adsabs.harvard.edu/abs/2014MNRAS.440..387K>`_  in the text and add the
following paragraph to the Acknowledgement section:

.. parsed-literal::

    This research made use of \\textsc{tardis}, a community-developed software package for spectral
    synthesis in supernovae \\citep{2014MNRAS.440..387K, |CITATION|}. The
    development of \\textsc{tardis} received support from GitHub, the Google Summer of Code
    initiative, and from ESA's Summer of Code in Space program. \\textsc{tardis} is a fiscally
    sponsored project of NumFOCUS. \\textsc{tardis} makes extensive use of Astropy and Pyne.

If you use any of the full relativity treatments or use TARDIS for modelling
Type II supernovae, also add `Spectral modeling of type II supernovae. I. Dilution factors <https://ui.adsabs.harvard.edu/abs/2019A%26A...621A..29V>`_
to the Acknowledgement.

.. parsed-literal::

    \citep{2019A&A...621A..29V}

The following BibTeX entries are needed for the references:

.. code-block:: bibtex

    @ARTICLE{2014MNRAS.440..387K,
           author = {{Kerzendorf}, W.~E. and {Sim}, S.~A.},
            title = "{A spectral synthesis code for rapid modelling of supernovae}",
          journal = {\mnras},
    archivePrefix = "arXiv",
           eprint = {1401.5469},
     primaryClass = "astro-ph.SR",
         keywords = {radiative transfer, methods: numerical, supernovae: general},
             year = 2014,
            month = may,
           volume = 440,
            pages = {387-404},
              doi = {10.1093/mnras/stu055},
           adsurl = {http://adsabs.harvard.edu/abs/2014MNRAS.440..387K},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

.. code-block:: bibtex

    @ARTICLE{2019A&A...621A..29V,
           author = {{Vogl}, C. and {Sim}, S.~A. and {Noebauer}, U.~M. and {Kerzendorf}, W.~E. and {Hillebrandt}, W.},
            title = "{Spectral modeling of type II supernovae. I. Dilution factors}",
          journal = {\aap},
         keywords = {radiative transfer, methods: numerical, stars: distances, supernovae: general, supernovae: individual: SN1999em, Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Solar and Stellar Astrophysics},
             year = "2019",
            month = "Jan",
           volume = {621},
              eid = {A29},
            pages = {A29},
              doi = {10.1051/0004-6361/201833701},
    archivePrefix = {arXiv},
           eprint = {1811.02543},
     primaryClass = {astro-ph.HE},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2019A&A...621A..29V},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

.. |CITATION| replace:: kerzendorf_wolfgang_2023_8244935

.. |DOI_BADGE| image:: https://img.shields.io/badge/DOI-10.5281/zenodo.8244935-blue
                 :target: https://doi.org/10.5281/zenodo.8244935

.. code-block:: bibtex

    @software{kerzendorf_wolfgang_2023_8244935,
      author       = {Kerzendorf, Wolfgang and
                      Sim, Stuart and
                      Vogl, Christian and
                      Williamson, Marc and
                      Pássaro, Ezequiel and
                      Flörs, Andreas and
                      Camacho, Yssa and
                      Jančauskas, Vytautas and
                      Harpole, Alice and
                      Nöbauer, Ulrich and
                      Lietzau, Stefan and
                      Mishin, Mikhail and
                      Tsamis, Fotis and
                      Boyle, Aoife and
                      Shingles, Luke and
                      Gupta, Vaibhav and
                      Desai, Karan and
                      Klauser, Michael and
                      Beaujean, Frederik and
                      Suban-Loewen, Adam and
                      Heringer, Epson and
                      Barna, Barnabás and
                      Gautam, Gaurav and
                      Fullard, Andrew and
                      Smith, Isaac and
                      Cawley, Kevin and
                      Singhal, Jaladh and
                      Arya, Atharva and
                      O'Brien, Jack and
                      Barbosa, Talytha and
                      Sondhi, Dhruv and
                      Yu, Jenny and
                      Patel, Maryam and
                      Varanasi, Kaushik and
                      Gillanders, James and
                      Rathi, Shikha and
                      Chitchyan, Sona and
                      Savel, Arjun and
                      Reinecke, Martin and
                      Eweis, Youssef and
                      Bylund, Tomas and
                      Bentil, Laud and
                      Black, William and
                      Shields, Joshua and
                      Eguren, Jordi and
                      Alam, Arib and
                      Kumar, Ansh and
                      Bartnik, Matthew and
                      Magee, Mark and
                      Singh, Shreyas and
                      Varma Buddaraju, Rohith and
                      Livneh, Ran and
                      Kambham, Satwik and
                      Rajagopalan, Srinath and
                      Daksh, Ayushi and
                      Mishra, Sashank and
                      Jain, Rinkle and
                      Reichenbach, John and
                      Floers, Andreas and
                      Actions, GitHub and
                      Holas, Alexander and
                      Singh, Sourav and
                      Brar, Antreev and
                      Chaumal, Aarya and
                      Bhakar, Jayant and
                      Selsing, Jonatan and
                      Kowalski, Nathan and
                      Kumar, Aman and
                      Patidar, Abhishek and
                      Talegaonkar, Chinmay and
                      Sofiatti, Caroline and
                      Venkat, Shashank and
                      Sharma, Sampark and
                      Sarafina, Nance and
                      Patel, Pratik and
                      Singh Rathore, Parikshit and
                      Patra, Nilesh and
                      Lu, Jing and
                      Zaheer, Musabbiha and
                      Sandler, Morgan and
                      Truong, Le and
                      Yap, Kevin and
                      Buchner, Johannes and
                      Gupta, Suyash and
                      Prasad, Shilpi and
                      Kolliboyina, Chaitanya and
                      Lemoine, Thom and
                      Wahi, Ujjwal and
                      Aggarwal, Yash and
                      Matsumura, Yuki and
                      Gupta, Harshul and
                      Volodin, Dmitry and
                      PATIDAR, ABHISHEK and
                      Martinez, Laureano and
                      Kharkar, Atharwa and
                      Nayak U, Ashwin and
                      Dasgupta, Debajyoti and
                      Kumar, Atul},
      title        = {tardis-sn/tardis: TARDIS v2023.08.13},
      month        = aug,
      year         = 2023,
      publisher    = {Zenodo},
      version      = {release-2023.08.13},
      doi          = {10.5281/zenodo.8244935},
      url          = {https://doi.org/10.5281/zenodo.8244935}
    }

*******
License
*******

.. image:: https://img.shields.io/conda/l/conda-forge/tardis-sn
    :target: https://github.com/tardis-sn/tardis/blob/master/licenses/LICENSE.rst

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
|

This project is Copyright (c) TARDIS Collaboration and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause license. See the licenses folder for
more information.

