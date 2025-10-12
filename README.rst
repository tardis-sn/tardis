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

.. image:: https://github.com/tardis-sn/tardis/actions/workflows/docstr-cov.yml/badge.svg
    :target: https://github.com/tardis-sn/tardis/actions/workflows/docstr-cov.yml

.. image:: https://github.com/tardis-sn/tardis/actions/workflows/tests.yml/badge.svg
    :target: https://github.com/tardis-sn/tardis/actions/workflows/tests.yml

.. image:: https://github.com/tardis-sn/tardis/actions/workflows/build-docs.yml/badge.svg
    :target: https://tardis-sn.github.io/tardis/index.html

.. image:: https://img.shields.io/badge/benchmarks-passing-green
    :target: https://github.com/tardis-sn/tardis/actions/workflows/benchmarks.yml

.. image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/tardis-sn/tardis/master/docs/_static/ruff_badge.json
    :target: https://github.com/tardis-sn/tardis/actions/workflows/codestyle.yml
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

.. |CITATION| replace:: kerzendorf_2025_17330175

.. |DOI_BADGE| image:: https://img.shields.io/badge/DOI-10.5281/zenodo.17330175-blue
                 :target: https://doi.org/10.5281/zenodo.17330175

.. code-block:: bibtex

    @software{kerzendorf_2025_17330175,
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
                      Arya, Atharva and
                      Fullard, Andrew and
                      Smith, Isaac and
                      Shields, Joshua and
                      Cawley, Kevin and
                      Singhal, Jaladh and
                      Barbosa, Talytha and
                      O'Brien, Jack and
                      Sondhi, Dhruv and
                      Yu, Jenny and
                      Shields, Josh and
                      Patel, Maryam and
                      Rathi, Shikha and
                      Varanasi, Kaushik and
                      Chitchyan, Sona and
                      Gillanders, James and
                      Marie Lynn, Haille and
                      Gupta, Sumit and
                      Singh, Shreyas and
                      Savel, Arjun and
                      Reinecke, Martin and
                      Eweis, Youssef and
                      Shah, Swayam and
                      Holas, Alexander and
                      Visser, Erin and
                      Bylund, Tomas and
                      Bentil, Laud and
                      Black, William and
                      Lu, Jing and
                      Kumar, Asish and
                      Groneck, Ryan and
                      Dutta, Anirban and
                      Eguren, Jordi and
                      Bartnik, Matthew and
                      Srivastava, Sarthak and
                      Kumar, Ansh and
                      Alam, Arib and
                      Magee, Mark and
                      Varma Buddaraju, Rohith and
                      Daksh, Ayushi and
                      Kambham, Satwik and
                      Livneh, Ran and
                      Mishra, Sashank and
                      Rajagopalan, Srinath and
                      Bhakar, Jayant and
                      Powers, Cecelia and
                      Roldan, Israel and
                      Nitish, P and
                      Reichenbach, John and
                      Jain, Rinkle and
                      Actions, GitHub and
                      Gupta, Harshul and
                      Brar, Antreev and
                      Singh, Sourav and
                      Chaumal, Aarya and
                      Selsing, Jonatan and
                      Sofiatti, Caroline and
                      Talegaonkar, Chinmay and
                      Kowalski, Nathan and
                      Matsumura, Yuki and
                      Gangbhoj, Riddhi and
                      Patidar, Abhishek and
                      Perkins, Haille and
                      Aggarwal, Yash and
                      Patel, Pratik and
                      Singh Rathore, Parikshit and
                      L. Lim, P. and
                      Yap, Kevin and
                      Bhandari, Jhalak and
                      Buchner, Johannes and
                      Patra, Nilesh and
                      Martinez, Laureano and
                      Chen, Nutan and
                      Zingale, Michael and
                      Sandler, Morgan and
                      Zaheer, Musabbiha and
                      Sarafina, Nance and
                      Vieira, Nicholas and
                      Lemoine, Thom and
                      Wahi, Ujjwal and
                      Kumar, Aman and
                      Gupta, Suyash and
                      Truong, Le and
                      Nayak U, Ashwin and
                      Kumar, Atul and
                      Volodin, Dmitry and
                      Dasgupta, Debajyoti and
                      McClellan, Connor and
                      PATIDAR, ABHISHEK and
                      Diddige, Harshitha and
                      Prasad, Rohit and
                      Gajanan Nalbalwar, Rudraksh and
                      Sharma, Sampark and
                      Venkat, Shashank and
                      Prasad, Shilpi},
      title        = {tardis-sn/tardis: TARDIS v2025.10.12},
      month        = oct,
      year         = 2025,
      publisher    = {Zenodo},
      version      = {release-2025.10.12},
      doi          = {10.5281/zenodo.17330175},
      url          = {https://doi.org/10.5281/zenodo.17330175},
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

