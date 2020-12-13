******
TARDIS
******

.. image:: https://img.shields.io/badge/read-documentation-blue
  :target: https://tardis-sn.github.io/tardis

.. image:: https://dev.azure.com/tardis-sn/TARDIS/_apis/build/status/tardis-sn.tardis?branchName=master
  :target: https://dev.azure.com/tardis-sn/TARDIS/_build/latest?definitionId=1&branchName=master

.. image:: https://codecov.io/gh/tardis-sn/tardis/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/tardis-sn/tardis

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3893940.svg
   :target: https://doi.org/10.5281/zenodo.3893940

.. image:: https://badges.gitter.im/Join%20Chat.svg
  :target: https://gitter.im/tardis-sn/tardis

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black

.. image:: https://img.shields.io/badge/Donate-to%20TARDIS-brightgreen.svg
    :target: https://numfocus.salsalabs.org/donate-to-tardis/index.html
    
.. image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
    :target: http://numfocus.org
    :alt: Powered by NumFOCUS

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org

TARDIS is a tool that creates synthetic observations (spectra) for exploding
stars (supernovae).

******************************
Credits & Publication Policies
******************************

We provide TARDIS as a free, open-source tool. If you are using it, please
adhere to a few policies and acknowledge the TARDIS Team.

Publication Policies
====================

If you use this code for any publications or presentations please acknowledge
it.  Please cite `Kerzendorf & Sim 2014
<http://adsabs.harvard.edu/abs/2014MNRAS.440..387K>`_  in the text

Please add this paragraph to the Acknowledgement:

.. code-block:: none

    This research made use of \textsc{Tardis}, a community-developed software
    package for spectral synthesis in supernovae
    \citep{2014MNRAS.440..387K, kerzendorf_wolfgang_2019_2590539}.
    The development of \textsc{Tardis} received support from the
    Google Summer of Code initiative
    and from ESA's Summer of Code in Space program. \textsc{Tardis} makes
    extensive use of Astropy and PyNE.



If you use any of the full relativity treatments or use TARDIS for
modelling Type II supernovae you also add this citation to acknowledgement
`Spectral modeling of type II supernovae. I. Dilution factors
<https://ui.adsabs.harvard.edu/abs/2019A%26A...621A..29V>`_:

.. code-block:: none

    \citep{2019A&A...621A..29V}

The following bibtex entries are needed for the references.

.. code-block:: none

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

    @ARTICLE{2019A&A...621A..29V,
           author = {{Vogl}, C. and {Sim}, S.~A. and {Noebauer}, U.~M. and
             {Kerzendorf}, W.~E. and {Hillebrandt}, W.},
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


    @software{kerzendorf_wolfgang_2020_3893940,
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
                      Barbosa, Talytha and
                      Patel, Maryam and
                      Varanasi, Kaushik and
                      Eweis, Youssef and
                      Reinecke, Martin and
                      Bylund, Tomas and
                      Bentil, Laud and
                      Eguren, Jordi and
                      Livneh, Ran and
                      Singhal, Jaladh and
                      O'Brien, Jack and
                      Rajagopalan, Srinath and
                      Jain, Rinkle and
                      Reichenbach, John and
                      Mishra, Sashank and
                      Singh, Sourav and
                      Sofiatti, Caroline and
                      Selsing, Jonatan and
                      Kowalski, Nathan and
                      Savel, Arjun and
                      Talegaonkar, Chinmay and
                      Patel, Pratik and
                      Patra, Nilesh and
                      Nayak, Ashwin and
                      Kumar, Atul and
                      Sarafina, Nance and
                      Gillanders, James and
                      Sharma, Sampark and
                      Wahi, Ujjwal and
                      Dasgupta, Debajyoti and
                      Magee, Mark and
                      Yap, Kevin and
                      Gupta, Suyash},
      title        = {tardis-sn/tardis: TARDIS v3.0.dev3459},
      month        = jun,
      year         = 2020,
      publisher    = {Zenodo},
      version      = {v3.0.dev3459},
      doi          = {10.5281/zenodo.3893940},
      url          = {https://doi.org/10.5281/zenodo.3893940}
    }
