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

.. |CITATION| replace:: kerzendorf_wolfgang_2022_7331855

.. |DOI_BADGE| image:: https://img.shields.io/badge/DOI-10.5281/zenodo.7331855-blue
                 :target: https://doi.org/10.5281/zenodo.7331855

.. code-block:: bibtex

    @software{kerzendorf_wolfgang_2022_7331855,
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
                      Barbosa, Talytha and
                      Sondhi, Dhruv and
                      Yu, Jenny and
                      O'Brien, Jack and
                      Patel, Maryam and
                      Varanasi, Kaushik and
                      Gillanders, James and
                      Savel, Arjun and
                      Reinecke, Martin and
                      Eweis, Youssef and
                      Bylund, Tomas and
                      Bentil, Laud and
                      Chitchyan, Sona and
                      Eguren, Jordi and
                      Alam, Arib and
                      Bartnik, Matthew and
                      Varma Buddaraju, Rohith and
                      Magee, Mark and
                      Shields, Joshua and
                      Livneh, Ran and
                      Kambham, Satwik and
                      Mishra, Sashank and
                      Rajagopalan, Srinath and
                      Jain, Rinkle and
                      Reichenbach, John and
                      Floers, Andreas and
                      Holas, Alexander and
                      Bhakar, Jayant and
                      Brar, Antreev and
                      Singh, Sourav and
                      Kowalski, Nathan and
                      Kumar, Aman and
                      Talegaonkar, Chinmay and
                      Sofiatti, Caroline and
                      Selsing, Jonatan and
                      Venkat, Shashank and
                      Patra, Nilesh and
                      Patel, Pratik and
                      Volodin, Dmitry and
                      Singh Rathore, Parikshit and
                      Yap, Kevin and
                      Sarafina, Nance and
                      Zaheer, Musabbiha and
                      Sandler, Morgan and
                      Gupta, Suyash and
                      Prasad, Shilpi and
                      Kumar, Atul and
                      Lemoine, Thom and
                      Wahi, Ujjwal and
                      Aggarwal, Yash and
                      Sharma, Sampark and
                      Dasgupta, Debajyoti and
                      Martinez, Laureano and
                      Kolliboyina, Chaitanya and
                      Nayak U, Ashwin and
                      Kharkar, Atharwa},
      title        = {tardis-sn/tardis: TARDIS v2022.11.17},
      month        = nov,
      year         = 2022,
      publisher    = {Zenodo},
      version      = {release-2022.11.17},
      doi          = {10.5281/zenodo.7331855},
      url          = {https://doi.org/10.5281/zenodo.7331855}
    }

