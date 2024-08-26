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

.. |CITATION| replace:: kerzendorf_2024_13370472

.. |DOI_BADGE| image:: https://img.shields.io/badge/DOI-10.5281/zenodo.13370472-blue
                 :target: https://doi.org/10.5281/zenodo.13370472

.. code-block:: bibtex

    @software{kerzendorf_2024_13370472,
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
                      Arya, Atharva and
                      Smith, Isaac and
                      Cawley, Kevin and
                      Singhal, Jaladh and
                      Shields, Joshua and
                      O'Brien, Jack and
                      Barbosa, Talytha and
                      Sondhi, Dhruv and
                      Yu, Jenny and
                      Patel, Maryam and
                      Shields, Josh and
                      Varanasi, Kaushik and
                      Rathi, Shikha and
                      Chitchyan, Sona and
                      Gillanders, James and
                      Singh, Shreyas and
                      Savel, Arjun and
                      Gupta, Sumit and
                      Reinecke, Martin and
                      Holas, Alexander and
                      Eweis, Youssef and
                      Bylund, Tomas and
                      Black, William and
                      Bentil, Laud and
                      Kumar, Asish and
                      Eguren, Jordi and
                      Kumar, Ansh and
                      Bartnik, Matthew and
                      Alam, Arib and
                      Magee, Mark and
                      Dutta, Anirban and
                      Srivastava, Sarthak and
                      Varma Buddaraju, Rohith and
                      Visser, Erin and
                      Daksh, Ayushi and
                      Lu, Jing and
                      Livneh, Ran and
                      Kambham, Satwik and
                      Roldan, Israel and
                      Bhakar, Jayant and
                      Mishra, Sashank and
                      Rajagopalan, Srinath and
                      Reichenbach, John and
                      Jain, Rinkle and
                      Actions, GitHub and
                      Floers, Andreas and
                      Gupta, Harshul and
                      Chaumal, Aarya and
                      Brar, Antreev and
                      Singh, Sourav and
                      Kowalski, Nathan and
                      Patidar, Abhishek and
                      Matsumura, Yuki and
                      Selsing, Jonatan and
                      Sofiatti, Caroline and
                      Talegaonkar, Chinmay and
                      Kumar, Aman and
                      Sharma, Sampark and
                      Buchner, Johannes and
                      Yap, Kevin and
                      Martinez, Laureano and
                      Truong, Le and
                      Zingale, Michael and
                      Sandler, Morgan and
                      Zaheer, Musabbiha and
                      Sarafina, Nance and
                      Dasgupta, Debajyoti and
                      Patra, Nilesh and
                      Singh Rathore, Parikshit and
                      Patel, Pratik and
                      Volodin, Dmitry and
                      Venkat, Shashank and
                      Prasad, Shilpi and
                      Gupta, Suyash and
                      Lemoine, Thom and
                      Wahi, Ujjwal and
                      Aggarwal, Yash and
                      Chen, Nutan and
                      Kolliboyina, Chaitanya and
                      PATIDAR, ABHISHEK and
                      Nayak U, Ashwin and
                      Kumar, Atul and
                      Kharkar, Atharwa},
      title        = {tardis-sn/tardis: TARDIS v2024.08.25},
      month        = aug,
      year         = 2024,
      publisher    = {Zenodo},
      version      = {release-2024.08.25},
      doi          = {10.5281/zenodo.13370472},
      url          = {https://doi.org/10.5281/zenodo.13370472}
    }

