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
    sponsored project of NumFOCUS. \\textsc{tardis} makes extensive use of Astropy.

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

.. |CITATION| replace:: kerzendorf_2026_20396341

.. |DOI_BADGE| image:: https://img.shields.io/badge/DOI-10.5281/zenodo.396341-blue
                 :target: https://doi.org/10.5281/zenodo.396341

.. code-block:: bibtex

    @software{kerzendorf_2026_20396341,
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
                      Sondhi, Dhruv and
                      Yu, Jenny and
                      O'Brien, Jack and
                      Shields, Josh and
                      Patel, Maryam and
                      Rathi, Shikha and
                      Varanasi, Kaushik and
                      Gillanders, James and
                      Chitchyan, Sona and
                      Gupta, Sumit and
                      Singh, Shreyas and
                      Marie Lynn, Haille and
                      Savel, Arjun and
                      Reinecke, Martin and
                      Eweis, Youssef and
                      Shah, Swayam and
                      Holas, Alexander and
                      Bylund, Tomas and
                      Visser, Erin and
                      Bentil, Laud and
                      Black, William and
                      Groneck, Ryan and
                      Lu, Jing and
                      Kumar, Asish and
                      Dutta, Anirban and
                      Eguren, Jordi and
                      Kumar, Ansh and
                      Srivastava, Sarthak and
                      Bartnik, Matthew and
                      Magee, Mark and
                      Alam, Arib and
                      Varma Buddaraju, Rohith and
                      Livneh, Ran and
                      Kambham, Satwik and
                      Daksh, Ayushi and
                      Mishra, Sashank and
                      Bhakar, Jayant and
                      Powers, Cecelia and
                      Roldan, Israel and
                      Rajagopalan, Srinath and
                      McClellan, Connor and
                      Reichenbach, John and
                      Nitish, P and
                      Actions, GitHub and
                      Jain, Rinkle and
                      Dadu, Aaryan and
                      Brar, Antreev and
                      Chaumal, Aarya and
                      Singh, Sourav and
                      Gupta, Harshul and
                      Kowalski, Nathan and
                      Gangbhoj, Riddhi and
                      Sofiatti, Caroline and
                      Talegaonkar, Chinmay and
                      Perkins, Haille and
                      Selsing, Jonatan and
                      Matsumura, Yuki and
                      Patidar, Abhishek and
                      Wahi, Ujjwal and
                      Aggarwal, Yash and
                      Patel, Pratik and
                      Singh Rathore, Parikshit and
                      L. Lim, P. and
                      Nagadevi, Kona and
                      Buchner, Johannes and
                      Bhandari, Jhalak and
                      Patra, Nilesh and
                      Yap, Kevin and
                      Truong, Le and
                      Chen, Nutan and
                      Zingale, Michael and
                      Sandler, Morgan and
                      Zaheer, Musabbiha and
                      Sarafina, Nance and
                      Vieira, Nicholas and
                      Gupta, Suyash and
                      Lemoine, Thom and
                      Kumar, Atul and
                      Saraf, Shreyans and
                      Nayak U, Ashwin and
                      Dasgupta, Debajyoti and
                      Jaiswal, Abhayraj and
                      Watson, Clyde and
                      Kumar, Aman and
                      Volodin, Dmitry and
                      Martinez, Laureano and
                      PATIDAR, ABHISHEK and
                      Diddige, Harshitha and
                      Rao, Rishmita and
                      Prasad, Rohit and
                      Gajanan Nalbalwar, Rudraksh and
                      Sharma, Sampark and
                      Venkat, Shashank and
                      Prasad, Shilpi},
      title        = {tardis-sn/tardis: TARDIS v2026.05.26},
      month        = may,
      year         = 2026,
      publisher    = {Zenodo},
      version      = {release-2026.05.26},
      doi          = {10.5281/zenodo.20396341},
      url          = {https://doi.org/10.5281/zenodo.20396341},
    }
