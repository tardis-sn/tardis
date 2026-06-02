***********************
Bibliography and Glossary
***********************


Bibliography
==========

.. bibliography:: ../tardis.bib


Glossary
========

.. glossary:: 
        
    Chianti
        CHIANTI consists of a critically evaluated set of up-to-date atomic data, together with user-friendly programs written in Interactive Data Language (IDL) and Python to calculate the spectra from astrophysical plasmas.
    
    Estimator
        A statistical tool or algorithm used in Monte Carlo simulations to calculate physical quantities. In the case of TARDIS this is usually the estimation of radiation field properties based on the behavior and distribution of photon packets.
    
    Grey Opacity
        A simplified opacity model where the absorption coefficient is assumed to be independent of photon wavelength or frequency, making calculations more tractable.
    
    HDF (Hierarchical Data Format)
        A file format designed for storing and organizing large amounts of scientific data in a hierarchical structure, commonly used with the .h5 extension.
    
    Kurucz
        A comprehensive database of atomic line data compiled by Robert Kurucz, widely used in stellar and supernova atmosphere modeling for opacity calculations.
    
    Meta-Stable
        Metastability is the condition of a system where the system has stability, but is not as stable as in the system's state of least energy. 
        In atomic physics, this is term is usually used to describe excitation states that have long spontaneous emission timescales, which corresponds to low oscillator strengths of transitions away from the metastable state. 
    
    Monte Carlo
        A computational method that uses random sampling to solve complex physical problems, particularly useful for simulating particle transport and radiative transfer.
    
    Nebular
        Relating to the nebular phase of supernovae, typically occurring weeks to months after explosion when the ejecta becomes optically thin and emission lines dominate the spectrum. 
        In TARDIS, this may also refer to adjustments to LTE assumptions to handle lower density environments, which are often applicable to supernova ejecta. 
    
    Opacity
        A measure of how opaque or transparent a material is to electromagnetic radiation, quantifying the probability of photon absorption or scattering per unit path length.
    
    Packets
        In Monte Carlo radiative transfer simulations, discrete bundles of energy or photons that are tracked as they propagate through the computational domain and interact with matter.
    
    Seed
        A numerical value used to initialize random number generators, ensuring reproducibility of Monte Carlo simulations when the same seed is used.
    
    Synapps
        SYNAPPS is an open-source spectrum fitter embedding a highly parameterized synthetic SN spectrum calculation within a parallel asynchronous optimizer, created to systematically interpret large sets of SN spectroscopy data.
    
    TOML
        TOML (Tom's Obvious, Minimal Language) is a minimal configuration file format that is designed to be easy to read due to obvious semantics. It is designed to map unambiguously to a hash table and to be easy to parse into data structures in a wide variety of languages. TOML files have the ending ".toml".
    
    YAML
        YAML (YAML Ain't Markup Language) is a human friendly data serialization standard for all programming languages. It is commonly used for configuration files and in applications where data is being stored or transmitted. YAML files have the ending ".yml".