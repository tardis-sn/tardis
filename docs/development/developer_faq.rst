*************
Developer FAQ
*************

Constants in TARDIS are all taken from Astropy. The module tardis.constants import all constants currently
from astropy.constants.astropy13constants.

Class design and inheritance:

 * If only constructor changed -> use classmethod
 * if overriding other methods -> subclass

We use black to check PEP8 compliances. 
