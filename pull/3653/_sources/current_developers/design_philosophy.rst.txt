TARDIS Design Philosophy
=======================

Exception handling
^^^^^^^^^^^^^^^^^^

As research software, we expect the user to correctly input parameters and respond 
to failures. Thus, we do not want to raise exceptions for every possible edge case, 
as this would make the code unnecessarily verbose and difficult to read. 

Input validation should raise exceptions only when the input is clearly invalid, 
such as a negative number for a parameter that must be positive. Otherwise, we 
should raise warnings for inputs that may invalidate physical assumptions.

Only raise exceptions in general when the code is unable to continue with
computations, not when the result may just be unphysical.
For example, if tau Sobolevs are nan, inf or -inf, we should raise an exception.

General notes
^^^^^^^^^^^^^

Constants in TARDIS are taken from Astropy. The ``tardis.constants`` module
imports all constants from ``astropy.constants.astropy13constants``.

Class design and inheritance:

- If only the constructor changed, use a classmethod.
- If overriding other methods, use a subclass.
- Prefer composition over inheritance.

Use a classmethod when construction changes but behavior stays the same, such as
adding ``from_config(...)`` or ``from_csvy(...)`` constructors around an existing
class. Use a subclass when method behavior changes after construction, such as a
specialized solver or reader that overrides parsing or calculation methods.

Dos and Don'ts
^^^^^^^^^^^^^^

Don't write functions that are used once. It is okay to have a longer function to make
flow easy to follow.
Don't write functions that perform extremely simple tasks, especially if they are
for programming purposes rather than a specific calculation.
Don't write functions that just call one other function. This makes flow difficult
to follow.
Clean up after yourself- delete functions that are no longer used, and remove code that is commented out.
Avoid long if-else chains. Prefer dictionaries or lists.
Use dependencies that we already have wherever possible. For example, use ``numpy`` or ``pandas`` instead of writing your own array manipulation functions.
Don't write test helper functions. Use pytest fixtures to access data or states.
TARDIS does not have a public API. All functions and classes are internal. Methods do not need to be marked private with a leading underscore.
