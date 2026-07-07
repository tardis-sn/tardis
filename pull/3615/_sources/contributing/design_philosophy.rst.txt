An important caveat to the above is that as research software, we expect the user 
to correctly input parameters and respond to failures. Thus, we do not want to 
raise exceptions for every possible edge case, as this would make the code unnecessarily 
verbose and difficult to read. Instead, we should only raise exceptions for edge 
cases that are likely to occur and can be easily handled by the user. For example, in the above code snippet, it is likely that a user might accidentally pass an invalid value to packets_mode, and it is easy for them to fix this error by simply passing a valid value. Hence, it is appropriate to raise an exception in this case. However, if there is an edge case that is unlikely to occur or cannot be easily fixed by the user, then it may be better to simply let the code fail without raising an exception.



Constants in TARDIS are taken from Astropy. The ``tardis.constants`` module
imports all constants from ``astropy.constants.astropy13constants``.

Class design and inheritance:

- If only the constructor changed, use a classmethod.
- If overriding other methods, use a subclass.

TARDIS uses Ruff to check PEP 8 compliance.

Use a classmethod when construction changes but behavior stays the same, such as
adding ``from_config(...)`` or ``from_csvy(...)`` constructors around an existing
class. Use a subclass when method behavior changes after construction, such as a
specialized solver or reader that overrides parsing or calculation methods.
