import importlib.util

import pkg_resources


class ModuleChecker:

    def find_module(self, module):
        """Search for modules specification."""
        try:
            return importlib.util.find_spec(module)
        except ImportError:
            return None

    def find_distribution(self, dist):
        """Search for distribution with specified version (eg 'numpy>=1.15')."""
        try:
            return pkg_resources.require(dist)
        except Exception:
            return None

    def check(self, module):
        """
        Return True if module with specified version exists.
        >>> ModuleChecker().check('foo>=1.0')
        False
        >>> ModuleChecker().check('pytest>1.0')
        True
        """
        mods = self.find_module(module) or self.find_distribution(module)
        return bool(mods)
