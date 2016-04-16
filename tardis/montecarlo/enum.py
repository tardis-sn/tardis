from ctypes import c_int


class EnumerationType(type(c_int)):
    def __new__(metacls, name, bases, dictionary):
        if not "_members_" in dictionary:
            _members_ = {}
            for key, value in dictionary.items():
                if not key.startswith("_"):
                    _members_[key] = value

            dictionary["_members_"] = _members_
        else:
            _members_ = dictionary["_members_"]

        dictionary["_reverse_map_"] = {value: key for key, value in _members_.items()}
        cls = type(c_int).__new__(metacls, name, bases, dictionary)

        for key, value in cls._members_.items():
            globals()[key] = value
        return cls

    def __repr__(self):
        return "<Enumeration %s>" % self.__name__


class CEnumeration(c_int):
    __metaclass__ = EnumerationType
    _members_ = {}

    def __eq__(self, other):
        if isinstance(other, int):
            return self.value == other

        return type(self) == type(other) and self.value == other.value

    def __repr__(self):
        value = self.value
        return "<%s.%s: %d>" % (self.__class__.__name__,
                                self._reverse_map_.get(value, '(unknown)'),
                                value)


class TardisError(CEnumeration):
    OK = 0
    BOUNDS_ERROR = 1
    COMOV_NU_LESS_THAN_NU_LINE = 2


class RPacketStatus(CEnumeration):
    IN_PROCESS = 0
    EMITTED = 1
    REABSORBED = 2


class ContinuumProcessesStatus(CEnumeration):
    OFF = 0
    ON = 1
