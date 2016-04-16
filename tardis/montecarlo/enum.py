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
            _members_ = dict["_members_"]

        dictionary["_reverse_map_"] = {value: key for key, value in _members_.items()}
        cls = type(c_int).__new__(metacls, name, bases, dict)

        for key, value in cls._members_.items():
            globals()[key] = value
        return cls


class CEnumeration(c_int):
    __metaclass__ = EnumerationType
    _members_ = {}

    def __eq__(self, other):
        if isinstance(other, int):
            return self.value == other

        return type(self) == type(other) and self.value == other.value
