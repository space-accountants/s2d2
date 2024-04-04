from collections import UserDict
from .sentinel2_band import Sentinel2Band

class bandCollection(UserDict):
    def __setitem__(self, key, value):
        key = key
        super().__setitem__(key, value)

    def __getitem__(self, key):
        key = str(key)
        return self.data[key]

class spectralIndices:
    pass
