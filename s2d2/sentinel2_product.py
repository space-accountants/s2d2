from typing import Optional

from .handler.xml import get_root_of_table
from .sentinel2_datastrip import Sentinel2Datastrip
from .sentinel2_tile import Sentinel2Tile
from .typing import Path


class Sentinel2Product:
    def __init__(self, path: Optional[Path] = None) -> None:
        # add as optional paths of all files used here?
        # default calculate relative paths for datastrip and granule
        # and pass them to objects below
        self.path = path
        self.datastrip = Sentinel2Datastrip(path)
        self.tile = Sentinel2Tile(path)
        self.sensing_time = None

    def load_metadata(self) -> None:
        root = get_root_of_table()
        # read_sentinel2.read_sensing_time_s2
        self.sensing_time = ...

        self.tile.load_metadata()
        self.datastrip.load_metadata()

    def calculate_correct_mapping(self):
        # orbit_tools.calculate_correct_mapping
        pass
