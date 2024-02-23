
from .handler.xml import get_root_of_table
from .typing import Path


class Sentinel2Datastrip:
    def __init__(self, path: Path) -> None:
        # add as optional paths of all files used here?
        self.path = path
        self.detector_time = None
        self.sat_time = None
        self.sat_ang = None
        self.sat_quat = None
        self.sat_str = None
        self.integration_time = None
        self.sampling_time = None

    def load_metadata(self) -> None:
        root = get_root_of_table()
        # read_sentinel2.read_detector_time_s2
        self.detector_time = ...
        # read_sentinel2.get_raw_str_s2
        self.sat_time = ...  # full name would make it clearer
        self.sat_ang = ...  # full name would make it clearer
        self.sat_quat = ...  # full name would make it clearer
        self.sat_str = ...  # full name would make it clearer
        # read_sentinel2.get_integration_and_sampling_time_s2
        self.integration_time = ...
        self.sampling_time = ...
