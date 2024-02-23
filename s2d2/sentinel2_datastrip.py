
import numpy as np

from .handler.xml import get_root_of_table
from .mapping_tools import ecef2llh
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
        # also in read_sentinel2.get_flight_orientation_s2 ?
        self.sat_time = ...  # full name would make it clearer
        self.sat_ang = ...  # full name would make it clearer
        self.sat_quat = ...  # full name would make it clearer
        self.sat_str = ...  # full name would make it clearer
        # read_sentinel2.get_integration_and_sampling_time_s2
        self.integration_time = ...
        self.sampling_time = ...
        # sensor_readings_sentinel2.get_flight_path_s2
        self.sat_time = ... # full name would make it clearer - is it same as from read_sentinel2.get_raw_str_s2?
        self.sat_xyz = ... # full name would make it clearer
        self.sat_err = ... # full name would make it clearer
        self.sat_uvw = ... # full name would make it clearer

    def get_altitude(self):
        # sensor_readings_sentinel2.get_flight_path_s2
        # estimate the altitude above the ellipsoid
        llh = ecef2llh(self.sat_xyz)
        return np.squeeze(llh[:,-1])

    def get_velocity(self):
        # sensor_readings_sentinel2.get_flight_path_s2
        # estimate platform speed
        velo = np.linalg.norm(self.sat_uvw, axis=1)
        return np.squeeze(velo)

    #TODO def get_intrinsic_temperatures_s2(self):