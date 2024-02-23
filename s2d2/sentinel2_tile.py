from typing import Iterable, Optional

from .handler.xml import get_root_of_table
from .typing import Path
from .sentinel2_instrument import CENTRAL_WAVELENGTH_MSI, dn_to_toa


class Sentinel2Tile:
    def __init__(self, path: Path) -> None:
        # add as optional paths of all files used here?
        self.path = path
        self.geotransform = None
        self.orbit_number = None
        self.sun_azimuth = None
        self.sun_zenith = None
        self.sun_azimuth_mean = None
        self.sun_zenith_mean = None
        self.view_azimuth = None
        self.view_zenith = None

    def load_metadata(self) -> None:
        # only metadata loading/parsing happening here
        root = get_root_of_table()
        # read_sentinel2.read_geotransform_s2
        self.geotransform = ...
        # read_sentinel2.read_orbit_number_s2
        self.orbit_number = ...
        # read_sentinel2.read_sun_angles_s2
        self.sun_azimuth = ...
        self.sun_zenith = ...
        # read_sentinel2.read_mean_sun_angles_s2
        self.sun_azimuth_mean = ...
        self.sun_zenith_mean = ...
        # read_sentinel2.read_view_angles_s2
        self.view_azimuth = ...
        self.view_zenith = ...

    def read_band(self, band: str, toa: bool = False):
        pass

    def read_bands(self, bands: Optional[Iterable] = None, toa: bool = False):
        # if "bands" not given (default), read all bands
        bands = ...
        if toa:
            bands = dn_to_toa(bands))
        return bands

    def read_detector_mask(self, bands: Optional[Iterable] = None):
        # read_sentinel2.read_detector_mask
        # if "bands" not given (default), read all bands
        pass

    def read_cloud_mask(self):
        # read_sentinel2.read_cloud_mask
        # should this be here?
        pass

    def get_sun_angle(self, angle: str, res: int = 10):
        # here is actually where the interpolation happens
        pass

    def get_view_angle(self, angle: str, bands: Optional[Iterable] = None,
                       res: int = 10):
        # here is actually where the interpolation happens
        pass

    def get_flight_bearing(self, detector_mask):
        # read_sentinel2.get_flight_bearing_from_detector_mask_s2
        # here the calculationn from the detector mask
        pass
