import numpy as np

from typing import Optional

from .handler.xml import get_root_of_table
from .orbit_tools import calculate_correct_mapping
from .sentinel2_datastrip import Sentinel2Datastrip
from .sentinel2_tile import Sentinel2Tile
from .sentinel2_instrument import S2A_PLATFORM_SPECS, S2B_PLATFORM_SPECS
from .typing import Path


class Sentinel2Product:
    """

    Notes
    -----
    The metadata is scattered over the file structure of Sentinel-2, Level
    1C meta data files. The folders and files are as follows:

    .. code-block:: text

        * S2X_MSIL1C_20XX...
        ├ AUX_DATA
        ├ DATASTRIP
        │  └ DS_XXX_XXXX...
        │     └ QI_DATA
        │        └ MTD_DS.xml <- metadata about the data-strip
        ├ GRANULE
        │  └ L1C_TXXXX_XXXX...
        │     ├ AUX_DATA
        │     ├ IMG_DATA
        │     ├ QI_DATA
        │     └ MTD_TL.xml <- metadata about the tile
        ├ HTML
        ├ rep_info
        ├ manifest.safe
        ├ INSPIRE.xml
        └ MTD_MSIL1C.xml <- metadata about the product
    """

    def __init__(self, path: Optional[Path] = None) -> None:
        # add as optional paths of all files used here?
        # default calculate relative paths for datastrip and granule
        # and pass them to objects below
        self.path = path
        self.datastrip = Sentinel2Datastrip(path)
        self.tile = Sentinel2Tile(path)
        self.sensing_time = None
        self.spacecraft = None

    def load_metadata(self) -> None:
        """ load meta-data from the product file (MTD_TL.xml)

        Notes
        -----
        The metadata structure of the Sentinel2 product is as follows:

        .. code-block:: text

            * MTD_TL.xml
            └ n1:Level-1C_Tile_ID
               ├ n1:General_Info
               │  ├ Product_Info
               │  │  ├ PRODUCT_URI <- .SAFE folder name
               │  │  ├ PROCESSING_LEVEL
               │  │  └ PROCESSING_BASELINE
               │  ├ SPACECRAFT_NAME
               │  ├ SENSING_ORBIT_NUMBER
               │  └ SENSING_ORBIT_DIRECTION
               ├ n1:Geometric_Info
               │  ├ Tile_Geocoding
               │  │  ├ HORIZONTAL_CS_NAME
               │  │  ├ HORIZONTAL_CS_CODE
               │  │  ├ Size : resolution={"10","20","60"}
               │  │  │  ├ NROWS
               │  │  │  └ NCOLS
               │  │  └ Geoposition
               │  │     ├ ULX
               │  │     ├ ULY
               │  │     ├ XDIM
               │  │     └ YDIM
               │  └ Tile_Angles
               └ n1:Quality_Indicators_Info
        """
        root = get_root_of_table()
        # read_sentinel2.read_sensing_time_s2
        self.sensing_time = ...
        self.spacecraft = 'A' # _get_spacecraft_s2_from_root(root)

        self.tile.load_metadata()
        self.datastrip.load_metadata()

    def get_flight_bearing_from_gnss(self) -> np.ndarray:
        # sensor_readings_sentinel2.get_flight_bearing_from_gnss_s2
        # is the name self-explanatory?
        pass

    def calculate_correct_mapping(self):
        # orbit_tools.calculate_correct_mapping

        if self.spacecraft == 'A':
            inclination = S2A_PLATFORM_SPECS['inclination']
            revolutions_per_day = S2A_PLATFORM_SPECS['revolutions_per_day']
        else:
            inclination = S2B_PLATFORM_SPECS['inclination']
            revolutions_per_day = S2B_PLATFORM_SPECS['revolutions_per_day']

        lat, lon, radius, inclination, period, time_para, combos = \
            calculate_correct_mapping(zn_grd, az_grd, bnd, det, grdtransform, crs,
                                  inclination, revolutions_per_day)

        pass
