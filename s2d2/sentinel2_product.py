import os
import numpy as np

from typing import Optional

from .handler.xml import get_root_of_table, get_branch
from .orbit_tools import calculate_correct_mapping
from .sentinel2_datastrip import Sentinel2Datastrip
from .sentinel2_tile import Sentinel2Tile
from .sentinel2_instrument import S2A_PLATFORM_SPECS, S2B_PLATFORM_SPECS
from .typing import Path


class Sentinel2Product:
    """

    Notes
    -----
    For Sentinel-2 Level 1C the metadata is scattered over the files and folders.
    In order to make this organization more understandable, the following visualization is made:

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

    datastrip: Sentinel2Datastrip
    tile: Sentinel2Tile
    path: Path
    rel_img_dir: Path
    rel_ds_dir: Path
    spacecraft: str
    nanval: int
    satval: int

    def __init__(self, path: Optional[Path] = None) -> None:
        # add as optional paths of all files used here?
        # default calculate relative paths for datastrip and granule
        # and pass them to objects below
        self.path = path
        self.datastrip = Sentinel2Datastrip(path)
        self.tile = Sentinel2Tile(path)
        self.sensing_time = None
        self.spacecraft = None
        self.nanval = None
        self.satval = None
        self.rel_img_dir = None
        self.rel_ds_dir = None

    def load_metadata(self) -> None:
        """ load meta-data from the product file (MTD_TL.xml)

        Notes
        -----
        The metadata structure of the Sentinel2 product is as follows:

        .. code-block:: text

            * MTD_MSIL1C.xml
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

        assert os.path.exists(self.path), 'folder does not seem to exist'

        root = get_root_of_table(self.path, fname='MTD_MSIL1C.xml')

        gnrl_info = get_branch(root, 'General_Info')
        prod_info = get_branch(gnrl_info, 'Product_Info')
        data_take = get_branch(prod_info, 'Datatake')
        self._get_spacecraft_from_xmlstruct(data_take)

        imag_spc = get_branch(gnrl_info, 'Product_Image_Characteristics')
        self._get_special_image_values(imag_spc)
        imag_org = get_branch(prod_info, 'Product_Organisation')
        gran_org = get_branch(imag_org, 'Granule_List')

        # get location where imagery is situated
        self._get_image_dir(gran_org)
        self.tile.path = os.path.join(self.path, os.path.dirname(self.rel_img_dir))
        self.tile.load_metadata()

        # read_sentinel2.read_sensing_time_s2
        self.sensing_time = ...


        self.datastrip.load_metadata()

    def _get_spacecraft_from_xmlstruct(self, data_take) -> None:
        platform = None
        for field in data_take:
            if field.tag == 'SPACECRAFT_NAME':
                platform = field.text[-1].upper()
        if platform is None: return
        self.spacecraft = platform

    def _get_special_image_values(self, imag_spc) -> None:
        for spec in imag_spc:
            if spec.tag == 'Special_Values':
                if spec[0].text == 'NODATA':
                    self.nanval = int(spec[1].text)
                elif spec[0].text == 'SATURATED':
                    self.satval =  int(spec[1].text)

    def _get_image_dir(self, gran_org) -> None:
        rel_path = os.path.dirname(gran_org[0][0].text)
        self.rel_img_dir = rel_path


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
