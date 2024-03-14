from typing import Iterable, Optional
from osgeo import osr

import xml.etree.ElementTree as ElementTree
import numpy as np
import pandas as pd

from .checking.naming import check_mgrs_code
from .handler.xml import get_root_of_table, get_branch, get_array_from_xml
from .typing import Path
from .sentinel2_instrument import MSI_SPECIFICS, dn_to_toa
from .sentinel2_grid import Sentinel2Anglegrid

class Sentinel2Tile:
    def __init__(self, path: Path) -> None:
        # add as optional paths of all files used here?
        self.path = path
        self.tile_id = None
        self.mgrs_id = None
        self.datastrip_id = None

        # image specifics
        self.resolution = [10, 20, 60]
        self.rows = dict.fromkeys(self.resolution)
        self.columns = dict.fromkeys(self.resolution)

        # mapping specifics
        self.epsg = None
        self.geotransforms = dict.fromkeys(self.resolution, tuple([None] * 6))

        # acquisition and solar angles
        self.sun_angle = Sentinel2Anglegrid()
        self.view_angle = Sentinel2Anglegrid()
        self.sun_azimuth_mean = None
        self.sun_zenith_mean = None

        self.sensing_time = None

    def __str__(self):
        return f"{self.path}({self.epsg})"

    def load_metadata(self) -> None:
        """

        Notes
        -----
        The metadata structure of the xml-file is as follows:

        .. code-block:: text

            * MTD_TL.xml
            └ n1:Level-1C_Tile_ID
               ├ n1:General_Info
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

        # only metadata loading/parsing happening here
        root = get_root_of_table(self.path, fname='MTD_TL.xml')

        gen_info = get_branch(root, 'General_Info')
        self._get_tile_id_from_xmltree(gen_info)

        geom_info = get_branch(root, 'Geometric_Info')
        geocoding = get_branch(geom_info, 'Tile_Geocoding')

        self._get_crs_s2_from_xmltree(geocoding)
        self._get_image_dimensions_from_xmltree(geocoding)
        self._get_geotransforms_from_xmltree(geocoding)

        tile_ang = get_branch(geom_info, 'Tile_Angles')
        self._get_mean_sunangles_from_xmltree(tile_ang)
        self._get_sunangles_from_xmltree(tile_ang)
        self._get_viewangles_from_xmltree(tile_ang)

    def get_upperleft(self) -> list[float]:
        return list(self.geotransforms.values())[0][slice(0, -1, 3)]

    def read_band(self, band: str, toa: bool = False):
        pass

    def read_bands(self, bands: Optional[Iterable] = None, toa: bool = False):
        # if "bands" not given (default), read all bands
        bands = ...
        if toa:
            bands = dn_to_toa(bands)
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

    def get_utmzone_from_tile_code(self):
        """

        Returns
        -------
        self.utmzone : integer
            code used to denote the number of the projection column of UTM

        See Also
        --------
        .get_epsg_from_mgrs_tile, .get_crs_from_mgrs_tile

        Notes
        -----
        The tile structure is a follows "AABCC"
            * "AA" utm zone number, starting from the East, with steps of 8 degrees
            * "B" latitude zone, starting from the South, with steps of 6 degrees
        """
        tile_code = check_mgrs_code(tile_code)
        return int(tile_code[:2])

    def get_epsg_from_mgrs_tile(self):
        """

        Returns
        -------
        self.epsg : integer
            code used to denote a certain database entry

        See Also
        --------
        get_utmzone_from_mgrs_tile
        get_crs_from_mgrs_tile

        Notes
        -----
        The tile structure is a follows "AABCC"
            * "AA" utm zone number, starting from the East, with steps of 8 degrees
            * "B" latitude zone, starting from the South, with steps of 6 degrees
        """
        tile_code = check_mgrs_code(tile_code)
        self.get_utmzone_from_tile_code(tile_code)
        epsg_code = 32600 + utm_num

        # N to X are in the Northern hemisphere
        if tile_code[2] < 'N': epsg_code += 100
        return epsg_code

    def get_crs_from_mgrs_tile(tile_code):
        """

        Parameters
        ----------
        tile_code : string
            US Military Grid Reference System (MGRS) tile code

        Returns
        -------
        crs : osgeo.osr.SpatialReference
            target projection system

        See Also
        --------
        .get_utmzone_from_mgrs_tile, .get_utmzone_from_mgrs_tile

        Notes
        -----
        The tile structure is a follows "AABCC"
            * "AA" utm zone number, starting from the East, with steps of 8 degrees
            * "B" latitude zone, starting from the South, with steps of 6 degrees
        """
        tile_code = check_mgrs_code(tile_code)
        epsg_code = get_epsg_from_mgrs_tile(tile_code)

        crs = osr.SpatialReference()
        crs.ImportFromEPSG(epsg_code)
        return crs

    def _get_tile_id_from_xmltree(self,
                                 general_info: ElementTree.Element) -> None:
        for field in general_info:
            if field.tag == 'TILE_ID':
                self.tile_id = field.text
            elif field.tag == 'DATASTRIP_ID':
                self.datastrip_id = field.text
            elif field.tag == 'SENSING_TIME':
                self.sensing_time = pd.Timestamp(field.text)
        self.mgrs_id = self.tile_id.split('_')[-2]

    def _get_crs_s2_from_xmltree(self,
                                 geocoding: ElementTree.Element) -> None:
        epsg = None
        for field in geocoding:
            if field.tag == 'HORIZONTAL_CS_CODE':
                epsg = int(field.text.split(':')[1])
        if epsg is None: return
        crs = osr.SpatialReference()
        crs.ImportFromEPSG(epsg)
        self.epsg = epsg
        self.crs = crs


    def _get_image_dimensions_from_xmltree(self,
                                           geocoding: ElementTree.Element) -> None:
        for box in geocoding:
            if not (box.tag == 'Size'): continue
            for field in box:
                if field.tag == 'NROWS':
                    self.rows[int(box.attrib['resolution'])] = int(field.text)
                elif field.tag == 'NCOLS':
                    self.columns[int(box.attrib['resolution'])] = int(field.text)


    def _get_geotransforms_from_xmltree(self,
                                        geocoding: ElementTree.Element) -> None:

        for box in geocoding:
            if not (box.tag == 'Geoposition'): continue
            gt = self.geotransforms[int(box.attrib['resolution'])]
            for field in box:
                if field.tag == 'ULX':
                    gt = tuple([float(field.text) if idx == 0 else val for idx, val in enumerate(gt)])
                elif field.tag == 'XDIM':
                    gt = tuple([float(field.text) if idx == 1 else val for idx, val in enumerate(gt)])
                    gt = tuple([0. if idx == 2 else val for idx, val in enumerate(gt)])
                elif field.tag == 'ULY':
                    gt = tuple([float(field.text) if idx == 3 else val for idx, val in enumerate(gt)])
                elif field.tag == 'YDIM':
                    gt = tuple([0. if idx == 4 else val for idx, val in enumerate(gt)])
                    gt = tuple([float(field.text) if idx == 5 else val for idx, val in enumerate(gt)])
            self.geotransforms = {**self.geotransforms,
                                  int(box.attrib['resolution']): gt
                                  }


    def _get_sunangles_from_xmltree(self,
                                    tileangles: ElementTree.Element) -> None:
        for grids in tileangles:
            if not (grids.tag == 'Sun_Angles_Grid'): continue

            angles, col_step, row_step = None, None, None
            for instance in grids:
                for field in instance:
                    if field.tag == 'Values_List':
                        angle = get_array_from_xml(field)
                    elif field.tag == 'COL_STEP':
                        col_step = float(field.text)
                    elif field.tag == 'ROW_STEP':
                        row_step = float(field.text)
                # update grids
                self.sun_angle.add_raster_layer(instance.tag.lower(), angles)
            ul = self.get_upperleft()
            self.sun_angle.geotransform = (ul[0], col_step, 0., ul[1], 0., -1*row_step)
            self.sun_angle.unit = 'deg'
            self.sun_angle.epsg = self.epsg


    def _get_mean_sunangles_from_xmltree(self,
                                         tileangles: ElementTree.Element) -> None:
        for grids in tileangles:
            if not (grids.tag == 'Mean_Sun_Angle'): continue

            for instance in grids:
                if instance.tag == 'ZENITH_ANGLE':
                    self.sun_zenith = float(instance.text)
                elif instance.tag == 'AZIMUTH_ANGLE':
                    self.sun_azimuth = float(instance.text)

    def _get_viewangles_from_xmltree(self,
                                         tileangles: ElementTree.Element) -> None:
        for grids in tileangles:
            if not (grids.tag == 'Viewing_Incidence_Angles_Grids'): continue

            self.view_angle.band.append(int(grids.attrib['bandId']))
            self.view_angle.detector.append(int(grids.attrib['detectorId']))

            angles, col_step, row_step = None, None, None
            for instance in grids:
                for field in instance:
                    if field.tag == 'Values_List':
                        angles = get_array_from_xml(field)
                    elif field.tag == 'COL_STEP':
                        col_step = float(field.text)
                    elif field.tag == 'ROW_STEP':
                        row_step = float(field.text)
                # update grids
                self.view_angle.add_raster_layer(instance.tag.lower(), angles)
            ul = self.get_upperleft()
            self.view_angle.geotransform = (ul[0], col_step, 0., ul[1], 0., -1 * row_step)
            self.view_angle.unit = 'deg'
            self.view_angle.epsg = self.epsg
