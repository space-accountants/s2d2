import os
import logging
import pyproj

from typing import Iterable, Optional
from osgeo import osr
from shapely import ops
from PIL import Image, ImageDraw

import xml.etree.ElementTree as ElementTree
import numpy as np
import pandas as pd

from .checking.naming import check_mgrs_code
from .handler.xml import get_root_of_table, get_branch, get_array_from_xml
from .typing import Path
from .image_coordinate_tools import map2pix, pix2map, get_max_pixel_spacing
from .sentinel2_instrument import MSI_SPECIFICS, dn_to_toa
from .sentinel2_band import Sentinel2Band
from .sentinel2_grid import Sentinel2Anglegrid
from .sentinel2_platform import S2_PLATFORM_SPECS
from .eo_imagery import bandCollection
from .orbit_tools import calculate_correct_mapping, remap_observation_angles

class Sentinel2Tile:
    def __init__(self, path: Path) -> None:
        # add as optional paths of all files used here?
        self.path = path
        self.file_dict = None
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

        self.bands = bandCollection()

    def __str__(self):
        return f"{self.path}({self.epsg})"

    def load_metadata(self) -> None:
        """
        Load basic meta-data of the tile, such as mapping specifics, as well as, rough view- and sun-angles

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

        self._update_bands_metadata()

    def _update_bands_metadata(self) -> None:
        for band_name, specs in MSI_SPECIFICS.iterrows():
            band = Sentinel2Band(band_name,
                                 self.epsg,
                                 self.geotransforms[specs.gsd],
                                 self.rows[specs.gsd],
                                 self.columns[specs.gsd])
            self.bands[band_name] = band

    def get_upperleft(self) -> list[float]:
        return list(self.geotransforms.values())[0][slice(0, -1, 3)]

    def read_bands(self,
                   bands: Optional[Iterable] = None,
                   toa: bool = False):
        # if "bands" not given (default), read all bands
        if bands is None: bands = MSI_SPECIFICS['bandid']

        for band_name, band_index in bands.items():
            self.bands[band_name].read_band(os.path.dirname(self.file_dict[band_name]),
                                            os.path.basename(self.file_dict[band_name]))
            if toa:
                self.bands[band_name].dn_to_toa()
        logging.info("Bands read")
        return

    def read_detector_masks(self,
                           bands: Optional[Iterable] = None):
        # if "bands" not given (default), read all bands
        if bands is None: bands = MSI_SPECIFICS['bandid']

        for band_name, band_index in bands.items():
            self.bands[band_name].read_detector_mask(os.path.join(self.path, 'QI_DATA'))

        logging.info("Detector masks read")
        return

    def read_cloud_mask(self):
        # read_sentinel2.read_cloud_mask
        # should this be here?
        pass

    def read_cloud_gml(self, fname='MSK_CLOUDS_B00.gml'):
        """
        The processing is performed at a spatial resolution of 60 m (the lower resolution of the three spectral bands).
        """
        f_meta = os.path.join(path_meta, 'MSK_CLOUDS_B00.gml')
        root = get_root_of_table(f_meta)

        if len(geoTransform) > 6:  # also image size is given
            msk_dim = (geoTransform[-2], geoTransform[-1])
            msk_clouds = np.zeros(msk_dim, dtype='int8')
        else:
            msk_dim = get_msk_dim_from_gml(root)
            msk_clouds = np.zeros(msk_dim, dtype='int8')  # create stack

        if len(root) > 2:  # look into meta-data for cloud polygons
            mask_members = root[2]
            for k in range(len(mask_members)):
                pos_arr = get_xy_poly_from_gml(mask_members, k)[0]

                # transform to image coordinates
                i_arr, j_arr = map2pix(geoTransform, pos_arr[:, 0], pos_arr[:, 1])
                ij_arr = np.hstack((j_arr[:, np.newaxis], i_arr[:, np.newaxis]))

                # make mask
                msk = Image.new("L", [msk_dim[1], msk_dim[0]],
                                0)  # in [width, height] format
                ImageDraw.Draw(msk).polygon(tuple(map(tuple, ij_arr[:, 0:2])),
                                            outline=1,
                                            fill=1)
                msk = np.array(msk)
                msk_clouds = np.maximum(msk_clouds, msk)
            return msk_clouds



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

    def refine_view_angles(self, chunking=False):
        # orbit_tools.calculate_correrevolutions_per_dayct_mapping
        platform = S2_PLATFORM_SPECS[self.tile_id[2]]  # 'A' or 'B'

        lat, lon, radius, inclination, period, time_para, combos = \
            calculate_correct_mapping(self.view_angle,
                                      inclination=platform.inclination,
                                      revolutions_per_day=platform.revolutions_per_day)
        logging.info("Observation angles estimated")
        self.bands = remap_observation_angles(self.view_angle, self.bands,
                                              lat, lon, radius, inclination, period, time_para, combos,
                                              chunking=chunking)
        logging.info("View angles refined")
        return

    def clip(self, polygon, epsg=4326):
        # clip to polygon and adjust the geotransform accordingly
        assert polygon.geom_type == 'Polygon', ('please provide a Polygon')

        s2_epsg = int(self.crs.GetAuthorityCode(None))
        if s2_epsg != epsg: # transform to same coordinate system
            # specify coordinate systems
            s2_proj = pyproj.CRS.from_epsg(s2_epsg)
            poly_proj = pyproj.CRS.from_epsg(epsg)
            transformer = pyproj.Transformer.from_crs(crs_from=poly_proj,
                                                      crs_to=s2_proj,
                                                      always_xy=True)
            polygon = ops.transform(transformer.transform, polygon)

        x_poly, y_poly = polygon.exterior.coords.xy

        # transform to image coordinates, take the largest image as reference
        roi = max(list(self.geotransforms.keys()))
        i_arr, j_arr = map2pix(self.geotransforms[roi], np.array(x_poly), np.array(y_poly))
        i_min, i_max = (np.maximum(np.floor(np.min(i_arr)), 0),
                        np.minimum(np.ceil(np.max(i_arr)), self.rows[roi]-1))
        j_min, j_max = (np.maximum(np.floor(np.min(j_arr)), 0),
                        np.minimum(np.ceil(np.max(j_arr)), self.columns[roi]-1))
        i_rng, j_rng = i_max-i_min, j_max-j_min
        x_min, y_max = pix2map(self.geotransforms[roi], i_min, j_min)


        #for res in geotransforms.keys():
        for key in self.geotransforms.keys():
            self.geotransforms[key] = tuple([x_min, *list(self.geotransforms[key][1:3]),
                                            y_max, *list(self.geotransforms[key][4:6])])
            self.rows[key] = int(roi/key * i_rng)
            self.columns[key] = int(roi/key * j_rng)

        for band_id, band in self.bands.items():
            # resolution of the band
            rob = get_max_pixel_spacing(band.geotransform)
            self.bands[band_id].geotransform = self.geotransforms[rob]
            self.bands[band_id].rows = self.rows[rob]
            self.bands[band_id].columns = self.columns[rob]

            if type(band.digitalnumbers) is type(None): continue

            # reduce size
            scaling = roi/rob
            row_min, row_max = int(i_min * scaling), int(i_max * scaling)
            col_min, col_max = int(j_min * scaling), int(j_max * scaling)
            img = band.digitalnumbers[row_min:row_max,col_min:col_max]

            # clip based on polygon
            i_arr, j_arr = map2pix(self.geotransforms[rob], np.array(x_poly), np.array(y_poly))
            ij_arr = np.hstack((j_arr[:, np.newaxis], i_arr[:, np.newaxis]))
            msk = Image.new("L", [self.columns[rob], self.rows[rob]], 0)

            ImageDraw.Draw(msk).polygon(tuple(map(tuple, ij_arr[:, 0:2])), outline=1,fill=1)
            msk = np.invert(np.array(msk, dtype=bool))
            img[msk] = 0

            self.bands[band_id].digitalnumbers = img

            # update detector, but do not clip to polygon
            self.bands[band_id].detector = band.detector[row_min:row_max,col_min:col_max]

        logging.info("Bands clipped")
        return

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
        if epsg is None:
            return
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
