from collections import namedtuple
from PIL import Image, ImageDraw

import os
import numpy as np

from osgeo import gdal

from typing import Optional
from .typing import Path

from .image_coordinate_tools import map2pix
from .handler.xml import get_root_of_table, get_branch
from .handler.gml import get_xy_polygon_from_gml

from .sentinel2_instrument import MSI_SPECIFICS, dn_to_toa


class Sentinel2Band():
    def __init__(self,
                 index: str,
                 epsg: int,
                 geotransform: tuple,
                 rows: int,
                 columns: int) -> None:
        self.index = index

        # mapping specifics
        self.epsg = epsg
        self.geotransform = geotransform

        # image specifics
        self.rows = rows
        self.columns = columns

        self.unit = None
        self.digitalnumbers = None
        self.detector = None
        self.zenith = None
        self.azimuth = None
        self.timing = None

    def read_detector_gml(self, path: Path):
        det_msk = np.zeros((self.rows, self.columns), dtype='int8')

        root = get_root_of_table(path, fname=f'MSK_DETFOO_{self.index}.gml')
        mask_members = get_branch(root, 'maskMembers')
        for k in range(len(mask_members)):
            pos_arr, det_num = get_xy_polygon_from_gml(mask_members, k)

            # transform to image coordinates
            i_arr, j_arr = map2pix(self.geotransform, pos_arr[:, 0], pos_arr[:, 1])
            ij_arr = np.hstack((j_arr[:, np.newaxis], i_arr[:, np.newaxis]))
            # make mask
            msk = Image.new("L", [np.size(det_msk, 1), np.size(det_msk, 0)], 0)
            ImageDraw.Draw(msk).polygon(tuple(map(tuple, ij_arr[:, 0:2])),
                                        outline=det_num, fill=det_num)
            msk = np.array(msk)
            det_msk = np.maximum(det_msk, msk)
        self.detector = det_msk

    def read_band(self,
                  path: Path,
                  fname: Optional[str]):
        # todo: what to do with full integration...
        img_path = os.path.join(path, fname)
        img = gdal.Open(img_path)
        assert img is not None, ('could not open dataset ' + fname)

        band = np.array(img.GetRasterBand(1).ReadAsArray())
        no_dat = img.GetRasterBand(1).GetNoDataValue()
        np.putmask(band, band == no_dat, 0)

        self.digitalnumbers = band
        self.unit = 'DN'

    def dn_to_toa(self):
        if self.unit == 'TOA':
            return
        assert (self.unit == 'DN'), 'input should be digital numbers'

        self.digitalnumbers = dn_to_toa(self.digitalnumbers)
        self.unit = 'TOA'
