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
        self.path = None

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
        """ Reads the image into the Sentinel2Band class as a masked array

        Parameters
        ----------
        path : str
            Folder where the data is located
        fname : str
            Name of the image
        """
        if fname.find('.') == -1: # no extension is given
            img_formats = ['.jp2', '.tif', '.tiff']
            for suffix in img_formats:
                if os.path.exists(os.path.join(path, fname+suffix)):
                    fname += suffix
                    continue

        img_path = os.path.join(path, fname)
        img = gdal.Open(img_path)
        assert img is not None, ('could not open dataset ' + fname)

        #todo: get satval and nanval from Sentinel2Product

        band = np.array(img.GetRasterBand(1).ReadAsArray())
        no_dat = img.GetRasterBand(1).GetNoDataValue()
        np.putmask(band, band == no_dat, 0)

        self.digitalnumbers = band
        self.unit = 'DN'

    def dn_to_toa(self):
        """
        convert the digital numbers of Sentinel-2 data to top of atmosphere (TOA), see
        for more details [wwwS2L1C]_.

        Notes
        -----
        .. [wwwS2L1C] https://sentinel.esa.int/web/sentinel/technical-guides/sentinel-2-msi/level-1c/algorithm
        """
        if self.unit == 'TOA':
            return
        assert (self.unit == 'DN'), 'input should be digital numbers'

        self.digitalnumbers = dn_to_toa(self.digitalnumbers)
        self.unit = 'TOA'
