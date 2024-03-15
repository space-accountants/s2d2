from collections import namedtuple
from PIL import Image, ImageDraw

import os
import numpy as np

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

        self.units = None
        self.digitalnumbers = None
        self.detector = None
        self.zenith = None
        self.azimuth = None
        self.timing = None

    def read_detector_gml(self, path: Path):
        det_msk = np.zeros((self.rows, self.columns), dtype='int8')

        root = get_root_of_table(self.path, fname=f'MSK_DETFOO_{self.index}.gml')
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
        setattr(self, 'detector', det_msk)

    def read_band(self, path: Path):

        setattr(self, 'unit', 'DN')

    def dn_to_toa(self):
        if getattr(self, 'unit') == 'TOA': return
        assert(getattr(self, 'unit') == 'DN'), 'input should be digital numbers'

        TOA = dn_to_toa(getattr(self, 'digitalnumbers'))

        setattr(self, 'digitalnumbers', TOA)
        setattr(self, 'unit', 'TOA')
