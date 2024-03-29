import numpy as np

from .checking.array import are_two_arrays_equal
from .checking.mapping import correct_geotransform


class Sentinel2Anglegrid:
    def __init__(self) -> None:
        # mapping specifics
        self.epsg = None

        # image specifics
        self.unit = None
        self.geotransform = tuple([None] * 6)

        self.zenith = None
        self.azimuth = None

        self.band = []
        self.detector = []

        # grid dimensions
        self.rows = None
        self.columns = None
        self.depth = None

    def add_raster_layer(self,
                         angle_type: str,
                         angles: np.ndarray):
        grid = getattr(self, angle_type)
        if grid is None:
            grid = np.atleast_3d(angles)
            self.rows = grid.shape[0]
            self.columns = grid.shape[1]
        else:
            grid = np.dstack((grid, angles))
        self.depth = grid.shape[2]
        setattr(self, angle_type, grid)

    def check_consistency(self) -> None:
        are_two_arrays_equal(self.zenith, self.azimuth)
        correct_geotransform(self.geotransform)
        if self.band is not None:
            are_two_arrays_equal(self.band, self.detector)
        assert self.rows == self.zenith.shape[0], 'dimensions are not consistent with data'
        assert self.columns == self.zenith.shape[1], 'dimensions are not consistent with data'
        assert self.depth == self.zenith.shape[2], 'dimensions are not consistent with data'
