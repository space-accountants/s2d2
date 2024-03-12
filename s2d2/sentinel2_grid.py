import numpy as np

class Sentinel2Anglegrid:
    unit: str
    epsg: str
    geotransform: list
    zentih: np.ndarray
    azimuth: np.ndarray
    band: list
    detector: list

    def __init__(self) -> None:
        # mapping specifics
        self.epsg = None

        # image specifics
        self.unit = None
        self.geotransform = [None] * 6

        self.zenith = None
        self.azimuth = None

        self.band = None
        self.detector = None

