import numpy as np

class Sentinel2Anglegrid:
    def __init__(self) -> None:
        # mapping specifics
        self.epsg = None

        # image specifics
        self.unit = None
        self.geotransform = [None] * 6

        self.zenith = None
        self.azimuth = None

        self.band = []
        self.detector = []

