#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import numpy as np

# geospatial libaries
from osgeo import gdal, osr

from .unit_conversion import deg2compass

def read_geo_info(fname):
    """ This function takes as input the geotiff name and the path of the
    folder that the images are stored, reads the geographic information of
    the image

    Parameters
    ----------
    fname : string
        path and file name of a geotiff image

    Returns
    -------
    spatialRef : string
        osr.SpatialReference in well known text
    geotransform : tuple, size=(8,1)
        affine transformation coefficients, but also giving the image dimensions
    targetprj : osgeo.osr.SpatialReference() object
        coordinate reference system (CRS)
    rows : integer, {x ∈ ℕ | x ≥ 0}
        number of rows in the image, that is its height
    cols : integer, {x ∈ ℕ | x ≥ 0}
        number of collumns in the image, that is its width
    bands : integer, {x ∈ ℕ | x ≥ 1}
        number of bands in the image, that is its depth

    See Also
    --------
    read_geo_image : basic function to import geographic imagery data
    """
    assert os.path.exists(fname), ('file does not seem to be present')

    img = gdal.Open(fname)
    spatialRef = img.GetProjection()
    geotransform = img.GetGeoTransform()
    geotransform = tuple(float(x) for x in geotransform)
    targetprj = osr.SpatialReference(wkt=img.GetProjection())
    rows = int(img.RasterYSize)
    cols = int(img.RasterXSize)
    bands = img.RasterCount

    geotransform += (rows, cols,)
    return spatialRef, geotransform, targetprj, rows, cols, bands

def read_geo_image(fname, boi=None, no_dat=None):
    """ This function takes as input the geotiff name and the path of the
    folder that the images are stored, reads the image and returns the data as
    an array

    Parameters
    ----------
    fname : string
        geotiff file name and path.
    boi : numpy.array, size=(k,1)
        bands of interest, if a multispectral image is read, a selection can
        be specified
    no_dat : {integer,float}
         no data value

    Returns
    -------
    data : {numpy.array, numpy.masked.array}, size=(m,n), ndim=2
        data array of the band
    spatialRef : string
        osr.SpatialReference in well known text
    geotransform : tuple, size=(8,)
        affine transformation coefficients.
    targetprj : osgeo.osr.SpatialReference() object
        coordinate reference system (CRS)

    See Also
    --------
    read_geo_info : basic function to get meta data of geographic imagery

    Examples
    --------
    >>> import os
    >>> from s2d2.mapping_input import read_geo_image
    >>> fpath = os.path.join(os.getcwd(), "data.jp2" )
    >>> I, spatialRef, geotransform, targetPrj = read_geo_image(fpath)
    """
    assert os.path.isfile(fname), ('file does not seem to be present')

    img = gdal.Open(fname)
    assert img is not None, ('could not open dataset ' + fname)

    # imagery can consist of multiple bands
    num_bands = img.RasterCount
    if boi is None: boi = np.arange(num_bands, dtype=int)

    assert (np.max(boi)+1)<=num_bands, 'bands of interest is out of range'
    for band_id, counter in enumerate(boi):
        band = np.array(img.GetRasterBand(band_id+1).ReadAsArray())
        no_dat = img.GetRasterBand(counter+1).GetNoDataValue()
        np.putmask(band, band==no_dat, np.nan)
        data = band if counter==0 else np.dstack((data,
                                                  np.atleast_3d(band)))
    spatialRef = img.GetProjection()
    geotransform = tuple(float(x) for x in img.GetGeoTransform())+data.shape[:2]
    targetprj = osr.SpatialReference(wkt=img.GetProjection())
    return data, spatialRef, geotransform, targetprj
