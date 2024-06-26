#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from s2d2.unit_conversion import deg2arg

def zenit_angle_check(theta):
    """ check if the solar zenit angle is within range

    Parameters
    ----------
    theta : np.ndarray, dtype=float, unit=degrees, range=-90...+90
        solar zenith angle

    Notes
    -----
    The angles related to the sun are as follows:

        .. code-block:: text

                *                     * sun
          ^    /                ^    /|
          |   /                 |   / | nadir
          |-- zenith angle      |  /  v
          | /                   | /|
          |/                    |/ | elevation angle
          └----- surface        └------

    """
    theta = np.maximum(theta, 0.)
    theta = np.minimum(theta, 90.)
    return theta

def lat_lon_angle_check(lat,lon):
    """

    Parameters
    ----------
    lat : np.ndarray, dtype=float, unit=degrees, range=-90...+90
        latitude
    lon : float, unit=degrees, unit=degrees, range=-180...+180
        longitude

    Returns
    -------
    lat,lon : float
        latitude and longitude

    See Also
    --------
    .deg2arg : transform angle to range of -180...+180
    """
    lat = np.maximum(lat, -90.)
    lat = np.minimum(lat, +90.)

    lon = deg2arg(lon)
    return lat, lon

def is_crs_an_srs(crs):
    """ is the coordinate reference system given in degrees, or a spatial
    reference system (with a metric scale).

    Parameters
    ----------
    crs : string
        osr.SpatialReference in well known text

    Returns
    -------
    verdict : bool
        is the reference system a mapping system

    See Also
    --------
    s2d2.image_coordinate_tools.create_local_crs
    """
    if not isinstance(crs, str): crs = crs.ExportToWkt()
    return crs.find('"metre"')!=-1

def correct_geotransform(geoTransform):
    """ check if image information in the form of its size is also provided

    Parameters
    ----------
    geoTransform : tuple, size={(6,), (8,)}
        affine transformation coefficients.

    Returns
    -------
    geoTransform : tuple, size=(8,)
        affine transformation coefficients.
    """
    assert isinstance(geoTransform, tuple), 'geoTransform should be a tuple'
    assert len(geoTransform) == 6, 'geoTransform not of the correct size'
    assert all(isinstance(n, float) for n in geoTransform)
    return geoTransform
