#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from s2d2.unit_conversion import deg2arg

def zenit_angle_check(θ):
    """ check if the solar zenit angle is within range

    Parameters
    ----------
    θ : np.ndarray, dtype=float, unit=degrees, range=-90...+90
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
    θ = np.maximum(θ, 0.)
    θ = np.minimum(θ, 90.)
    return θ

def lat_lon_angle_check(ϕ,λ):
    """

    Parameters
    ----------
    ϕ : np.ndarray, dtype=float, unit=degrees, range=-90...+90
        latitude
    λ : float, unit=degrees, unit=degrees, range=-180...+180
        longitude

    Returns
    -------
    ϕ,λ : float
        latitude and longitude

    See Also
    --------
    .deg2arg : transform angle to range of -180...+180
    """
    ϕ = np.maximum(ϕ, -90.)
    ϕ = np.minimum(ϕ, +90.)

    λ = deg2arg(λ)
    return ϕ, λ

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

def correct_geoTransform(geoTransform):
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
    assert len(geoTransform) in (6,8,), 'geoTransform not of the correct size'
    if len(geoTransform)==6:
        assert all(isinstance(n, float) for n in geoTransform)
    else:
        assert all(isinstance(n, float) for n in geoTransform[:-2])
        assert all(isinstance(n, int) for n in geoTransform[-2:])
        assert all(n>=0 for n in geoTransform[-2:])
    return geoTransform
