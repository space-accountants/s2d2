#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Space Accountants"
__license__ = "MIT License - You must cite this source"
__version__ = "202311"
__maintainer__ = "B. Altena"
__email__ = "info at space hyphen accountants dot eu"

import numpy as np

from s2d2.unit_conversion import deg2arg

def zenit_angle_check(zn):
    zn = np.maximum(zn, 0.)
    zn = np.minimum(zn, 90.)
    return zn

def lat_lon_angle_check(ϕ,λ):
    """

    Parameters
    ----------
    ϕ : float, unit=degrees, range=-90...+90
        latitude
    λ : float, unit=degrees
        longitude

    Returns
    -------
    ϕ : float
        latitude
    λ : float
        longitude

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
    dhdt.testing.mapping_tools.create_local_crs
    """
    if not isinstance(crs, str): crs = crs.ExportToWkt()
    return crs.find('"metre"')!=-1

def correct_geoTransform(geoTransform):
    assert isinstance(geoTransform, tuple), 'geoTransform should be a tuple'
    assert len(geoTransform) in (6,8,), 'geoTransform not of the correct size'
    if len(geoTransform)==6:
        assert all(isinstance(n, float) for n in geoTransform)
    else:
        assert all(isinstance(n, float) for n in geoTransform[:-2])
        assert all(isinstance(n, int) for n in geoTransform[-2:])
        assert all(n>=0 for n in geoTransform[-2:])
    return geoTransform
