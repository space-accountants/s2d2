#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from osgeo import osr

def ll2map(ll, spatialRef):
    """ transforms angles to map coordinates (that is 2D) in a projection frame

    Parameters
    ----------
    llh : np.array, size=(m,2), unit=(deg,deg)
        np.array with spherical coordinates. In the following form:
        [[lat, lon], [lat, lon], ... ]
    spatialRef : osgeo.osr.SpatialReference
        target projection system

    Returns
    -------
    xyz : np.array, size=(m,2), unit=meter
        np.array with 2D coordinates. In the following form:
        [[x, y], [x, y], ... ]

    Examples
    --------
    Get the Universal Transverse Mercator (UTM) cooridnates from spherical
    coordinates:

    >>> import numpy as np
    >>> from osgeo import osr
    >>> proj = osr.SpatialReference()
    >>> proj.SetWellKnownGeogCS('WGS84')
    >>> lat, lon = 52.09006426183974, 5.173794246145571# Utrecht University
    >>> NH = True if lat>0 else False
    >>> proj.SetUTM(32, True)

    >>> xy = ll2map(np.array([[lat, lon]]), proj)
    >>> xy
    array([[ 237904.03625329, 5777964.65056734,       0.        ]])
    """
    if isinstance(spatialRef, str):
        spatialStr = spatialRef
        spatialRef = osr.SpatialReference()
        spatialRef.ImportFromWkt(spatialStr)
    llSpatialRef = osr.SpatialReference()
    llSpatialRef.ImportFromEPSG(4326)

    coordTrans = osr.CoordinateTransformation(llSpatialRef, spatialRef)
    xy = coordTrans.TransformPoints(list(ll))
    xy = np.stack(xy, axis=0)
    return xy

def map2ll(xy, spatialRef):
    """ transforms map coordinates (that is 2D) in a projection frame to angles

    Parameters
    ----------
    xy : np.array, size=(m,2), unit=meter
        np.array with 2D coordinates. In the following form:
        [[x, y], [x, y], ... ]
    spatialRef : osgeo.osr.SpatialReference
        target projection system

    Returns
    -------
    ll : np.array, size=(m,2), unit=(deg,deg)
        np.array with spherical coordinates. In the following form:
        [[lat, lon], [lat, lon], ... ]
    """
    if isinstance(spatialRef, str):
        spatialStr = spatialRef
        spatialRef = osr.SpatialReference()
        spatialRef.ImportFromWkt(spatialStr)
    llSpatialRef = osr.SpatialReference()
    llSpatialRef.ImportFromEPSG(4326)

    coordTrans = osr.CoordinateTransformation(spatialRef, llSpatialRef)
    ll = coordTrans.TransformPoints(list(xy))
    ll = np.stack(ll, axis=0)
    return ll[:, :-1]

def ecef2llh(xyz):
    """ transform 3D cartesian Earth Centered Earth fixed coordinates, to
    spherical angles and height above the ellipsoid

    Parameters
    ----------
    xyz : np.array, size=(m,3), unit=meter
        np.array with 3D coordinates, in WGS84. In the following form:
        [[x, y, z], [x, y, z], ... ]

    Returns
    -------
    llh : np.array, size=(m,2), unit=(deg,deg,meter)
        np.array with angles and height. In the following form:
        [[lat, lon, height], [lat, lon, height], ... ]
    """

    ecefSpatialRef = osr.SpatialReference()
    ecefSpatialRef.ImportFromEPSG(4978)

    llhSpatialRef = osr.SpatialReference()
    llhSpatialRef.ImportFromEPSG(4979)

    coordTrans = osr.CoordinateTransformation(ecefSpatialRef, llhSpatialRef)
    llh = coordTrans.TransformPoints(list(xyz))
    llh = np.stack(llh, axis=0)
    return llh

def ecef2map(xyz, spatialRef):
    """ transform 3D cartesian Earth Centered Earth fixed coordinates, to
    map coordinates (that is 2D) in a projection frame

    Parameters
    ----------
    xyz : np.array, size=(m,3), float
        np.array with 3D coordinates, in WGS84. In the following form:
        [[x, y, z], [x, y, z], ... ]
    spatialRef : osgeo.osr.SpatialReference
        target projection

    Returns
    -------
    xyz : np.array, size=(m,2), float
        np.array with planar coordinates, within a given projection frame
    """
    if isinstance(spatialRef, str):
        spatialStr = spatialRef
        spatialRef = osr.SpatialReference()
        spatialRef.ImportFromWkt(spatialStr)

    llh = ecef2llh(xyz)  # get spherical coordinates and height
    xy = ll2map(llh[:, :-1], spatialRef)
    return xy

def get_utm_zone(ϕ, λ):
    """ get the UTM zone for a specific location

    Parameters
    ----------
    ϕ_lim, λ_lim : float, unit=degrees
        latitude and longitude of a point of interest

    Returns
    -------
    utm_zone : string
        string specifying the UTM zone
    """
    ϕ_zones = [chr(i) for i in list(range(67,73)) +
               list(range(74,79)) + list(range(80,89))] # tile letters
    ϕ_cen = np.append(np.arange(-80, 72 + 1, 8), 84)
    λ_cen = np.arange(-180, 180+1, 6)

    ϕ_idx, λ_num = np.argmin(np.abs(ϕ_cen-ϕ)), np.argmin(np.abs(λ_cen-λ))
    λ_num += 1 # OBS: not a python index, but a numbering
    ϕ_idx, λ_num = np.minimum(ϕ_idx, 20), np.minimum(λ_num, 60)

    # in Southern Norway and Svalbard, the utM zones are merged
    if np.all((λ_num>31, λ_num<37, np.mod(λ_num,2)!=1, ϕ_idx==19)):
        # is location situated on Svalbard?
        if λ_num==32:
            λ_num = 31 if ϕ<9 else 33
        elif λ_num==34:
            λ_num = 33 if ϕ<21 else 35
        elif λ_num == 36:
            λ_num = 35 if ϕ<33 else 37
    elif np.all((λ_num==31, ϕ_idx==17, λ>3)):
        # is location situated in Southern Norway?
        λ_num = 32

    utm_zone = str(λ_num).zfill(2) + ϕ_zones[ϕ_idx]
    return utm_zone
