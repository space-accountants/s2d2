#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from osgeo import osr

def ll2map(ll, spatial_ref):
    """ transforms angles to map coordinates (that is 2D) in a projection frame

    Parameters
    ----------
    llh : np.array, size=(m,2), unit=(deg,deg)
        np.array with spherical coordinates. In the following form:
        [[lat, lon], [lat, lon], ... ]
    spatial_ref : osgeo.osr.SpatialReference
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
    >>> proj.SetUTM(32, lat>0)

    >>> xy = ll2map(np.array([[lat, lon]]), proj)
    >>> xy
    array([[ 237904.03625329, 5777964.65056734,       0.        ]])
    """
    if isinstance(spatial_ref, str):
        spatial_str = spatial_ref
        spatial_ref = osr.SpatialReference()
        spatial_ref.ImportFromWkt(spatial_str)
    spatial_ref_ll = osr.SpatialReference()
    spatial_ref_ll.ImportFromEPSG(4326)

    coord_trans = osr.CoordinateTransformation(spatial_ref_ll, spatial_ref)
    xy = coord_trans.TransformPoints(list(ll))
    xy = np.stack(xy, axis=0)
    return xy

def map2ll(xy, spatial_ref):
    """ transforms map coordinates (that is 2D) in a projection frame to angles

    Parameters
    ----------
    xy : np.array, size=(m,2), unit=meter
        np.array with 2D coordinates. In the following form:
        [[x, y], [x, y], ... ]
    spatial_ref : osgeo.osr.SpatialReference
        target projection system

    Returns
    -------
    ll : np.array, size=(m,2), unit=(deg,deg)
        np.array with spherical coordinates. In the following form:
        [[lat, lon], [lat, lon], ... ]
    """
    if isinstance(spatial_ref, str):
        spatial_str = spatial_ref
        spatial_ref = osr.SpatialReference()
        spatial_ref.ImportFromWkt(spatial_str)
    spatial_ref_ll = osr.SpatialReference()
    spatial_ref_ll.ImportFromEPSG(4326)

    coord_trans = osr.CoordinateTransformation(spatial_ref, spatial_ref_ll)
    ll = coord_trans.TransformPoints(list(xy))
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

    spatial_ref_ecef = osr.SpatialReference()
    spatial_ref_ecef.ImportFromEPSG(4978)

    spatial_ref_llh = osr.SpatialReference()
    spatial_ref_llh.ImportFromEPSG(4979)

    coord_trans = osr.CoordinateTransformation(spatial_ref_ecef,
                                               spatial_ref_llh)
    llh = coord_trans.TransformPoints(list(xyz))
    llh = np.stack(llh, axis=0)
    return llh

def ecef2map(xyz, spatial_ref):
    """ transform 3D cartesian Earth Centered Earth fixed coordinates, to
    map coordinates (that is 2D) in a projection frame

    Parameters
    ----------
    xyz : np.array, size=(m,3), float
        np.array with 3D coordinates, in WGS84. In the following form:
        [[x, y, z], [x, y, z], ... ]
    spatial_ref : osgeo.osr.SpatialReference
        target projection

    Returns
    -------
    xyz : np.array, size=(m,2), float
        np.array with planar coordinates, within a given projection frame
    """
    if isinstance(spatial_ref, str):
        spatial_str = spatial_ref
        spatial_ref = osr.SpatialReference()
        spatial_ref.ImportFromWkt(spatial_str)

    llh = ecef2llh(xyz)  # get spherical coordinates and height
    xy = ll2map(llh[:, :-1], spatial_ref)
    return xy

def get_utm_zone(lat, lon):
    """ get the UTM zone for a specific location

    Parameters
    ----------
    lat_lim, lon_lim : float, unit=degrees
        latitude and longitude of a point of interest

    Returns
    -------
    utm_zone : string
        string specifying the UTM zone
    """
    lat_zones = [chr(i) for i in list(range(67,73)) +
               list(range(74,79)) + list(range(80,89))] # tile letters
    lat_cen = np.append(np.arange(-80, 72 + 1, 8), 84)
    lon_cen = np.arange(-180, 180+1, 6)

    lat_idx, lon_num = np.argmin(np.abs(lat_cen-lat)), np.argmin(np.abs(lon_cen-lon))
    lon_num += 1 # OBS: not a python index, but a numbering
    lat_idx, lon_num = np.minimum(lat_idx, 20), np.minimum(lon_num, 60)

    # in Southern Norway and Svalbard, the utM zones are merged
    if np.all((lon_num>31, lon_num<37, np.mod(lon_num,2)!=1, lat_idx==19)):
        # is location situated on Svalbard?
        if lon_num==32:
            lon_num = 31 if lat<9 else 33
        elif lon_num==34:
            lon_num = 33 if lat<21 else 35
        elif lon_num == 36:
            lon_num = 35 if lat<33 else 37
    elif np.all((lon_num==31, lat_idx==17, lon>3)):
        # is location situated in Southern Norway?
        lon_num = 32

    utm_zone = str(lon_num).zfill(2) + lat_zones[lat_idx]
    return utm_zone
