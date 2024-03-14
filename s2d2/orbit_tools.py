#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from osgeo import osr

from typing import Optional

import numpy as np
import pandas as pd

from .image_coordinate_tools import pix_centers
from .mapping_tools import map2ll, ecef2llh
from .checking.mapping import (
    correct_geotransform, lat_lon_angle_check, is_crs_an_srs)
from .checking.array import (
    are_two_arrays_equal, are_three_arrays_equal)
from .sentinel2_grid import Sentinel2Anglegrid


def wgs84_param() -> tuple[float, float]:
    """ get paramters of the WGS84 ellipsoid

    Returns
    -------
    major_axis : float, unit=meter
        largest axis of the ellipsoid
    flattening : float
        amount of compression of the ellipsoid

    See Also
    --------
    earth_eccentricity, earth_axes

    Notes
    -----
    The flattening (f) stem from the major- (a) en minor-axis (b) via:

    .. math:: f = (a - b)/a

    or via the eccentricity (e):

    .. math:: f = 1 - \sqrt{ 1 - e^2}
    """
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    major_axis = wgs84.GetSemiMajor()
    flattening = 1 / wgs84.GetInvFlattening()
    return major_axis, flattening

def earth_eccentricity() -> float:
    """ get the eccentricity of the WGS84 ellipsoid

    Returns
    -------
    eccentricity : float
        amount of deviation from circularity

    See Also
    --------
    wgs84_param, earth_axes
    """
    major_axis, flattening = wgs84_param()
    eccentricity = (2 * flattening) - (flattening**2)
    return eccentricity

def earth_axes() -> tuple[float, float]:
    """ get axis length of the WGS84 ellipsoid

    Returns
    -------
    major_axis : float, unit=meter
        the largest axis of the ellipsoid
    minor_axis : float, unit=meter
        the shortest axis of the ellipsoid

    See Also
    --------
    earth_eccentricity, wgs84_param
    """
    major_axis, flattening = wgs84_param()
    eccentricity = earth_eccentricity()
    minor_axis = major_axis * np.sqrt(1 - eccentricity)
    return major_axis, minor_axis

def standard_gravity() -> float:
    """
    provide value for the standard gravitational parameter of the Earth

    Returns
    -------
    mu : float, unit= m**3 * s**-2
        standard gravity
    """
    mu = 3.986004418E14
    return mu

def transform_rpy_xyz2ocs(xyz, uvw, roll, pitch, yaw, xyz_time, ang_time):
    x,y,z = np.interp(ang_time, xyz_time, xyz[...,0]), \
            np.interp(ang_time, xyz_time, xyz[...,1]), \
            np.interp(ang_time, xyz_time, xyz[...,2])
    u,v,w = np.interp(ang_time, xyz_time, uvw[...,0]), \
            np.interp(ang_time, xyz_time, uvw[...,1]), \
            np.interp(ang_time, xyz_time, uvw[...,2])
    xyz, uvw = np.stack((x,y,z), axis=1), np.stack((u,v,w), axis=1)

    m,n = xyz.shape
    z_ocs = np.divide(xyz, np.tile(np.linalg.norm(xyz, axis=-1), (n,1)).T)
    x_ocs = np.cross(uvw,xyz, axisa=-1, axisb=-1)
    x_ocs = np.divide(x_ocs, np.tile(np.linalg.norm(x_ocs, axis=-1), (n, 1)).T)
    y_ocs = np.cross(z_ocs, x_ocs, axisa=-1, axisb=-1)
    y_ocs = np.divide(y_ocs, np.tile(np.linalg.norm(y_ocs, axis=-1), (n, 1)).T)

    r_stack = np.dstack((x_ocs,y_ocs,z_ocs))
    rpy = np.stack((roll, pitch, yaw), axis=1)

    np.einsum('...i,...i->...i', r_stack, rpy)
    return roll, pitch, yaw

def estimate_inclination_via_xyz_uvw(xyz, uvw):
    xyz, uvw = np.atleast_2d(xyz), np.atleast_2d(uvw)
    m,n = xyz.shape
    h = np.cross(xyz, uvw)
    h_mag = np.linalg.norm(h, axis=-1)
    h = np.divide(h, np.tile(h_mag, (n,1)).T)

    z = np.array([0,0,1])
    i = 90 - np.rad2deg(np.arccos(np.dot(h, z)))
    i += 90
    return i

def calculate_correct_mapping(grid: Sentinel2Anglegrid,
                              inclination:float = 98.5621,
                              revolutions_per_day:float = 14.30824258387262,
                              radius: Optional[float] = None,
                              mean_altitude: Optional[float] = None):
    """

    Parameters
    ----------
    grid : Sentinel2Anglegrid
        with the following entries:
            - zenith : {numpy.ndarray, numpy.masked.array}, size=(k,l,h)
                observation angles of the different detectors/bands
            - azimuth : {numpy.ndarray, numpy.masked.array}, size=(k,l,h)
                observation angles of the different detectors/bands
            - band : numpy.ndarray, size=(h,)
                number of the band, corresponding to the third dimension of 'zenith'
            - detector : numpy.ndarray, size=(h,)
                number of the detector, corresponding to the third dimension of 'zenith'
            - geotransform : tuple, size=(6,)
                geotransform of the grid of 'zenit' and 'azimuth'
            - crs : osgeo.osr.SpatialReference() object
                coordinate reference system (CRS)
    inclination : float, unit=degrees
        angle of the orbital plane in relation to the equatorial plane
    revolutions_per_day : float
        amount of revolutions a satellite platform performs around the Earth
    mean_altitude : float, unit=meter
        estimate of mean satellite altitude above the Earth surface

    Returns
    -------
    l_time : numpy.ndarray, size=(p,), unit=seconds
        asd
    lat, lon : float, unit=degrees
        ground location of the satelite at time 0
    radius : float, unit=metre
        distance away from Earths' center
    inclination : unit=degrees
        tilt of the orbital plane in relation to the equator
    period : float, unit=seconds
        time it takes to complete one revolution
    time_para : numpy.ndarray, size=(p,b)
        polynomial fitting parameters for the different bands (b)
    combos : numpy.ndarray, size=(h,2)
        combinations of band and detector pairs,
        corresponding to the third dimension of 'grid.zenith'
    """
    crs = osr.SpatialReference()
    crs.ImportFromEPSG(grid.epsg)

    if is_crs_an_srs(crs):
        x_grd, y_grd = pix_centers(grid.geotransform, rows=grid.rows, cols=grid.columns, make_grid=True)
        ll_grd = map2ll(np.stack((x_grd.ravel(), y_grd.ravel()), axis=1), crs)
        lat_arr_grd,lon_arr_grd = ll_grd[:,0].reshape((grid.columns, grid.rows)), \
                          ll_grd[:,1].reshape((grid.columns, grid.rows))
        del ll_grd
    else:
        lon_arr_grd, lat_arr_grd = pix_centers(grid.geotransform, rows=grid.rows, cols=grid.columns, make_grid=True)

    # remove NaN's, and create vectors
    ok = np.invert(np.isnan(grid.azimuth))
    lat_arr_grd = np.tile(np.atleast_3d(lat_arr_grd), (1, 1, grid.depth))
    lon_arr_grd = np.tile(np.atleast_3d(lon_arr_grd), (1, 1, grid.depth))
    lat_arr, lon_arr = lat_arr_grd[ok], lon_arr_grd[ok]
    az_arr, zn_arr = grid.azimuth[ok], grid.zenith[ok]

    sat, g_x = line_of_sight(lat_arr, lon_arr, zn_arr, az_arr)
    del lat_arr, lon_arr, lat_arr_grd, lon_arr_grd

    l_time, lat, lon, radius, inclination, period = orbital_fitting(sat, g_x, inclination,
            revolutions_per_day=revolutions_per_day, mean_altitude=mean_altitude, radius=radius)

    # vectorize band and detector indicators
    Bnd, Det = (np.tile(grid.band,(grid.columns, grid.rows, 1))[ok],
                np.tile(grid.detector,(grid.columns, grid.rows, 1))[ok])
    x_arr, y_arr = np.tile(np.atleast_3d(x_grd), (1, 1, grid.depth))[ok], \
           np.tile(np.atleast_3d(y_grd), (1, 1, grid.depth))[ok]
    time_para, combos = time_fitting(l_time, az_arr, zn_arr, Bnd, Det,
                                     x_arr, y_arr, grid.geotransform)

    return lat, lon, radius, inclination, period, time_para, combos

def remap_observation_angles(lat, lon, radius, inclination, period,
                             time_para, combos, x_grd, y_grd, det_stack,
                             bnd_list, geotransform, crs):
    if type(bnd_list) in (pd.core.frame.DataFrame,):
        bnd_list = np.asarray(bnd_list['bandid'])
        bnd_list -= 1 # numbering of python starts at 0
    lat,lon = lat_lon_angle_check(lat,lon)
    are_two_arrays_equal(x_grd,y_grd)
    geotransform = correct_geotransform(geotransform)

    m,n = x_grd.shape

    omega_0,lon_0 = _omega_lon_calculation(np.deg2rad(lat), np.deg2rad(lon),
                                           inclination)

    # reconstruct observation angles and sensing time
    b = bnd_list.size
    if type(det_stack) in (np.ma.core.MaskedArray,):
        T, zn_arr, az_arr = np.ma.zeros((m,n,b)), np.ma.zeros((m,n,b)), \
                    np.ma.zeros((m,n,b))
    else:
        T, zn_arr, az_arr = np.zeros((m,n,b)), np.zeros((m,n,b)), \
                    np.zeros((m,n,b))
    for idx, bnd in enumerate(bnd_list):
        doi = combos[:,0]==bnd
        if type(det_stack) in (np.ma.core.MaskedArray,):
            dt_bnd, zn_bnd, az_bnd = -9999.*np.ma.ones((m,n)), \
                                     -9999.*np.ma.ones((m,n)),\
                                     -9999.*np.ma.ones((m,n))
        else:
            dt_bnd, zn_bnd, az_bnd = np.zeros((m,n)), np.zeros((m,n)), \
                                     np.zeros((m,n))

        for sca in combos[doi,1]:
            ok = (det_stack[...,idx] == sca)
            if not np.any(ok): continue
            dx, dy = x_grd[ok]-geotransform[0], geotransform[3]-y_grd[ok]

            # time stamps
            coef_id = np.where(np.logical_and(combos[:,0]==bnd,
                                              combos[:,1]==sca))[0][0]
            coeffs = time_para[coef_id,:]
            dt = coeffs[0] + coeffs[1]*dx + coeffs[2]*dy + coeffs[3]*dx*dy
            dt_bnd[ok] = dt
            del dx, dy, coeffs, coef_id
            # acquisition angles
            p_x = orbital_calculation(dt, radius, inclination, period,
                                     omega_0, lon_0) # satellite vector

            ll_pix = map2ll(np.stack((x_grd[ok], y_grd[ok]), axis=1), crs)
            g_x = np.transpose(ground_vec(ll_pix[:, 0], ll_pix[:, 1]))  # ground vector
            del ll_pix, dt

            zn, az = acquisition_angles(p_x,g_x)
            zn_bnd[ok], az_bnd[ok] = zn, az
            del p_x,g_x
        # put estimates in stack

        if type(det_stack) in (np.ma.core.MaskedArray,):
            dt_bnd = np.ma.array(dt_bnd, mask=dt_bnd == -9999.)
            zn_bnd = np.ma.array(zn_bnd, mask=zn_bnd == -9999.)
            az_bnd = np.ma.array(az_bnd, mask=az_bnd == -9999.)
        T[...,idx], zn_arr[...,idx], az_arr[...,idx] = dt_bnd, zn_bnd, az_bnd
        del zn_bnd, az_bnd
    return zn_arr, az_arr, T

def get_absolute_timing(lat,lon,sat_dict):
    assert isinstance(sat_dict, dict), 'please provide a dictionary'
    assert ('gps_xyz' in sat_dict.keys()), 'please include GPS metadata'
    assert ('gps_tim' in sat_dict.keys()), 'please include GPS metadata'
    lat, lon = lat_lon_angle_check(lat, lon)

    cen_ll = np.array([lat, lon])
    sat_llh = ecef2llh(sat_dict['gps_xyz'])
    ll_diff = sat_llh[:,:2] - cen_ll[None, :]
    ll_dist = np.linalg.norm(ll_diff, axis=1)

    min_idx = np.argmin(ll_dist)
    # sub-grid interpolation
    dd = np.divide(ll_dist[min_idx+1]-ll_dist[min_idx-1],
                   2*((2*ll_dist[min_idx])
                      - ll_dist[min_idx-1] - ll_dist[min_idx+1]))

    T = sat_dict['gps_tim']

    dT = np.abs(T[min_idx]-T[int(min_idx+np.sign(dd))])
    t_bias = T[min_idx] + dd*dT
    return t_bias

def acquisition_angles(p_x: np.ndarray,
                       g_x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """ given satellite and ground coordinates, estimate observation angles"""
    are_two_arrays_equal(p_x, g_x)

    major_axis,minor_axis = earth_axes()
    v_x = p_x - g_x  # observation vector
    del p_x
    v_dist = np.linalg.norm(v_x, axis=1)  # make unit length
    v_x = np.einsum('i...,i->i...', v_x, np.divide(1, v_dist))
    del v_dist

    e_z = np.einsum('...i,i->...i', g_x,
                    1 / np.array([major_axis, major_axis, minor_axis]))
    e_e = np.zeros_like(e_z)
    e_e[..., 0], e_e[..., 1] = -e_z[:, 1].copy(), e_z[:, 0].copy()
    e_plan = np.linalg.norm(e_z[:, :2], axis=1)
    e_e = np.einsum('i...,i->i...', e_e, np.divide(1, e_plan))
    del e_plan
    e_n = np.array([np.multiply(e_z[:, 1], e_e[:, 2]) -
                    np.multiply(e_z[:, 2], e_e[:, 1]),
                    np.multiply(e_z[:, 2], e_e[:, 0]) -
                    np.multiply(e_z[:, 0], e_e[:, 2]),
                    np.multiply(e_z[:, 0], e_e[:, 1]) -
                    np.multiply(e_z[:, 1], e_e[:, 0])]).T

    los = np.zeros_like(e_z)
    e_str = '...i,...i->...'
    los[..., 0] = np.einsum(e_str, v_x, e_e)
    del e_e
    los[..., 1] = np.einsum(e_str, v_x, e_n)
    del e_n
    los[..., 2] = np.einsum(e_str, v_x, e_z)
    del e_z

    az = np.rad2deg(np.arctan2(los[..., 0], los[..., 1]))
    zn = np.rad2deg(np.arccos(los[...,2]))
    return zn, az

def line_of_sight(lat_arr: np.ndarray,
                  lon_arr: np.ndarray,
                  zn_arr: np.ndarray,
                  az_arr: np.ndarray,
                  eccentricity: Optional[float] = None,
                  major_axis: Optional[float] = None) -> tuple[np.ndarray, np.ndarray]:
    """

    Parameters
    ----------
    lat_arr, lon_arr : numpy.ndarray, size=(m,n), unit=degrees
        location of observation angles
    zn_arr, az_arr : numpy.ndarray, size=(m,n), unit=degrees
        polar angles of line of sight to the satellite, also known as,
        delcination and right ascension.
    eccentricity: float
        eccentricity squared, if None default is WGS84 (6.69E-10)
    major_axis: float, unit=meters
        length of the major axis, if None default is WGS84 (6.37E6)

    Returns
    -------
    sat : numpy.ndarray, size=(m*n,3)
        observation vector from ground location towards the satellite
    g_x : numpy.ndarray, size=(m*n,3)
        ground location in Cartesian coordinates
    """
    lat_arr, lon_arr = lat_lon_angle_check(lat_arr, lon_arr)
    assert len(set({lat_arr.shape[0],
                    lon_arr.shape[0],
                    zn_arr.shape[0],
                    az_arr.shape[0]}))==1,\
        ('please provide arrays of the same size')

    if (major_axis is None) or (eccentricity is None):
        major_axis = wgs84_param()[0]
        eccentricity = earth_eccentricity()

    lat_arr, lon_arr = np.deg2rad(lat_arr.flatten()), np.deg2rad(lon_arr.flatten())
    zn_arr, az_arr = np.deg2rad(zn_arr.flatten()), np.deg2rad(az_arr.flatten())

    # local tangent coordinate system
    e_e = np.stack((-np.sin(lon_arr),
                    +np.cos(lon_arr),
                    np.zeros_like(lon_arr)), axis=1)
    e_n = np.stack((-np.sin(lat_arr)*np.cos(lon_arr),
                    -np.sin(lat_arr)*np.sin(lon_arr),
                    +np.cos(lat_arr)), axis=1)
    e_z = np.stack((np.cos(lat_arr)*np.cos(lon_arr),
                    np.cos(lat_arr)*np.sin(lon_arr),
                    np.sin(lat_arr)), axis=1)
    los = np.stack((np.sin(zn_arr)*np.sin(az_arr),
                    np.sin(zn_arr)*np.cos(az_arr),
                    np.cos(zn_arr)), axis=1)

    sat = np.array([los[:,0]*e_e[:,0] + los[:,1]*e_n[:,0] + los[:,2]*e_z[:,0],
                    los[:,0]*e_e[:,1] + los[:,1]*e_n[:,1] + los[:,2]*e_z[:,1],
                    los[:,0]*e_e[:,2] + los[:,1]*e_n[:,2] + los[:,2]*e_z[:,2]]
                   ).T
    radi = np.divide(major_axis,
                     np.sqrt( 1.0 - eccentricity*np.sin(lat_arr)*np.sin(lat_arr)))
    g_x = np.array([radi * np.cos(lat_arr) * np.cos(lon_arr),
                    radi * np.cos(lat_arr) * np.sin(lon_arr),
                    radi * (1-eccentricity) * np.sin(lat_arr)]).T
    return sat, g_x

def ground_vec(lat_arr: np.ndarray,
               lon_arr: np.ndarray,
               eccentricity: Optional[float] = None,
               major_axis: Optional[float] = None) -> np.ndarray:
    """ get ground coordinates in Cartesian system

    Parameters
    ----------
    lat_arr, lon_arr : numpy.ndarray, size=(m,n), unit=degrees
        location of observation angles
    eccentricity: float
        eccentricity squared, if None default is WGS84 (6.69E-10)
    major_axis: float, unit=meters
        length of the major axis, if None default is WGS84 (6.37E6)

    Returns
    -------
    g_x : numpy.ndarray, size=(m*n,3)
        ground location in Cartesian coordinates

    """
    lat_arr, lon_arr = lat_lon_angle_check(lat_arr, lon_arr)
    are_two_arrays_equal(lat_arr,lon_arr)
    if (major_axis is None) or (eccentricity is None):
        major_axis, flattening = wgs84_param()
        eccentricity = (2*flattening) - (flattening**2)

    lat_arr, lon_arr = np.deg2rad(lat_arr.flatten()), np.deg2rad(lon_arr.flatten())
    radi = np.divide(major_axis,
                       np.sqrt(1.0 - eccentricity*np.sin(lat_arr)*np.sin(lat_arr)))
    g_x = np.array([radi*np.cos(lat_arr)*np.cos(lon_arr),
                   radi*np.cos(lat_arr)*np.sin(lon_arr),
                   radi* (1 - eccentricity) * np.sin(lat_arr)])
    return g_x

def _make_timing_system(ok,x,y,l_time):
    are_three_arrays_equal(x, y, l_time)
    x, y, t = x[ok], y[ok], l_time[ok]

    lt = np.array([t, x*t, y*t, x*y*t])
    l = np.sum(lt, axis=1)

    dx, dy, n = np.sum(x), np.sum(y), np.sum(ok)
    dx2, dxy, dy2 = np.sum(x**2), np.sum(x*y), np.sum(y**2)
    dx2y, dxy2, dx2y2 = np.sum(x*x*y), np.sum(x*y*y), np.sum(x*x*y*y)

    A = np.array([[  n,   dx,   dy,   dxy],
                  [ dx,  dx2,  dxy,  dx2y],
                  [ dy,  dxy,  dy2,  dxy2],
                  [dxy, dx2y, dxy2, dx2y2]])
    return A, l

def time_fitting(l_time, az_arr, zn_arr, bnd, det, x, y, geotransform):
    are_two_arrays_equal(zn_arr,az_arr)
    are_two_arrays_equal(bnd, det)
    are_two_arrays_equal(x,y)
    geotransform = correct_geotransform(geotransform)

    # translate coordinates
    dx, dy = x - geotransform[0], geotransform[3] - y

    # look per band and detector
    combos, idx_inv = np.unique(np.stack((bnd,det), axis=1), axis=0,
                               return_inverse=True)
    x_hat = np.zeros((combos.shape[0],4))
    for idx, pair in enumerate(combos):
        ok = np.logical_and(bnd==pair[0], det==pair[1])
        a_sca, l_sca = _make_timing_system(ok, dx, dy, l_time)
        try:
            np.linalg.inv(a_sca)
        except:
            if 'a_bnd' not in locals():
                ok = bnd == pair[0]
                a_bnd, l_bnd = _make_timing_system(ok, dx, dy, l_time)
            a_sca,l_sca = a_bnd,l_bnd
        x_tilde = np.transpose(np.linalg.inv(a_sca)@l_sca[:,np.newaxis])
        x_hat[idx,...] = x_tilde
        if 'a_bnd' in locals():
            del a_bnd,l_bnd
    return x_hat, combos

def orbital_fitting(sat, g_x, inclination, lat=None, lon=None, radius=None,
                    period=None, revolutions_per_day=None, mean_altitude=None,
                    gps=None, convtol = 0.001, orbtol=1.0, maxiter=20,
                    printing=False):
    """

    Parameters
    ----------
    sat : numpy.ndarray, size=(m*n,3)
        observation vector from ground location towards the satellite
    g_x : numpy.ndarray, size=(m*n,3)
        ground location in Cartesian coordinates
    lat : float, unit=degrees
        location of satellite within orbital track
    lon : float, unit=degrees
        location of satellite within orbital track
    radius : float, unit=meter
        radius towards the orbiting satellite
    inclination : float, unit=degrees, range=-180...+180
        inclination of the orbital plane with the equatorial plane
    period : float, unit=seconds
        time it takes to revolve one time around the Earth

    Returns
    -------

    """
    are_two_arrays_equal(sat, g_x)
    inclination = np.deg2rad(inclination)
    if (period is None) and (revolutions_per_day is not None):
        period = (24*60*60) / revolutions_per_day
    if radius is None:
        major_axis = wgs84_param()[0]
        if mean_altitude is not None:
            radius = major_axis + mean_altitude
        else:
            radius = np.cbrt(np.divide(period, 2*np.pi)**2 * standard_gravity())

    if (lat is None) or (lon is None):
        v_dist = np.sqrt( radius**2 +
                          np.einsum('...i,...i', sat, g_x)**2 -
                          np.linalg.norm(g_x, axis=1)**2 )
        p_x = g_x + np.einsum('i,ij->ij', v_dist, sat)
        if gps is not None: # use GPS trajectory
            poi = np.argmin(np.linalg.norm(gps -
                                           np.mean(p_x, axis=0), axis=1))
            s_x = gps[poi,:]
            lat_bar = np.arctan2( s_x[2], np.linalg.norm(s_x[0:2]))
            lon_bar = np.arctan2( s_x[1], s_x[0] )
            del poi, s_x
        else: # use center of the scene as approximation
            lon_arr = np.arctan2( p_x[...,1], p_x[...,0] )
            lat_arr = np.arctan2( p_x[...,2], np.linalg.norm(p_x[...,0:2], axis=1))
            lat_bar, lon_bar = np.mean(lat_arr), np.mean(lon_arr)
            del lat_arr, lon_arr
        del p_x, v_dist

    numobs = sat.shape[0]
    l_time = np.zeros((numobs,))

    rmstime, orbrss = 15.0, 1000.0
    counter = 0
    while (rmstime > convtol) or (orbrss > orbtol) or (counter < maxiter):
        omega_0, lon_0 = _omega_lon_calculation(lat_bar, lon_bar, inclination)

        A_0, L_0 = np.zeros((4,4,numobs)), np.zeros((4,numobs))

        v_x = observation_calculation(l_time, sat, g_x, radius, inclination,
                                     period, omega_0, lon_0)

        # Calculate the partial derivatives w.r.t. the orbit parameters
        p_0 = partial_obs(l_time, sat, g_x, lat_bar, lon_bar, radius,
                          inclination, period)
        A_0 += np.einsum('ij...,ik...->jk...', p_0, p_0)
        L_0 += np.einsum('...i,ij...->j...', v_x, p_0)

        p_1 = partial_tim(l_time, sat, g_x,
                          lat_bar, lon_bar, radius,
                          inclination, period)
        p_1 = np.squeeze(p_1)
        M_1 = np.einsum('ij...,i...->j...', p_0, p_1)
        A_1 = np.reciprocal(np.einsum('i...,i...->...', p_1, p_1))
        L_1 = np.multiply(np.einsum('i...,...i->...',p_1, v_x), A_1)
        del p_0, p_1, v_x

        A_0 += np.einsum('i...,j...->ij...', A_1*M_1, M_1)
        L_0 += np.multiply(M_1,L_1)

        A_0, L_0 = np.sum(A_0, axis=2), np.sum(L_0, axis=1)[:,np.newaxis]
        N_1 = np.multiply(M_1,A_1)
        del M_1, A_1

        if counter!=0:
            X_0 = np.einsum('ij,i...->j...', np.linalg.inv(A_0), L_0)
        else:
            X_0 = np.zeros((4,))

        # back substitute for time corrections
        dtime = L_1 - np.einsum('i...,i...->...', N_1, X_0)
        l_time -= dtime
        rmstime = np.sum(dtime**2)
        del dtime, N_1, L_1

        # update orbit parameters
        lat_bar -= X_0[0]
        lon_bar -= X_0[1]
        radius -= X_0[2]
        inclination -= X_0[3]

        # evaluate convergence
        rmstime = np.sqrt(rmstime / numobs)
        # orbit convergence, transform to metric units
        X_0[0] *= 6378137.0
        X_0[1] *= 6378137.0
        X_0[3] *= radius
        orbrss = np.linalg.norm(X_0)
        counter += 1
        if printing:
            print('RMS Orbit Fit (meters): ', orbrss)
            print('RMS Time Fit (seconds): ', rmstime)
    lat_bar, lon_bar = np.rad2deg(lat_bar[0]), np.rad2deg(lon_bar[0])
    return l_time, lat_bar, lon_bar, radius, inclination, period

def _omega_lon_calculation(lat, lon, inclination):
    """

    Parameters
    ----------
    lat, lon : {float, numpy.array}, unit=radians
        location on the Earth
    inclination : {float, numpy.array}, unit=radians
        angle of the orbital plane i.r.t. the equator
    """
    omega_0 = np.arcsin(np.divide(np.sin(lat), np.sin(inclination)))
    lon_0 = lon - np.arcsin(np.divide(np.tan(lat), -np.tan(inclination)))
    return omega_0, lon_0

def _gc_calculation(ltime,period,inclination, omega_0, lon_0):
    cta = omega_0 - np.divide(2 * np.pi * ltime, period)
    gclat = np.arcsin(np.sin(cta) * np.sin(inclination))
    gclon = lon_0 + np.arcsin(np.tan(gclat) / -np.tan(inclination)) - \
        2*np.pi*ltime / (24*60*60)
    return cta, gclat, gclon

def orbital_calculation(ltime,radius,inclination,period,
                        omega_0, lon_0):
    ltime, omega_0, lon_0 = np.squeeze(ltime), np.squeeze(omega_0), np.squeeze(lon_0)
    cta, gclat, gclon = _gc_calculation(ltime, period, inclination,
                                        omega_0, lon_0)
    p_x = np.stack((np.multiply(np.cos(gclat), np.cos(gclon)),
                   np.multiply(np.cos(gclat), np.sin(gclon)),
                   np.sin(gclat)), axis=1)
    p_x *= radius
    return p_x

def observation_calculation(ltime, sat, g_x, radius, inclination,
                            period, omega_0, lon_0):
    """

    Parameters
    ----------
    ltime : numpy.array, size=(m,1)
    sat : numpy.ndarray, size=(m,3)
        observation vector from ground location towards the satellite
    g_x : numpy.ndarray, size=(m,3)
        ground location in Cartesian coordinates
    radius : float, unit=meter
        radius towards the orbiting satellite
    inclination : float, unit=degrees, range=-180...+180
        inclination of the orbital plane with the equator
    period : float, unit=seconds
        time it takes to revolve one time around the Earth
    omega_0 : float, unit=degrees
        angle towards ascending node, one of the orbit Euler angles
    lon_0 : float, unit=degrees
        ephemeris longitude

    Returns
    -------
    v_x : numpy.array, size=(m,3)
    """
    are_two_arrays_equal(sat, g_x)

    cta, gclat, gclon = _gc_calculation(ltime, period, inclination, omega_0, lon_0)
    v_x = np.atleast_2d(np.zeros_like(sat))
    g_x = np.atleast_2d(g_x)
    v_x[...,0] = np.squeeze(radius * np.multiply(np.cos(gclat), np.cos(gclon)))
    v_x[...,1] = np.squeeze(radius * np.multiply(np.cos(gclat), np.sin(gclon)))
    v_x[...,2] = np.squeeze(radius * np.sin(gclat))
    v_x -= g_x

    v_dist = np.linalg.norm(v_x, axis=1) # make unit length
    v_x = np.einsum('i...,i->i...', v_x, np.divide(1, v_dist))
    v_x -= sat
    return v_x

def _pert_param(idx, pert, *args):
    """ perterp one of the arguments

    Parameters
    ----------
    idx : integer, range=0...m-1
        index which of the arguments should be perturbed
    pert : {float, numpy.array}
        amount of pertubation
    args : tuple, size=m
        tuple with different arguments, with datatypes like arrays or floats

    Returns
    -------
    args : tuple, size=m
        tuple with different arguments, with datatypes like arrays or floats
    """
    assert idx<len(args), 'please provide correct index'
    args = list(args)
    args[idx] += pert
    args = tuple(args)
    return args

def partial_obs(ltime, sat, g_x, lat, lon, radius, inclination, period):
    """ numerical differentiation, via pertubation of the observation vector
    """
    are_two_arrays_equal(sat, g_x)
    p_0 = np.zeros((3, 4, ltime.size))
    omega_0, lon_0 = _omega_lon_calculation(lat, lon, inclination)
    dx = observation_calculation(ltime, sat, g_x, radius, inclination,
                                 period, omega_0, lon_0)
    pert_var = ['lat', 'lon', 'radius', 'inclination']
    pert = np.array([1E-5, 1E-5, 1E+1, 1E-4])
    for idx,pert in enumerate(pert):
        (lat, lon, radius, inclination) = _pert_param(idx, +pert, lat, lon,
                                                      radius, inclination)

        omega_0, lon_0 = _omega_lon_calculation(lat, lon, inclination)
        dp = observation_calculation(ltime, sat, g_x, radius,
                                     inclination, period, omega_0, lon_0)
        p_0[0,idx,:] = np.divide(dp[:,0] - dx[:,0],pert)
        p_0[1,idx,:] = np.divide(dp[:,1] - dx[:,1],pert)
        p_0[2,idx,:] = np.divide(dp[:,2] - dx[:,2],pert)
        (lat, lon, radius, inclination) = _pert_param(idx, -pert, lat, lon,
                                                      radius, inclination)
    return p_0

def partial_tim(ltime, sat, g_x, lat, lon, radius, inclination, period,
                pertubation=.1):
    are_two_arrays_equal(sat, g_x)
    p_1 = np.zeros((3, 1, ltime.size))
    omega_0, lon_0 = _omega_lon_calculation(lat, lon, inclination)
    dx = observation_calculation(ltime, sat, g_x, radius, inclination,
                            period, omega_0, lon_0)

    # pertubation in the time domain
    ltime += pertubation
    dp = observation_calculation(ltime, sat, g_x, radius, inclination,
                            period, omega_0, lon_0)
    ltime -= pertubation

    p_1[0,0,...] = np.divide(dp[...,0] - dx[...,0], pertubation)
    p_1[1,0,...] = np.divide(dp[...,1] - dx[...,1], pertubation)
    p_1[2,0,...] = np.divide(dp[...,2] - dx[...,2], pertubation)
    return p_1
