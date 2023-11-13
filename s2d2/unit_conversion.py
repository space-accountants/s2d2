#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Space Accountants"
__license__ = "MIT License - You must cite this source"
__version__ = "202311"
__maintainer__ = "B. Altena"
__email__ = "info at space hyphen accountants dot eu"

import numpy as np

from datetime import datetime, timedelta

# time conversions
def datetime2doy(dt):
    """ Convert array of datetime64 to a day of year.

    Parameters
    ----------
    dt : np.datetime64
        times of interest

    Returns
    -------
    year : integer, {x ∈ ℕ}
        calender year
    doy : integer, {x ∈ ℕ}
        day of year
    """
    year = dt.astype('datetime64[Y]').astype(int)+1970
    doy = (dt.astype('datetime64[D]') -
          dt.astype('datetime64[Y]') + 1).astype('float64')
    return year, doy

def doy2dmy(doy,year):
    """ given day of years, calculate the day and month

    Parameters
    ----------
    year : integer, {x ∈ ℕ}
        calender year
    doy : integer, {x ∈ ℕ}
        day of year

    Returns
    -------
    day, month, year
    """

    if year.size < doy.size:
        year = year*np.ones_like(doy)

    # sometimes the doy is longer than 365, than put this to the year counter
    xtra_years = np.floor(doy / 365).astype(int)
    year += xtra_years

    doy = doy % 365

    dates = (year-1970).astype('datetime64[Y]') + doy.astype('timedelta64[D]')
    month = dates.astype('datetime64[M]').astype(int) % 12 + 1
    day = dates - dates.astype('datetime64[M]') + 1
    day = day.astype(int)
    if day.size==1:
        return day[0], month[0], year[0]
    else:
        return day, month, year

def datetime2calender(dt):
    """ Convert array of datetime64 to a calendar year, month, day.

    Parameters
    ----------
    dt : np.datetime64
        times of interest

    Returns
    -------
    cal : numpy array
        calendar array with last axis representing year, month, day
    """

    # decompose calendar floors
    Y, M, D = [dt.astype(f"M8[{x}]") for x in "YMD"]
    year = (Y + 1970).astype('timedelta64[Y]').astype(int)
    month = ((M - Y) + 1).astype('timedelta64[M]').astype(int)
    day = ((D - M) + 1).astype('timedelta64[D]').astype(int)
    return year, month, day

def datenum2datetime(t):
    """ some write a date as a number "YYYYMMDD" conver this to numpy.datetime

    Parameters
    ----------
    t : {int,numpy.array}
        time, in format YYYYMMD

    Returns
    -------
    dt : {numpy.datetime64, numpy.array}, type=numpy.datetime64[D]
        times in numpy time format

    See Also
    --------
    datetime2datenum
    """
    if type(t) in (np.ndarray, ):
        dt = np.array([np.datetime64(str(e)[:4]+'-'+str(e)[4:6]+'-'+str(e)[-2:])
              for e in t.astype(str)])
    else:
        dt = np.datetime64(str(t)[:4] + '-' + str(t)[4:6] + '-' + str(t)[-2:])
    return dt

def datetime2datenum(dt):
    """ convert numpy.datetime to a date as a number "YYYYMMDD"

    Parameters
    ----------
    dt : {numpy.datetime64, numpy.array}, type=numpy.datetime64[D]
        times in numpy time format

    Returns
    -------
    t : {int,numpy.array}
        time, in format YYYYMMD

    See Also
    --------
    datenum2datetime
    """
    if type(dt) in (np.ndarray, ):
        t = np.array([int(e.astype(str)[:4] +
                          e.astype(str)[5:7]+
                          e.astype(str)[-2:])
              for e in dt.astype(str)], dtype=np.int64)
    else:
        dt = int(dt.astype(str)[:4]+dt.astype(str)[5:7]+dt.astype(str)[-2:])
    return t

# geometric and angular scales and formats
def deg2dms(ang):
    """ convert decimal degrees to degree minutes seconds format

    Parameters
    ----------
    ang : {float,np.array}, unit=decimal degrees
        angle(s) of interest

    Returns
    -------
    deg : {integer,np.array}
        degrees
    min : {integer,np.array}, range=0...60
        angular minutes
    sec : {float,np.array}, range=0...60
        angular seconds
    """
    min,sec = np.divmod(ang*3600, 60)
    deg,min = np.divmod(min, 60)
    deg,min = deg.astype(int), min.astype(int)
    return deg,min,sec

def dms2deg(deg,min,sec):
    """ convert degree minutes seconds format to decimal degrees

    Parameters
    ----------
    deg : {integer,np.array}
        degrees
    min : {integer,np.array}, range=0...60
        angular minutes
    sec : {float,np.array}, range=0...60
        angular seconds

    Returns
    -------
    ang : {float,np.array}, unit=decimal degrees
        angle(s) of interest
    """
    ang = deg.astype(float) + (min.astype(float)/60) + (sec.astype(float)/3600)
    return ang

def deg2compass(θ):
    """ adjust angle to be in bounds of a positive argument angle,like a compass

    Parameters
    ----------
    θ : unit=degrees

    Returns
    -------
    θ : unit=degrees, range=0...+360

    See Also
    --------
    deg2arg

    Notes
    -----
    The angle is declared in the following coordinate frame:

        .. code-block:: text

                 ^ North & y
                 |
            - <--|--> +
                 |
                 +----> East & x
    """
    return θ % 360

def deg2arg(θ):
    """ adjust angle to be in bounds of an argument angle

    Parameters
    ----------
    θ : unit=degrees

    Returns
    -------
    θ : unit=degrees, range=-180...+180

    See Also
    --------
    deg2compass

    Notes
    -----
    The angle is declared in the following coordinate frame:

        .. code-block:: text

                 ^ North & y
                 |
            - <--|--> +
                 |
                 +----> East & x
    """
    return ((θ + 180) % 360) -180
