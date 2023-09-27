#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Space Accountants"
__license__ = "MIT License - You must cite this source"
__version__ = "202309"
__maintainer__ = "B. Altena"
__email__ = "info at space hyphen accountants dot eu"

import numpy as np

from osgeo import osr

from naming_check import check_mgrs_code
from mapping_tools import ll2map, get_utm_zone

def _get_mgrs_abc():
    mgrsABC = [chr(i) for i in list(range(65,73)) +
               list(range(74,79)) + list(range(80,91))]
    return mgrsABC

def get_mgrs_tile(ϕ,λ):
    """ return a military grid reference system zone designation string.
    This zoning is used in the tiling of Sentinel-2.

    Parameters
    ----------
    ϕ : float, unit=degrees, range=-90...+90
        latitude
    λ : float, unit=degrees
        longitude

    Returns
    -------
    mgrs_code : string

    Notes
    -----
    The MGRS tile structure is a follows "AABCC"
        * "AA" utm zone number, starting from the East, with steps of 8 degrees
        * "B" latitude zone, starting from the South, with steps of 6 degrees
        * "CC" 100 km square identifier

    The following acronyms are used:
        - CRS : coordinate reference system
        - MGRS : US military grid reference system
        - WKT : well known text
        - UTM: universal transverse mercator projection
    """

    # numbering goes with the alphabet, excluding "O" and "I"
    mgrsABC = _get_mgrs_abc()
    tile_size = 100. # [km]
    tile_size *= 1E3

    utm_zone = get_utm_zone(ϕ, λ)
    utm_no = int(utm_zone[:-1])

    # transform to UTM coordinates
    proj = osr.SpatialReference()
    NH = True if ϕ>0 else False
    proj.SetUTM(utm_no, NH)
    falseEasting = proj.GetProjParm( osr.SRS_PP_FALSE_EASTING)

    xyz = np.squeeze(ll2map(np.array([[ϕ,λ]]), proj))

    # λ letter
    shift_letter = np.floor( np.divide(xyz[0] - falseEasting,
                                       tile_size)).astype(int)
    center_letter = int(np.mod(utm_no - 1, 3) * 8) + 4
    λ_letter = mgrsABC[center_letter + shift_letter]

    # ϕ letter
    tile_shift = np.fix(xyz[1])
    ϕ_num = np.mod(tile_shift, 20)
    if np.mod(utm_no,2)==0: # even get an additional five letter shift
        if np.all((xyz[1] < 0, np.abs(tile_shift) > 5)):
            ϕ_num -= 4 # counts up to "R"
        else:
            ϕ_num += 5
    else:
        if xyz[1] < 0: # a leap hole is present in the southern hemisphere
            ϕ_num -= 4 # counts up to "R"
    ϕ_idx = np.mod(ϕ_num, 20)
    ϕ_idx -= 1
    ϕ_letter = mgrsABC[ϕ_idx]

    mgrs_code = utm_zone + λ_letter + ϕ_letter
    return mgrs_code

def get_geom_mgrs_tile(mgrs_code):
    utm_zone = mgrs_code[:2]
    λ_letter = mgrs_code[-2]
    ϕ_letter = mgrs_code[-1]

    return
