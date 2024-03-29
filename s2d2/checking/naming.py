#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import re


# sentinel relative orbit numbers go from 1 to 143

def check_mgrs_code(code_str):
    """ validate a MGRS tile code and make it uppercase

    Parameters
    ----------
    tile_code : str
        MGRS tile code

    Returns
    -------
    code_str : str
        validated and normalized tile code

    Notes
    -----
    The following acronyms are used:

        - MGRS : US Military Grid Reference System
    """

    if not isinstance(code_str, str):
        raise TypeError("please provide a string")

    tile_code = code_str.upper()
    pattern = "[0-9][0-9][A-Z][A-Z][A-Z]"

    if not bool(re.match(pattern, tile_code)):
        raise ValueError("please provide a correct MGRS tile code")

    return tile_code

def check_sen2_safe_code(code_str):
    """ check if the .SAFE code string is correct

    Parameters
    ----------
    code_str: str, len=65
        folder name where L1C Sentinel-2 data is situated
    """
    if not isinstance(code_str, str):
        raise TypeError("please provide a string")

    code_str = code_str.upper()
    pattern = "S2A_MSIL1C_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]" + \
              "T[0-9][0-9][0-9][0-9][0-9][0-9]_" + \
              "N[0-9][0-9][0-9][0-9]_" + \
              "R[0-9][0-9][0-9]_" + \
              "T[0-9][0-9][A-Z][A-Z][A-Z]_" + \
              "[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]" + \
              "T[0-9][0-9][0-9][0-9][0-9][0-9].SAFE"

    if not bool(re.match(pattern, code_str)):
        raise ValueError("please provide a correct Sentinel-2 tile code")
    return code_str
