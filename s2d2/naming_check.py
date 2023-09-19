#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Space Accountants"
__license__ = "MIT License - You must cite this source"
__version__ = "202309"
__maintainer__ = "B. Altena"
__email__ = "info at space hyphen accountants dot eu"

import re

def check_mgrs_code(tile_code):
    """
    Validate a MGRS tile code and make it uppercase

    Parameters
    ----------
    tile_code : str
        MGRS tile code

    Returns
    -------
    tile_code : str
        validated and normalized tile code

    Notes
    -----
    The following acronyms are used:

        - MGRS : US Military Grid Reference System
    """

    if not isinstance(tile_code, str):
        raise TypeError("please provide a string")

    tile_code = tile_code.upper()

    if not bool(re.match("[0-9][0-9][A-Z][A-Z][A-Z]", tile_code)):
        raise ValueError("please provide a correct MGRS tile code")

    return tile_code
