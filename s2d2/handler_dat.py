#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Space Accountants"
__license__ = "MIT License - You must cite this source"
__version__ = "202309"
__maintainer__ = "B. Altena"
__email__ = "info at space hyphen accountants dot eu"

import os
import re

def get_list_files(dat_path, ending):
    """ generate a list with files, with a given ending/extension

    Parameters
    ----------
    dat_path : string
        path of interest.
    ending : string
        ending of a file.

    Returns
    -------
    file_list : list
        collection of file names.

    """
    if not os.path.isdir(dat_path):
        dat_path = os.path.dirname(dat_path)
    file_list = [x for x in os.listdir(dat_path) if x.endswith(ending)]
    return file_list
