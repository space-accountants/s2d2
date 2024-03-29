#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import os

import xml.etree.ElementTree as ElementTree
import numpy as np


def get_array_from_xml(tree_struc: ElementTree.Element) -> np.ndarray:
    """ transform data in xml-string to numpy array

    Parameters
    ----------
    tree_struc : string
        xml structure

    Returns
    -------
    np.ndarray
        data array
    """
    for i in range(0, len(tree_struc)):
        t_row = [float(s) for s in tree_struc[i].text.split(' ')]
        if i == 0:
            t_n = t_row
        elif i == 1:
            t_row = np.stack((t_row, t_row), 0)
            t_n = np.stack((t_n, t_row[1, :]), 0)
        else:
            t_n = np.concatenate((t_n, [t_row]), 0)
    return t_n


def get_root_of_table(path: str,
                      fname: str = None) -> ElementTree:
    if fname is None:
        full_name = path
    else:
        full_name = os.path.join(path, fname)
    if '*' not in full_name:
        assert os.path.isfile(full_name), \
            ('please provide correct path and file name')
    dom = ElementTree.parse(glob.glob(full_name)[0])
    root = dom.getroot()
    return root


def get_branch(tree_struct: ElementTree,
               syntax_ending: str) -> ElementTree.Element:
        branch_struct = None
        for child in tree_struct:
            if child.tag.endswith(syntax_ending):
                return child
        assert (branch_struct is not None), ('metadata not in this part of the xml tree')
        return branch_struct
