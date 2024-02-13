#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def are_two_arrays_equal(A, B):
    """ check if two arrays have the same dimensions
    """
    assert type(A) in (np.ma.core.MaskedArray, np.ndarray), \
        ('please provide an array')
    assert type(B) in (np.ma.core.MaskedArray, np.ndarray), \
        ('please provide an array')
    assert A.ndim==B.ndim, ('please provide arrays of same dimension')
    assert A.shape==B.shape, ('please provide arrays of equal shape')
    return

def are_three_arrays_equal(A, B, C):
    """ check if three arrays have the same dimensions
    """
    assert type(A) in (np.ma.core.MaskedArray, np.ndarray), \
        ('please provide an array')
    assert type(B) in (np.ma.core.MaskedArray, np.ndarray), \
        ('please provide an array')
    assert type(C) in (np.ma.core.MaskedArray, np.ndarray), \
        ('please provide an array')
    assert len(set({A.shape[0], B.shape[0], C.shape[0]})) == 1, \
         ('please provide arrays of the same size')
    if np.logical_and.reduce((A.ndim>1, B.ndim>1, C.ndim>1)):
        assert len(set({A.shape[1], B.shape[1], C.shape[1]})) == 1, \
             ('please provide arrays of the same size')
    assert len(set({A.ndim, B.ndim, C.ndim})) == 1, \
         ('please provide arrays of the same dimension')
    return

def correct_floating_parameter(a):
    """ sometimes a float is asked for, but this is given in a list, array or
    tuple. This function provides the float that is within
    """
    if type(a) in (np.ma.core.MaskedArray, np.ndarray):
        assert a.size == 1, 'please provide one parameter'
        a = a[0]
    if type(a) in (list, tuple):
        assert len(a) == 1, 'please provide one parameter'
        a = a[0]

    assert isinstance(a, (int, float)), 'please provide an integer'
    if isinstance(a, int):
        a = float(a)
    return a
