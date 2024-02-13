#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from s2d2.checking.mapping import correct_geoTransform

def are_two_arrays_equal(A, B):
    """ check if two arrays have the same dimensions

    Parameters
    ----------
    A, B : {np.ma.array, np.ndarray}
        arrays of interest

    Returns
    -------
    bool, provides True if all sizes are equal
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

    Parameters
    ----------
    A, B, C : {np.ma.array, np.ndarray}
        arrays of interest

    Returns
    -------
    bool, provides True if all sizes are equal
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

    Parameters
    ----------
    a : {np.ma.array, np.ndarray, tuple, list, float, int}
        float of interest

    Returns
    -------
    float, int
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

def make_same_size(Old,
                   geoTransform_old,
                   geoTransform_new,
                   rows_new=None,
                   cols_new=None):
    """ clip array to the same size as another array

    Parameters
    ----------
    Old : np.array, size=(m,n), dtype={float, complex}
        data array to be clipped.
    geoTransform_old : tuple, size={(6,), (8,)}
        georeference transform of the old image.
    geoTransform_new : tuple, size={(6,), (8,)}
        georeference transform of the new image.
    rows_new : integer, {x ∈ ℕ | x ≥ 0}
        amount of rows of the new image.
    cols_new : integer, {x ∈ ℕ | x ≥ 0}
        amount of collumns of the new image.

    Returns
    -------
    New : np.array, size=(k,l), dtype={float,complex}
        clipped data array.
    """
    geoTransform_old = correct_geoTransform(geoTransform_old)
    geoTransform_new = correct_geoTransform(geoTransform_new)

    if len(geoTransform_new) == 8:
        rows_new, cols_new = geoTransform_new[-2], geoTransform_new[-1]

    # look at upper left coordinate
    dj = np.round(
        (geoTransform_new[0] - geoTransform_old[0]) / geoTransform_new[1])
    di = np.round(
        (geoTransform_new[3] - geoTransform_old[3]) / geoTransform_new[1])

    if np.sign(dj) == -1:  # extend array by simple copy of border values
        Old = np.concatenate((np.repeat(
            np.expand_dims(Old[:, 0], axis=1), abs(dj), axis=1), Old),
                             axis=1)
    elif np.sign(dj) == 1:  # reduce array
        Old = Old[:, abs(dj).astype(int):]

    if np.sign(di) == -1:  # reduce array
        Old = Old[abs(di).astype(int):, :]
    elif np.sign(di) == 1:  # extend array by simple copy of border values
        Old = np.concatenate((np.repeat(
            np.expand_dims(Old[0, :], axis=1).T, abs(di), axis=0), Old),
                             axis=0)

    # as they are now alligned, look at the lower right corner
    di, dj = rows_new - Old.shape[0], cols_new - Old.shape[1]

    if np.sign(dj) == -1:  # reduce array
        Old = Old[:, :dj]
    elif np.sign(dj) == 1:  # extend array by simple copy of border values
        Old = np.concatenate((np.repeat(
            Old, np.expand_dims(Old[:, -1], axis=1), abs(dj), axis=1)),
                             axis=1)

    if np.sign(di) == -1:  # reduce array
        Old = Old[:di, :]
    elif np.sign(di) == 1:  # extend array by simple copy of border values
        Old = np.concatenate((np.repeat(
            Old, np.expand_dims(Old[-1, :], axis=1).T, abs(di), axis=0)),
                             axis=0)

    New = Old
    return New
