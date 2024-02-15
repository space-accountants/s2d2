#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from s2d2.checking.mapping import correct_geotransform

def are_two_arrays_equal(a, b):
    """ check if two arrays have the same dimensions

    Parameters
    ----------
    a, b : {np.ma.array, np.ndarray}
        arrays of interest

    Returns
    -------
    bool, provides True if all sizes are equal
    """
    no_arr_str = 'please provide an array'
    assert type(a) in (np.ma.core.MaskedArray, np.ndarray), no_arr_str
    assert type(b) in (np.ma.core.MaskedArray, np.ndarray), no_arr_str
    assert (a.ndim == b.ndim), ('please provide arrays of same dimension')
    assert (a.shape == b.shape), ('please provide arrays of equal shape')

def are_three_arrays_equal(a, b, c):
    """ check if three arrays have the same dimensions

    Parameters
    ----------
    a, b, c : {np.ma.array, np.ndarray}
        arrays of interest

    Returns
    -------
    bool, provides True if all sizes are equal
    """
    no_arr_str = 'please provide an array'
    no_size_str = 'please provide arrays of the same size'
    assert type(a) in (np.ma.core.MaskedArray, np.ndarray), no_arr_str
    assert type(b) in (np.ma.core.MaskedArray, np.ndarray), no_arr_str
    assert type(c) in (np.ma.core.MaskedArray, np.ndarray), no_arr_str
    assert len(set({a.shape[0], b.shape[0], c.shape[0]})) == 1, no_size_str
    if np.logical_and.reduce((a.ndim>1, b.ndim>1, c.ndim>1)):
        assert len(set({a.shape[1], b.shape[1], c.shape[1]})) == 1, no_size_str
    assert len(set({a.ndim, b.ndim, c.ndim})) == 1, no_size_str

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
        return float(a)

def make_same_size(old,
                   geotransform_old,
                   geotransform_new,
                   rows_new=None,
                   cols_new=None):
    """ clip array to the same size as another array

    Parameters
    ----------
    old : np.array, size=(m,n), dtype={float, complex}
        data array to be clipped.
    geotransform_old : tuple, size={(6,), (8,)}
        georeference transform of the old image.
    geotransform_new : tuple, size={(6,), (8,)}
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
    geotransform_old = correct_geotransform(geotransform_old)
    geotransform_new = correct_geotransform(geotransform_new)

    if len(geotransform_new) == 8:
        rows_new, cols_new = geotransform_new[-2], geotransform_new[-1]

    # look at upper left coordinate
    dj = np.round(
        (geotransform_new[0] - geotransform_old[0]) / geotransform_new[1])
    di = np.round(
        (geotransform_new[3] - geotransform_old[3]) / geotransform_new[1])

    if np.sign(dj) == -1:  # extend array by simple copy of border values
        old = np.concatenate((np.repeat(
            np.expand_dims(old[:, 0], axis=1), abs(dj), axis=1), old),
                             axis=1)
    elif np.sign(dj) == 1:  # reduce array
        old = old[:, abs(dj).astype(int):]

    if np.sign(di) == -1:  # reduce array
        old = old[abs(di).astype(int):, :]
    elif np.sign(di) == 1:  # extend array by simple copy of border values
        old = np.concatenate((np.repeat(
            np.expand_dims(old[0, :], axis=1).T, abs(di), axis=0), old),
                             axis=0)

    # as they are now alligned, look at the lower right corner
    di, dj = rows_new - old.shape[0], cols_new - old.shape[1]

    if np.sign(dj) == -1:  # reduce array
        old = old[:, :dj]
    elif np.sign(dj) == 1:  # extend array by simple copy of border values
        old = np.concatenate((np.repeat(
            old, np.expand_dims(old[:, -1], axis=1), abs(dj), axis=1)),
                             axis=1)

    if np.sign(di) == -1:  # reduce array
        old = old[:di, :]
    elif np.sign(di) == 1:  # extend array by simple copy of border values
        old = np.concatenate((np.repeat(
            old, np.expand_dims(old[-1, :], axis=1).T, abs(di), axis=0)),
                             axis=0)

    new = old
    return new
