import numpy as np

from .checking.mapping import correct_geoTransform
from .checking.array import are_two_arrays_equal, correct_floating_parameter

def pix2map(geoTransform, i, j):
    """ transform local image coordinates to map coordinates

    Parameters
    ----------
    geoTransform : tuple, size={(6,1), (8,1)}
        georeference transform of an image.
    i : np.array, ndim={1,2,3}, dtype=float
        row coordinate(s) in local image space
    j : np.array, ndim={1,2,3}, dtype=float
        column coordinate(s) in local image space

    Returns
    -------
    x : np.array, size=(m), ndim={1,2,3}, dtype=float
        horizontal map coordinate.
    y : np.array, size=(m), ndim={1,2,3}, dtype=float
        vertical map coordinate.

    See Also
    --------
    map2pix

    Notes
    -----
    Two different coordinate system are used here:

        .. code-block:: text

          indexing   |           indexing    ^ y
          system 'ij'|           system 'xy' |
                     |                       |
                     |       i               |       x
             --------+-------->      --------+-------->
                     |                       |
                     |                       |
          image      | j         map         |
          based      v           based       |

    """
    geoTransform = correct_geoTransform(geoTransform)
    if type(i) in (np.ma.core.MaskedArray, np.ndarray):
        are_two_arrays_equal(i, j)
    else:  # if only a float is given
        i, j = correct_floating_parameter(i), correct_floating_parameter(j)

    x = geoTransform[0] + \
        np.multiply(geoTransform[1], j) + np.multiply(geoTransform[2], i)

    y = geoTransform[3] + \
        np.multiply(geoTransform[4], j) + np.multiply(geoTransform[5], i)
    return x, y


def map2pix(geoTransform, x, y):
    """ transform map coordinates to local image coordinates

    Parameters
    ----------
    geoTransform : tuple, size={(6,1), (8,1)}
        georeference transform of an image.
    x : np.array, size=(m), ndim={1,2,3}, dtype=float
        horizontal map coordinate.
    y : np.array, size=(m), ndim={1,2,3}, dtype=float
        vertical map coordinate.

    Returns
    -------
    i : np.array, ndim={1,2,3}, dtype=float
        row coordinate(s) in local image space
    j : np.array, ndim={1,2,3}, dtype=float
        column coordinate(s) in local image space

    See Also
    --------
    pix2map

    Notes
    -----
    Two different coordinate system are used here:

        .. code-block:: text

          indexing   |           indexing    ^ y
          system 'ij'|           system 'xy' |
                     |                       |
                     |       i               |       x
             --------+-------->      --------+-------->
                     |                       |
                     |                       |
          image      | j         map         |
          based      v           based       |

    """
    geoTransform = correct_geoTransform(geoTransform)
    if type(x) in (np.ma.core.MaskedArray, np.ndarray):
        are_two_arrays_equal(x, y)
    else:  # if only a float is given
        x, y = correct_floating_parameter(x), correct_floating_parameter(y)

    A = np.array(geoTransform[:-2]).reshape(2, 3)[:, 1:]
    A_inv = np.linalg.inv(A)

    x_loc = x - geoTransform[0]
    y_loc = y - geoTransform[3]

    j = np.multiply(x_loc, A_inv[0, 0]) + np.multiply(y_loc, A_inv[0, 1])
    i = np.multiply(x_loc, A_inv[1, 0]) + np.multiply(y_loc, A_inv[1, 1])
    return i, j

def pix_centers(geoTransform, rows=None, cols=None, make_grid=True):
    """ provide the pixel coordinate from the axis, or the whole grid

    Parameters
    ----------
    geoTransform : tuple, size={(6,), (8,)}
        georeference transform of an image.
    rows : integer, {x ∈ ℕ | x ≥ 0}, default=None
        amount of rows in an image.
    cols : integer, {x ∈ ℕ | x ≥ 0}, default=None
        amount of collumns in an image.
    make_grid : bool, optional
        Should a grid be made. The default is True.

    Returns
    -------
    X : np.array, dtype=float
         * "make_grid" : size=(m,n)
         otherwise : size=(1,n)
    Y : np.array, dtype=float
         * "make_grid" : size=(m,n)
         otherwise : size=(m,1)

    See Also
    --------
    map2pix, pix2map

    Notes
    -----
    Two different coordinate system are used here:

        .. code-block:: text

          indexing   |           indexing    ^ y
          system 'ij'|           system 'xy' |
                     |                       |
                     |       j               |       x
             --------+-------->      --------+-------->
                     |                       |
                     |                       |
          image      | i         map         |
          based      v           based       |

    """
    geoTransform = correct_geoTransform(geoTransform)
    if rows is None:
        assert len(geoTransform) == 8, (
            'please provide the dimensions of the ' +
            'imagery, or have this included in the ' + 'geoTransform.')
        rows, cols = int(geoTransform[-2]), int(geoTransform[-1])
    i, j = np.linspace(0, rows - 1, rows), np.linspace(0, cols - 1, cols)

    if make_grid:
        jj, ii = np.meshgrid(j, i)
        X, Y = pix2map(geoTransform, ii, jj)
        return X, Y
    else:
        x, y_dummy = pix2map(geoTransform, np.repeat(i[0], len(j)), j)
        x_dummy, y = pix2map(geoTransform, i, np.repeat(j[0], len(i)))
        return x, y
