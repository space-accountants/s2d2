import numpy as np

from osgeo import osr

from .checking.mapping import correct_geotransform
from .checking.array import are_two_arrays_equal, correct_floating_parameter

def pix2map(geotransform, i, j):
    """ transform local image coordinates to map coordinates

    Parameters
    ----------
    geotransform : tuple, size={(6,1), (8,1)}
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
    geotransform = correct_geotransform(geotransform)
    if type(i) in (np.ma.core.MaskedArray, np.ndarray):
        are_two_arrays_equal(i, j)
    else:  # if only a float is given
        i, j = correct_floating_parameter(i), correct_floating_parameter(j)

    x = geotransform[0] + \
        np.multiply(geotransform[1], j) + np.multiply(geotransform[2], i)

    y = geotransform[3] + \
        np.multiply(geotransform[4], j) + np.multiply(geotransform[5], i)
    return x, y


def map2pix(geotransform, x, y):
    """ transform map coordinates to local image coordinates

    Parameters
    ----------
    geotransform : tuple, size={(6,1), (8,1)}
        georeference transform of an image.
    x : np.array, size=(m), ndim={1,2,3}, dtype=float
        horizontal map coordinate.
    y : np.array, size=(m), ndim={1,2,3}, dtype=float
        vertical map coordinate.

    Returns
    -------
    i : np.ndarray, ndim={1,2,3}, dtype=float
        row coordinate(s) in local image space
    j : np.ndarray, ndim={1,2,3}, dtype=float
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
    geotransform = correct_geotransform(geotransform)
    if type(x) in (np.ma.core.MaskedArray, np.ndarray):
        are_two_arrays_equal(x, y)
    else:  # if only a float is given
        x, y = correct_floating_parameter(x), correct_floating_parameter(y)

    A = np.array(geotransform[:-2]).reshape(2, 3)[:, 1:]
    P = np.linalg.inv(A)

    x_loc = x - geotransform[0]
    y_loc = y - geotransform[3]

    j = np.multiply(x_loc, P[0, 0]) + np.multiply(y_loc, P[0, 1])
    i = np.multiply(x_loc, P[1, 0]) + np.multiply(y_loc, P[1, 1])
    return i, j

def pix_centers(geotransform, rows=None, cols=None, make_grid=True):
    """ provide the pixel coordinate from the axis, or the whole grid

    Parameters
    ----------
    geotransform : tuple, size={(6,), (8,)}
        georeference transform of an image.
    rows : integer, {x ∈ ℕ | x ≥ 0}, default=None
        amount of rows in an image.
    cols : integer, {x ∈ ℕ | x ≥ 0}, default=None
        amount of collumns in an image.
    make_grid : bool, optional
        Should a grid be made. The default is True.

    Returns
    -------
    X : np.ndarray, dtype=float
         * "make_grid" : size=(m,n)
         * otherwise : size=(1,n)
    Y : np.ndarray, dtype=float
         * "make_grid" : size=(m,n)
         * otherwise : size=(m,1)

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
    geotransform = correct_geotransform(geotransform)
    if rows is None:
        assert len(geotransform) == 8, (
            'please provide the dimensions of the ' +
            'imagery, or have this included in the ' + 'geotransform.')
        rows, cols = int(geotransform[-2]), int(geotransform[-1])
    i, j = np.linspace(0, rows - 1, rows), np.linspace(0, cols - 1, cols)

    if make_grid:
        jj, ii = np.meshgrid(j, i)
        X, Y = pix2map(geotransform, ii, jj)
        return X, Y
    x, _ = pix2map(geotransform, np.repeat(i[0], len(j)), j)
    _, y = pix2map(geotransform, i, np.repeat(j[0], len(i)))
    return x, y

def get_bbox(geotransform, rows=None, cols=None):
    """ given array meta data, calculate the bounding box

    Parameters
    ----------
    geotransform : tuple, size=(6,)
        georeference transform of an image.
    rows : integer, {x ∈ ℕ | x ≥ 0}
        amount of rows in an image.
    cols : integer, {x ∈ ℕ | x ≥ 0}
        amount of collumns in an image.

    Returns
    -------
    bbox : np.ndarray, size=(4,), dtype=float
        bounding box, in the following order: min max X, min max Y

    See Also
    --------
    map2pix, pix2map, pix_centers

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
    geotransform = correct_geotransform(geotransform)
    if rows is None:
        assert len(geotransform) >= 8, ('please provide raster information')
        rows, cols = geotransform[6], geotransform[7]

    X = geotransform[0] + \
        np.array([0, cols]) * geotransform[1] + np.array([0, rows]) * \
        geotransform[2]

    Y = geotransform[3] + \
        np.array([0, cols]) * geotransform[4] + np.array([0, rows]) * \
        geotransform[5]

    bbox = np.hstack((np.sort(X), np.sort(Y)))
    return bbox

def create_local_crs():
    """ create spatial refence of local horizontal datum

    Returns
    -------
    crs : {osr.SpatialReference, string}
        the coordinate reference system via GDAL SpatialReference description

    See Also
    --------
    s2d2.checking.mapping.is_crs_an_srs
    """
    crs = osr.SpatialReference()
    crs.ImportFromEPSG(8377)

    return crs
