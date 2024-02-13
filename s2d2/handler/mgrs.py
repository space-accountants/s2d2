#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import os

import numpy as np
import geopandas as gpd

from fiona.drvsupport import supported_drivers
from osgeo import osr
from shapely import wkt
from shapely.geometry import Polygon, MultiPolygon

from ..checking.naming import check_mgrs_code
from ..mapping_tools import ll2map, get_utm_zone

supported_drivers['KML'] = 'rw'

MGRS_TILING_URL = ("https://sentinels.copernicus.eu/documents/247904/1955685/"
                   "S2A_OPER_GIP_TILPAR_MPC__20151209T095117_V20150622T000000"
                   "_21000101T000000_B00.kml")

MGRS_TILING_FILENAME = 'sentinel2_tiles_world.geojson'
MGRS_TILING_DIR_DEFAULT = os.path.join('.', 'data', 'MGRS')

def _kml_to_gdf(tile_path):
    """
    Read MGRS kml file as a GeoPandas DataFrame, keep only polygon geometries.

    Parameters
    ----------
    tile_path : str
        kml file name

    Returns
    -------
    geopandas.geodataframe.GeoDataFrame
        KML file read as a GeoDataFrame
    """
    # Load kml file
    gdf = gpd.read_file(tile_path, driver="KML")

    # Drop Description, whick is KML specific for visualization
    gdf = gdf.drop(columns=['Description'])

    # Unpack geometries
    gdf_exploded = gdf.explode(index_parts=False)

    # Select all entries which are Polygon
    mask = gdf_exploded["geometry"].apply(lambda x: isinstance(x, Polygon))
    gdf_out = gdf_exploded.loc[mask]

    return gdf_out

def _get_mgrs_abc():
    mgrsABC = [chr(i) for i in list(range(65,73)) +
               list(range(74,79)) + list(range(80,91))]
    return mgrsABC

def _mgrs_to_search_geometry(tile_code):
    """
    Get a geometry from the tile code. The search geometry is a 6-deg longitude
    stripe, augmented with the stripe across the antimeridian for the tiles
    that are crossing this line.

    Parameters
    ----------
    tile_code : str
        MGRS code of the tile

    Returns
    -------
    shapely.Geometry
        geometry to restrict the search area
    """
    nr_lon = int(tile_code[0:2])  # first two letters indicates longitude range
    min_lon = -180. + (nr_lon - 1) * 6.
    max_lon = -180. + nr_lon * 6.
    geom = Polygon.from_bounds(min_lon, -90.0, max_lon, 90.0)
    extra = None
    if nr_lon == 1:
        extra = Polygon.from_bounds(174, -90.0, 180, 90.0)
    elif nr_lon == 60:
        extra = Polygon.from_bounds(-180, -90.0, -174, 90.0)
    if extra is not None:
        geom = MultiPolygon(polygons=[geom, extra])
    return geom

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

def get_geom_from_tile_code(tile_code, tile_path=None):
    """
    Get the geometry of a certain MGRS tile

    Parameters
    ----------
    tile_code : string
        MGRS tile coding, e.g.: '05VMG'
    tile_path : string
        Path to the geometric metadata

    Returns
    -------
    shapely.geometry.polygon.Polygon
        Geometry of the MGRS tile, in lat/lon

    See Also
    --------
    .get_tile_code_from_geom, .get_bbox_from_tile_code
    """
    if tile_path is None:
        tile_path = os.path.join(MGRS_TILING_DIR_DEFAULT, MGRS_TILING_FILENAME)

    tile_code = check_mgrs_code(tile_code)

    # Derive a search geometry from the tile code
    search_geom = _mgrs_to_search_geometry(tile_code)

    # Load tiles
    mgrs_tiles = gpd.read_file(tile_path, mask=search_geom)

    geom = mgrs_tiles[mgrs_tiles['Name'] == tile_code]["geometry"]

    if len(geom) == 0:
        raise ValueError('MGRS tile code does not seem to exist')

    return geom.unary_union

def get_tile_codes_from_geom(geom, tile_path=None):
    """
    Get the codes of the MGRS tiles intersecting a given geometry

    Parameters
    ----------
    geom : {shapely.geometry, string, dict, GeoDataFrame, GeoSeries}
        geometry object with the given dict-like geojson geometry, GeoSeries,
        GeoDataFrame, shapely geometry or well known text, i.e.:
        'POLYGON ((x y, x y, x y))'
    tile_path : string
        Path to the geometric metadata

    Returns
    -------
    tuple
        MGRS tile codes

    See Also
    --------
    .get_geom_from_tile_code, .get_bbox_from_tile_code
    """

    if tile_path is None:
        tile_path = os.path.join(MGRS_TILING_DIR_DEFAULT, MGRS_TILING_FILENAME)

    # If a wkt str, convert to shapely geometry
    if isinstance(geom, str):
        geom = wkt.loads(geom)

    # Uniform CRS
    if isinstance(geom, gpd.GeoSeries) or isinstance(geom, gpd.GeoDataFrame):
        example = gpd.read_file(tile_path, rows=1)
        geom = geom.set_crs(example.crs)

    # Load tiles intersects the search box
    tiles = gpd.read_file(tile_path, mask=geom)

    # Get the codes in tuple
    codes = tuple(tiles["Name"])

    return codes

def get_bbox_from_tile_code(tile_code, tile_path=None):
    """
    Get the bounds of a certain MGRS tile

    Parameters
    ----------
    tile_code : string, e.g.: '05VMG'
        MGRS tile coding
    tile_path : string
        Path to the geometric metadata

    Returns
    -------
    numpy.ndarray, size=(1,4), dtype=float
        bounding box, in the following order: min max X, min max Y
    """

    geom = get_geom_from_tile_code(tile_code, tile_path=tile_path)

    toi = geom.bounds
    bbox = np.array([toi[0], toi[2], toi[1], toi[3]])
    return bbox

