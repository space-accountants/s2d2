import pytest
import urllib.request
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from shapely import wkt

from s2d2.checking.naming import check_mgrs_code
from s2d2.handler import mgrs

# A subset of original MGRS, sampling every 1500 entries
TESTDATA_KML = 'testdata/MGRS/mgrs_tiles_tiny.kml'

# The above kml in geojson
TESTDATA_GEOJSON = 'testdata/MGRS/mgrs_tiles_tiny.geojson'

# Example tiles to query
TILE_CODE = '51GWL'  # tile code
POINT_WITHIN = 'POINT (123.69324 -44.74430)'  # A point within the tile
POINT_OUTSIDE = 'POINT (-53.1 1.5)'  # A point outside tile
# A polygon within the tile
POLYGON_WITHIN = 'POLYGON ((123.199 -44.451, 124.179 -44.445, 124.193 -45.034, 123.199 -45.040, 123.199 -44.451))'
# A polygon overlaps but larger than the tile
POLYGON_LARGER = 'POLYGON ((0.0 -10, 125.0 -10, 125.0 -45, 0.0 -45, 0.0 -10))'
# A ploygon outside
POLYGON_OUTSIDE = 'POLYGON ((0.0 0.0, 1.0 0.0, 1.0 1.0, 0.0 1.0, 0.0 0.0))'


def geom_types(wkt_text):
    """
        Make all types of geom from wkt
    """
    return [wkt_text,
            wkt.loads(wkt_text),
            gpd.GeoSeries.from_wkt([wkt_text]),
            gpd.GeoDataFrame(geometry=gpd.GeoSeries.from_wkt([wkt_text]))]


# Data retrieval tests
def test_mgrs_kml_url_valid():
    assert urllib.request.urlopen(mgrs.MGRS_TILING_URL).code == 200


def test_kml_to_gdf_output_all_polygon():
    gdf = mgrs._kml_to_gdf(TESTDATA_KML)
    assert all([isinstance(geom, Polygon) for geom in gdf["geometry"]])


def test_kml_to_gdf_output_same_crs():
    gdf_kml = gpd.read_file(TESTDATA_KML, driver="KML")
    gdf = mgrs._kml_to_gdf(TESTDATA_KML)
    assert gdf.crs.equals(gdf_kml.crs)


# MGRS query tests
def test_normalize_mgrs_code_not_str():
    with pytest.raises(TypeError):
        check_mgrs_code(123)


def test_normalize_mgrs_code_wrong_code():
    with pytest.raises(ValueError):
        check_mgrs_code('WRONG')


def test_mgrs_to_search_geometry_range():
    geom = mgrs._mgrs_to_search_geometry(TILE_CODE)
    assert geom == Polygon.from_bounds(120.0, -90.0, 126.0, 90.0)


def test_get_bbox_from_tile_code_bbox_equal():
    bbox = mgrs.get_bbox_from_tile_code(
        TILE_CODE, tile_path=TESTDATA_GEOJSON)
    assert bbox == pytest.approx(
        np.array([122.999, 124.398, -45.241, -44.244]), 0.001)


def test_get_geom_from_tile_code_geometry():
    # geom from query equals to directly read geom
    geom = mgrs.get_geom_from_tile_code(
        TILE_CODE, tile_path=TESTDATA_GEOJSON)
    gdf = gpd.read_file(TESTDATA_GEOJSON)
    geom_read = gdf.loc[gdf["Name"] == TILE_CODE]["geometry"]
    assert geom.equals(geom_read.squeeze())


def test_get_geom_from_tile_code_crossing_antimeridian():
    # geom of tiles crossing the antimeridian should be multipolygons
    geom = mgrs.get_geom_from_tile_code(
        "01CCV", tile_path=TESTDATA_GEOJSON)
    assert isinstance(geom, MultiPolygon)


@pytest.mark.parametrize("geom", geom_types(POINT_WITHIN))
def test_get_tile_codes_from_geom_point_within(geom):
    code = mgrs.get_tile_codes_from_geom(
        geom, tile_path=TESTDATA_GEOJSON)
    assert len(code) == 1
    assert code[0] == TILE_CODE


@pytest.mark.parametrize("geom", geom_types(POINT_OUTSIDE))
def test_get_tile_codes_from_geom_point_outside(geom):
    code = mgrs.get_tile_codes_from_geom(
        geom, tile_path=TESTDATA_GEOJSON)
    assert len(code) == 0


@pytest.mark.parametrize("geom", geom_types(POLYGON_WITHIN))
def test_get_tile_codes_from_geom_polygon_within(geom):
    code = mgrs.get_tile_codes_from_geom(
        geom, tile_path=TESTDATA_GEOJSON)
    assert len(code) == 1
    assert code[0] == TILE_CODE


@pytest.mark.parametrize("geom", geom_types(POLYGON_LARGER))
def test_get_tile_codes_from_geom_polygon_overlap(geom):
    code = mgrs.get_tile_codes_from_geom(
        geom, tile_path=TESTDATA_GEOJSON)
    assert len(code) > 1
    assert TILE_CODE in code


@pytest.mark.parametrize("geom", geom_types(POLYGON_OUTSIDE))
def test_get_tile_codes_from_geom_polygon_outside(geom):
    code = mgrs.get_tile_codes_from_geom(
        geom, tile_path=TESTDATA_GEOJSON)
    assert len(code) == 0
