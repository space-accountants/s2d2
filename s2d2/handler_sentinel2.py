#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Space Accountants"
__license__ = "MIT License - You must cite this source"
__version__ = "202309"
__maintainer__ = "B. Altena"
__email__ = "info at space hyphen accountants dot eu"

import os

from osgeo import osr

import numpy as np
import pandas as pd
import pyproj
import shapely.ops

from handler_mgrs import get_geom_for_tile_code
from naming_check import check_mgrs_code
from handler_xml import get_root_of_table
from dhdt.generic.handler_dat import get_list_files


def get_s2_dict(s2_df):
    """ given a dataframe with a filename within, create a dictionary with scene
    and satellite specific metadata.

    Parameters
    ----------
    s2_df : pd.dataframe
        metadata and general multi spectral information about the MSI
        instrument that is onboard Sentinel-2

    Returns
    -------
    s2_dict : dictionary
        dictionary with scene specific meta data, giving the following
        information:
        * 'COSPAR' : string
            id of the satellite
        * 'NORAD' : string
            id of the satellite
        * 'instruments' : {'MSI'}
            specific instrument
        * 'launch_date': string, e.g.:'+2015-06-23'
            date when the satellite was launched into orbit
        * 'orbit': string
            specifies the type of orbit
        * 'full_path': string
            the location where the root folder is situated
        * 'MTD_DS_path': string
            the location where metadata about the datastrip is present
        * 'MTD_TL_path': string
            the location where metadata about the tile is situated
        * 'QI_DATA_path': string
            the location where metadata about detector readings are present
        * 'IMG_DATA_path': string
            the location where the imagery is present
        * 'date': string, e.g.: '+2019-10-25'
            date of acquisition
        * 'tile_code': string, e.g.: '05VMG'
            MGRS tile coding
        * 'relative_orbit': integer, {x ∈ ℕ}
            orbit from which the imagery were taken

    Notes
    -----
    The metadata is scattered over a folder structure for Sentinel-2, it has
    typically the following set-up:

    .. code-block:: text

        * S2X_MSIL1C_20XX...   <- path is given in s2_dict["full_path"]
        ├ AUX_DATA
        ├ DATASTRIP
        │  └ DS_XXX_XXXX...    <- path is given in s2_dict["MTD_DS_path"]
        │     └ QI_DATA
        │        └ MTD_DS.xml  <- metadata about the data-strip
        ├ GRANULE
        │  └ L1C_TXXXX_XXXX... <- path is given in s2_dict["MTD_TL_path"]
        │     ├ AUX_DATA
        │     ├ IMG_DATA
        │     │  ├ TXXX_X..XXX.jp2
        │     │  └ TXXX_X..XXX.jp2
        │     ├ QI_DATA
        │     └ MTD_TL.xml     <- metadata about the tile
        ├ HTML
        ├ rep_info
        ├ manifest.safe
        ├ INSPIRE.xml
        └ MTD_MSIL1C.xml       <- metadata about the product

    The following acronyms are used:

    - AUX : auxiliary
    - COSPAR : committee on space research international designator
    - DS : datastrip
    - IMG : imagery
    - L1C : product specification,i.e.: level 1, processing step C
    - MGRS : military grid reference system
    - MTD : metadata
    - MSI : multi spectral instrument
    - NORAD : north american aerospace defense satellite catalog number
    - TL : tile
    - QI : quality information
    - s2 : Sentinel-2

    See Also
    --------
    dhdt.generic.get_s2_image_locations

    """
    assert isinstance(s2_df, pd.DataFrame), ('please provide a dataframe')
    assert 'filepath' in s2_df, ('please first run "get_s2_image_locations"'+
                                ' to find the proper file locations')

    # search for folder with .SAFE
    im_path = os.path.dirname(s2_df.filepath[0])
    safe_folder = _get_safe_foldername(im_path)
    if safe_folder is not None:
        if safe_folder.startswith('S2A'):
            s2_dict = list_platform_metadata_s2a()
        else:
            s2_dict = list_platform_metadata_s2b()
        s2_dict = _get_safe_structure_s2(im_path, s2_dict=s2_dict)
    elif len(get_list_files(im_path, '.json'))!=0:
        s2_dict = _get_stac_structure_s2(im_path)
    return s2_dict

def _get_safe_foldername(im_path):
    safe_folder = None
    for x in im_path.split(os.sep):
        if x.endswith('SAFE'):
            safe_folder = x
            break
    return safe_folder

def _get_safe_path(im_path):
    safe_path = ''
    for x in im_path.split(os.sep):
        safe_path += x
        safe_path += os.sep
        if x.endswith('SAFE'):
            break
    return safe_path

def _get_safe_structure_s2(im_path, s2_dict=None):
    safe_path = _get_safe_path(im_path)
    ds_folder = list(filter(lambda x: x.startswith('DS'),
                            os.listdir(os.path.join(safe_path,
                                                    'DATASTRIP'))))[0]
    tl_path = os.sep.join(im_path.split(os.sep)[:-1])
    ds_path = os.path.join(safe_path, 'DATASTRIP', ds_folder)
    assert os.path.isfile(os.path.join(ds_path, 'MTD_DS.xml')), \
        'please make sure MTD_DS.xml is present in the directory'
    assert os.path.isfile(os.path.join(tl_path, 'MTD_TL.xml')), \
        'please make sure MTD_TL.xml is present in the directory'

    safe_folder = _get_safe_foldername(im_path)
    s2_time, s2_orbit, s2_tile = meta_s2string(safe_folder)
    s2_dict.update({'full_path': safe_path,
                    'MTD_DS_path': ds_path,
                    'MTD_TL_path': tl_path,
                    'IMG_DATA_path': os.path.join(tl_path, 'IMG_DATA'),
                    'QI_DATA_path': os.path.join(tl_path, 'QI_DATA'),
                    'date': s2_time,
                    'tile_code': s2_tile[1:],
                    'relative_orbit': int(s2_orbit[1:])})
    return s2_dict

def _get_stac_structure_s2(im_path, s2_dict=None):
    # make sure metadata files are present
    assert os.path.isfile(os.path.join(im_path, 'MTD_DS.xml')), \
        'please make sure MTD_DS.xml is present in the directory'
    assert os.path.isfile(os.path.join(im_path, 'MTD_TL.xml')), \
        'please make sure MTD_TL.xml is present in the directory'

    if s2_dict is None:
        json_file = get_list_files(im_path, '.json')[0]
        # Planet has a generic json-file
        if json_file.find('metadata')!=0:
            if json_file[:3] in 'S2A':
                s2_dict = list_platform_metadata_s2a()
            else:
                s2_dict = list_platform_metadata_s2b()
            s2_time, s2_orbit, s2_tile = meta_s2string(json_file)
        else: # look into the MTD_TL?x
            im_list = get_list_files(im_path, '.jp2') + \
                      get_list_files(im_path, '.tif')
            s2_time, s2_orbit, s2_tile = meta_s2string(im_list[0])
    s2_dict.update({'full_path': im_path,
                    'MTD_DS_path': im_path,
                    'MTD_TL_path': im_path,
                    'IMG_DATA_path': im_path,
                    'date': s2_time,
                    'tile_code': s2_tile[1:]})
    if s2_orbit is not None:
        s2_dict.update({'relative_orbit': int(s2_orbit[1:])})
    return s2_dict

def list_platform_metadata_s2a():
    s2a_dict = {
        'COSPAR': '2015-028A',
        'NORAD': 40697,
        'instruments': {'MSI'},
        'satellite': 'sentinel-2',
        'constellation': 'sentinel',
        'launch_date': '+2015-06-23',
        'orbit': 'sso',
        'mass': 1129.541, # [kg]
        'inclination': 98.5621, # https://www.n2yo.com/satellite/?s=40697
        'revolutions_per_day': 14.30824258387262,
        'J': np.array([[558, 30, -30],[30, 819, 30],[-30, 30, 1055]])}
    return s2a_dict

def list_platform_metadata_s2b():
    s2b_dict = {
        'COSPAR': '2017-013A',
        'NORAD': 42063,
        'instruments': {'MSI'},
        'satellite': 'sentinel-2',
        'constellation': 'sentinel',
        'launch_date': '+2017-03-07',
        'orbit': 'sso',
        'inclination': 98.5664,
        'revolutions_per_day': 14.30818491298178,
        'J': np.array([[558, 30, -30],[30, 819, 30],[-30, 30, 1055]])}
    return s2b_dict

def get_generic_s2_raster(tile_code, spac=10, tile_path=None):
    """
    Create spatial metadata of a Sentinel-2, so no downloading is needed.

    Parameters
    ----------
    tile_code : string
        mgrs tile coding, which is also given in the filename ("TXXXXX")
    spac : integer, {10,20,60}, unit=m
        pixel spacing of the raster
    tile_path : string
        path to the MGRS tiling file

    Returns
    -------
    tuple
        georeference transform of an image, including shape
    pyproj.crs.crs.CRS
        coordinate reference system (CRS)

    Notes
    -----
    The tile structure is a follows "AABCC"
        * "AA" utm zone number, starting from the East, with steps of 8 degrees
        * "B" latitude zone, starting from the South, with steps of 6 degrees

    The following acronyms are used:

    - CRS : coordinate reference system
    - MGRS : US military grid reference system
    - s2 : Sentinel-2
    """
    assert spac in (10, 20, 60,), 'please provide correct pixel resolution'

    tile_code = check_mgrs_code(tile_code)
    geom = get_geom_for_tile_code(tile_code, tile_path=tile_path)

    # specify coordinate systems
    crs = pyproj.CRS.from_epsg(4326)  # lat lon
    epsg = get_epsg_from_mgrs_tile(tile_code)
    crs_utm = pyproj.CRS.from_epsg(epsg)

    # reproject and round off
    transformer = pyproj.Transformer.from_crs(
        crs_from=crs, crs_to=crs_utm, always_xy=True
    )
    geom_utm = shapely.ops.transform(transformer.transform, geom)
    # geom can still be a multipolygon for tiles crossing the antimeridian
    geom_merged = shapely.unary_union(geom_utm, grid_size=0.001)
    xx, yy = geom_merged.exterior.coords.xy
    x, y = np.round(np.array(xx)/spac)*spac, np.round(np.array(yy)/spac)*spac

    spac = float(spac)
    nx = int(np.round(np.ptp(x)/spac))
    ny = int(np.round(np.ptp(y)/spac))
    geoTransform = (np.min(x), +spac, 0., np.max(y), 0., -spac, ny, nx)
    return geoTransform, crs_utm

def get_s2_image_locations(fname,s2_df):
    """
    The Sentinel-2 imagery are placed within a folder structure, where one
    folder has an ever changing name. Fortunately this function finds the path
    from the metadata

    Parameters
    ----------
    fname : string
        path string to the Sentinel-2 folder
    s2_df : pandas.dataframe
        index of the bands of interest

    Returns
    -------
    s2_df_new : pandas.dataframe
        dataframe series with relative folder and file locations of the bands
    datastrip_id : string
        folder name of the metadata

    Examples
    --------
    >>> import os
    >>> fpath = '/Users/Data/'
    >>> sname = 'S2A_MSIL1C_20200923T163311_N0209_R140_T15MXV_20200923T200821.SAFE'
    >>> fname = 'MTD_MSIL1C.xml'
    >>> full_path = os.path.join(fpath, sname, fname)
    >>> im_paths,datastrip_id = get_s2_image_locations(full_path)
    >>> im_paths
    ['GRANULE/L1C_T15MXV_A027450_20200923T163313/IMG_DATA/T15MXV_20200923T163311_B01',
     'GRANULE/L1C_T15MXV_A027450_20200923T163313/IMG_DATA/T15MXV_20200923T163311_B02']
    >>> datastrip_id
    'S2A_OPER_MSI_L1C_DS_VGS1_20200923T200821_S20200923T163313_N02.09'
    """
    if os.path.isdir(fname): fname = os.path.join(fname, 'MTD_MSIL1C.xml')
    assert os.path.isfile(fname), ('metafile does not seem to be present')
    root = get_root_of_table(fname)
    root_dir = os.path.split(fname)[0]

    for att in root.iter('Granule'):
        datastrip_full = att.get('datastripIdentifier')
        datastrip_id = '_'.join(datastrip_full.split('_')[4:-1])
    assert datastrip_id != None, ('metafile does not have required info')

    im_paths, band_id = [], []
    for im_loc in root.iter('IMAGE_FILE'):
        full_path = os.path.join(root_dir, im_loc.text)
        boi = im_loc.text[-3:] # band of interest
        if not os.path.isdir(os.path.dirname(full_path)):
            full_path = None
            for _,_,files in os.walk(root_dir):
                for file in files:
                    if file.find(boi)==-1: continue
                    if file.endswith(tuple(['.jp2', '.tif'])):
                        full_path = os.path.join(root_dir,
                                                 file.split('.')[0])
            if full_path is None: # file does not seem to be present
                continue

        im_paths.append(full_path)
        band_id.append(boi)

    band_path = pd.Series(data=im_paths, index=band_id, name="filepath")
    s2_df_new = pd.concat([s2_df, band_path], axis=1, join="inner")
    return s2_df_new, datastrip_id

def get_s2_granule_id(fname, s2_df):
    assert os.path.isfile(fname), ('metafile does not seem to be present')
    root = get_root_of_table(fname)

    for im_loc in root.iter('IMAGE_FILE'):
        granule_id = im_loc.text.split(os.sep)[1]
    return granule_id

def meta_s2string(s2_str):
    """ get meta information of the Sentinel-2 file name

    Parameters
    ----------
    s2_str : string
        filename of the L1C data

    Returns
    -------
    s2_time : string
        date "+YYYY-MM-DD"
    s2_orbit : string
        relative orbit "RXXX"
    s2_tile : string
        tile code "TXXXXX"

    Examples
    --------
    >>> s2_str = 'S2A_MSIL1C_20200923T163311_N0209_R140_T15MXV_20200923T200821.SAFE'
    >>> s2_time, s2_orbit, s2_tile = meta_s2string(s2_str)
    >>> s2_time
    '+2020-09-23'
    >>> s2_orbit
    'R140'
    >>> s2_tile
    'T15MXV'
    """
    assert isinstance(s2_str, str), ("please provide a string")

    if s2_str[0:2]=='S2': # some have info about orbit and sensor
        s2_split = s2_str.split('_')
        s2_time = s2_split[2][0:8]
        s2_orbit = s2_split[4]
        s2_tile = s2_split[5]
    elif s2_str[0:1]=='T': # others have no info about orbit, sensor, etc.
        s2_split = s2_str.split('_')
        s2_time = s2_split[1][0:8]
        s2_tile = s2_split[0]
        s2_orbit = None
    else:
        assert True, "please provide a Sentinel-2 file string"
    # convert to +YYYY-MM-DD string
    # then it can be in the meta-data of following products
    s2_time = '+' + s2_time[0:4] + '-' + s2_time[4:6] + '-' + s2_time[6:8]
    return s2_time, s2_orbit, s2_tile

def get_s2_folders(im_path):
    assert os.path.isdir(im_path), 'please specify a folder'
    s2_list = [x for x in os.listdir(im_path)
               if (os.path.isdir(os.path.join(im_path,x))) & (x[0:2]=='S2')]
    return s2_list

def get_tiles_from_s2_list(s2_list):
    tile_list = [x.split('_')[5] for x in s2_list]
    tile_list = list(set(tile_list))
    return tile_list

def get_utm_from_s2_tiles(tile_list):
    utm_list = [int(x[1:3]) for x in tile_list]
    utm_list = list(set(utm_list))
    return utm_list

def get_utmzone_from_tile_code(tile_code):
    """

    Parameters
    ----------
    tile_code : string, e.g.: '05VMG'
        MGRS tile coding

    Returns
    -------
    utmzone : integer
        code used to denote the number of the projection column of UTM

    See Also
    --------
    .get_epsg_from_mgrs_tile, .get_crs_from_mgrs_tile

    Notes
    -----
    The tile structure is a follows "AABCC"
        * "AA" utm zone number, starting from the East, with steps of 8 degrees
        * "B" latitude zone, starting from the South, with steps of 6 degrees

    The following acronyms are used:

    - CRS : coordinate reference system
    - MGRS : US military grid reference system
    - UTM : universal transverse mercator
    - WGS : world geodetic system
    """
    tile_code = check_mgrs_code(tile_code)
    return int(tile_code[:2])

def get_epsg_from_mgrs_tile(tile_code):
    """

    Parameters
    ----------
    tile_code : string, e.g.: '05VMG'
        MGRS tile coding

    Returns
    -------
    epsg : integer
        code used to denote a certain database entry

    See Also
    --------
    dhdt.generic.get_utmzone_from_mgrs_tile
    dhdt.generic.get_crs_from_mgrs_tile

    Notes
    -----
    The tile structure is a follows "AABCC"
        * "AA" utm zone number, starting from the East, with steps of 8 degrees
        * "B" latitude zone, starting from the South, with steps of 6 degrees

    The following acronyms are used:

    - CRS : coordinate reference system
    - EPSG : european petroleum survey group (a coordinate refenence database)
    - MGRS : US military grid reference system
    - UTM : universal transverse mercator
    - WGS : world geodetic system
    """
    tile_code = check_mgrs_code(tile_code)
    utm_num = get_utmzone_from_tile_code(tile_code)
    epsg_code = 32600 + utm_num

    # N to X are in the Northern hemisphere
    if tile_code[2] < 'N': epsg_code += 100
    return epsg_code

def get_crs_from_mgrs_tile(tile_code):
    """

    Parameters
    ----------
    tile_code : string
        US Military Grid Reference System (MGRS) tile code

    Returns
    -------
    crs : osgeo.osr.SpatialReference
        target projection system

    See Also
    --------
    .get_utmzone_from_mgrs_tile, .get_utmzone_from_mgrs_tile

    Notes
    -----
    The tile structure is a follows "AABCC"
        * "AA" utm zone number, starting from the East, with steps of 8 degrees
        * "B" latitude zone, starting from the South, with steps of 6 degrees
    """
    tile_code = check_mgrs_code(tile_code)
    epsg_code = get_epsg_from_mgrs_tile(tile_code)

    crs = osr.SpatialReference()
    crs.ImportFromEPSG(epsg_code)
    return crs
