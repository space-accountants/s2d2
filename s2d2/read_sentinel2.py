#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# generic libraries
import glob
import os

from xml.etree import ElementTree
from osgeo import osr

import numpy as np
import pandas as pd

# raster/image libraries
from PIL import Image, ImageDraw
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import label
from scipy.signal import convolve2d

from s2d2.handler.xml import get_array_from_xml, get_root_of_table
from s2d2.mapping_input import read_geo_image, read_geo_info
from s2d2.image_coordinate_tools import map2pix, get_bbox
from s2d2.checking.array import make_same_size

def list_central_wavelength_msi():
    """ create dataframe with metadata about Sentinel-2

    Returns
    -------
    pandas.dataframe
        metadata and general multispectral information about the MSI
        instrument that is onboard Sentinel-2, having the following collumns:

            * center_wavelength, unit=µm : central wavelength of the band
            * full_width_half_max, unit=µm : extent of the spectral sensativity
            * bandid : number for identification in the meta data
            * resolution, unit=m : spatial resolution of a pixel
            * along_pixel_size, unit=µm : physical size of the sensor
            * across_pixel_size, unit=µm : size of the photosensative sensor
            * focal_length, unit=m : lens characteristic of the telescope
            * field_of_view, unit=degrees : angle of swath of instrument
            * crossdetector_parallax, unit=degress : in along-track direction
            * name : general name of the band, if applicable
            * solar_illumination, unit=W m-2 µm-1 :
        mostly following the naming convension of the STAC EO extension [stac]_

    Notes
    -----
    The Multi Spectral instrument (MSI) has a Tri Mirror Anastigmat (TMA)
    telescope, which is opened at F/4 and has a focal length of 60 centimeter.

    The crossdetector parallex is given here as well, see [La15]_. While the
    crossband parallax is 0.018 for VNIR and 0.010 for the SWIR detector,
    respectively. The dimensions of the Silicon CMOS detector are given in
    [MG10]_. While the SWIR detector is based on a MCT sensor.

    The following acronyms are used:

        - MSI : multi-spectral instrument
        - MCT : mercury cadmium telluride, HgCdTe
        - TMA : trim mirror anastigmat
        - VNIR : very near infra-red
        - SWIR : short wave infra-red
        - CMOS : silicon complementary metal oxide semiconductor, Si

    Examples
    --------
    >>> from s2d2.read_sentinel2 import list_central_wavelength_msi

    make a selection by name:

    >>> boi = ['red', 'green', 'blue', 'near infrared']
    >>> s2_df = list_central_wavelength_msi()
    >>> s2_df = s2_df[s2_df['common_name'].isin(boi)]
    >>> s2_df
             wavelength  bandwidth  resolution       name  bandid
    B02         492         66          10           blue       1
    B03         560         36          10          green       2
    B04         665         31          10            red       3
    B08         833        106          10  near infrared       7

    similarly you can also select by pixel resolution:

    >>> s2_df = list_central_wavelength_msi()
    >>> tw_df = s2_df[s2_df['resolution']==20]
    >>> tw_df.index
    Index(['B05', 'B06', 'B07', 'B8A', 'B11', 'B12'], dtype='object')

    References
    -------
    .. [La15] Languille et al. "Sentinel-2 geometric image quality commissioning:
              first results" Proceedings of the SPIE, 2015.
    .. [MG10] Martin-Gonthier et al. "CMOS detectors for space applications:
              From R&D to operational program with large volume foundry",
              Proceedings of the SPIE conference on sensors, systems, and next
              generation satellites XIV, 2010.
    .. [stac] https://github.com/stac-extensions/eo
    """
    center_wavelength = {"B01": 443, "B02": 492, "B03": 560, "B04": 665,
                  "B05": 704, "B06": 741, "B07": 783, "B08": 833, "B8A": 865,
                  "B09": 945, "B10":1374, "B11":1614, "B12":2202,
                  }
    # convert from nm to µm
    center_wavelength = {k: v/1E3 for k, v in center_wavelength.items()}

    full_width_half_max = {"B01": 21, "B02": 66, "B03": 36, "B04": 31,
                 "B05": 15, "B06": 15, "B07": 20, "B08":106, "B8A": 21,
                 "B09": 20, "B10": 31, "B11": 91, "B12":175,
                 }
    # convert from nm to µm
    full_width_half_max = {k: v/1E3 for k, v in full_width_half_max.items()}

    bandid = {"B01": 0, "B02": 1, "B03": 2, "B04": 3,
              "B05": 4, "B06": 5, "B07": 6, "B08": 7, "B8A": 8,
              "B09": 9, "B10":10, "B11":11, "B12":12,
                  }
    gsd = {"B01": 60, "B02": 10, "B03": 10, "B04": 10,
           "B05": 20, "B06": 20, "B07": 20, "B08": 10, "B8A": 20,
           "B09": 60, "B10": 60, "B11": 20, "B12": 20,
           }

    along_pixel_size = {"B01": 45., "B02": 7.5, "B03": 7.5, "B04": 7.5,
                        "B05": 15., "B06": 15., "B07": 15., "B08": 7.5,
                        "B8A": 15., "B09": 45., "B10": 15., "B11": 15.,
                        "B12": 15.,}

    across_pixel_size = {"B01": 15., "B02": 7.5, "B03": 7.5, "B04": 7.5,
                        "B05": 15., "B06": 15., "B07": 15., "B08": 7.5,
                        "B8A": 15., "B09": 15., "B10": 15., "B11": 15.,
                        "B12": 15.,}

    focal_length = {"B01": .598, "B02": .598, "B03": .598, "B04": .598,
                    "B05": .598, "B06": .598, "B07": .598, "B08": .598,
                    "B8A": .598, "B09": .598, "B10": .598, "B11": .598,
                    "B12": .598}
    # along_track_view_angle =
    field_of_view = {"B01": 20.6, "B02": 20.6, "B03": 20.6, "B04": 20.6,
                     "B05": 20.6, "B06": 20.6, "B07": 20.6, "B08": 20.6,
                     "B8A": 20.6, "B09": 20.6, "B10": 20.6, "B11": 20.6,
                     "B12": 20.6}
    # given as base-to-height
    crossdetector_parallax = {"B01": 0.055, "B02": 0.022, "B03": 0.030,
                              "B04": 0.034, "B05": 0.038, "B06": 0.042,
                              "B07": 0.046, "B08": 0.026, "B8A": 0.051,
                              "B09": 0.059, "B10": 0.030, "B11": 0.040,
                              "B12": 0.050,
                          }
    # convert to angles in degrees
    crossdetector_parallax = \
        {k: 2*np.tan(v/2) for k,v in crossdetector_parallax.items()}

    common_name = {"B01": 'coastal',       "B02" : 'blue',
                   "B03" : 'green',        "B04" : 'red',
                   "B05" : 'rededge',      "B06" : 'rededge',
                   "B07" : 'rededge',      "B08" : 'nir',
                   "B8A" : 'nir08',        "B09" : 'nir09',
                   "B10": 'cirrus',        "B11" : 'swir16',
                   "B12": 'swir22'}
    solar_illumination = {"B01":1913, "B02":1942, "B03":1823, "B04":1513,
                  "B05":1426, "B06":1288, "B07":1163, "B08":1036, "B8A":955,
                  "B09": 813, "B10": 367, "B11": 246, "B12": 85,
                  } # these numbers are also given in the meta-data
    d = {
         "center_wavelength": pd.Series(center_wavelength,
                                        dtype=np.dtype('float')),
         "full_width_half_max": pd.Series(full_width_half_max,
                                          dtype=np.dtype('float')),
         "gsd": pd.Series(gsd,
                          dtype=np.dtype('float')),
         "across_pixel_size": pd.Series(across_pixel_size,
                                        dtype=np.dtype('float')),
         "along_pixel_size": pd.Series(along_pixel_size,
                                       dtype=np.dtype('float')),
         "focal_length": pd.Series(focal_length,
                                   dtype=np.dtype('float')),
         "common_name": pd.Series(common_name,
                                  dtype=np.dtype('str')),
         "bandid": pd.Series(bandid,
                             dtype=np.dtype('int64')),
         "field_of_view" : pd.Series(field_of_view,
                                     dtype=np.dtype('float')),
         "solar_illumination": pd.Series(solar_illumination,
                                         dtype=np.dtype('float')),
         "crossdetector_parallax": pd.Series(crossdetector_parallax,
                                             dtype=np.dtype('float'))
    }
    df = pd.DataFrame(d)
    return df


def read_stack_s2(s2_df):
    """ read imagery data of interest into an three dimensional np.array

    Parameters
    ----------
    s2_df : pandas.Dataframe
        metadata and general multispectral information about the MSI
        instrument that is onboard Sentinel-2

    Returns
    -------
    im_stack : numpy.ndarray, ndim=3
        array of the band image
    spatialRef : string
        projection
    geoTransform : tuple, size=(8,1)
        affine transformation coefficients
    targetprj : osr.SpatialReference object
        spatial reference

    See Also
    --------
    list_central_wavelength_msi :
        creates a dataframe for the MSI instrument
    get_s2_image_locations :
        provides dataframe with specific file locations
    read_band_s2 :
        reading a single Sentinel-2 band

    Examples
    --------
    >>> s2_df = list_central_wavelength_msi()
    >>> s2_df = s2_df[s2_df['gsd'] == 10]
    >>> s2_df,_ = get_s2_image_locations(IM_PATH, s2_df)

    >>> im_stack, spatialRef, geoTransform, targetprj = read_stack_s2(s2_df)
    """
    assert isinstance(s2_df, pd.DataFrame), ('please provide a dataframe')
    assert 'filepath' in s2_df, ('please first run "get_s2_image_locations"'+
                                ' to find the proper file locations')

    roi = np.min(s2_df['gsd'].array) # resolution of interest
    s2_df = s2_df['gsd']==roi
    # only import the highest resolution
    for val, idx in enumerate(s2_df.sort_values('gsd').index):
        full_path = s2_df['filepath'][idx] + '.jp2'
        if val==0:
            im_stack, spatialRef, geoTransform, targetprj = read_band_s2(full_path)
        else: # stack others bands
            if im_stack.ndim==2:
                im_stack = np.atleast_3d(im_stack)
            band = np.atleast_3d(read_band_s2(full_path)[0])
            im_stack = np.concatenate((im_stack, band), axis=2)

    # in the meta data files there is no mention of its size, hence include
    # it to the geoTransform
    if len(geoTransform)==6:
         geoTransform = geoTransform + (im_stack.shape[0], im_stack.shape[1])
    return im_stack, spatialRef, geoTransform, targetprj

def _get_geometric_info_from_root(root):
    geom_info = None
    for child in root:
        if child.tag.endswith('Geometric_Info'):
            geom_info = child
    assert(geom_info is not None), ('metadata not in xml file')
    return geom_info

def _get_frame_from_geometric_info(geom_info):
    frame = None
    for tile in geom_info:
        if tile.tag == 'Tile_Geocoding':
            frame = tile
    assert (frame is not None), ('metadata not in xml file')
    return frame

def _get_image_size_s2_from_root(root, resolution=10):
    geom_info = _get_geometric_info_from_root(root)
    frame = _get_frame_from_geometric_info(geom_info)

    m,n = None,None
    for box in frame:
        if (box.tag=='Size') and \
            (box.attrib['resolution']==str(resolution)):
            for field in box:
                if field.tag=='NROWS':
                    m = int(field.text)
                elif field.tag=='NCOLS':
                    n = int(field.text)
    return m, n

def _get_ul_coord_s2_from_root(root, resolution=10):
    geom_info = _get_geometric_info_from_root(root)
    frame = _get_frame_from_geometric_info(geom_info)

    for box in frame:
        if (box.tag=='Geoposition') and \
            (box.attrib['resolution']==str(resolution)):
            for field in box:
                if field.tag=='ULX':
                    ul_x = float(field.text)
                elif field.tag=='ULY':
                    ul_y = float(field.text)
    return ul_x, ul_y

def _get_geotransform_s2_from_root(root, resolution=10):
    geom_info = _get_geometric_info_from_root(root)
    frame = _get_frame_from_geometric_info(geom_info)

    for box in frame:
        if (box.tag=='Geoposition') and \
            (box.attrib['resolution']==str(resolution)):
            for field in box:
                if field.tag == 'ULX':
                    ul_X = float(field.text)
                elif field.tag == 'ULY':
                    ul_Y = float(field.text)
                elif field.tag == 'XDIM':
                    d_X = float(field.text)
                elif field.tag == 'YDIM':
                    d_Y = float(field.text)
        elif (box.tag=='Size') and \
            (box.attrib['resolution']==str(resolution)):
            for field in box:
                if field.tag == 'NROWS':
                    m = int(field.text)
                elif field.tag == 'NCOLS':
                    n = int(field.text)
    geoTransform = (ul_X, d_X, 0., ul_Y, 0., d_Y, m, n)
    return geoTransform

def _get_crs_s2_from_root(root):
    geom_info = _get_geometric_info_from_root(root)
    frame = _get_frame_from_geometric_info(geom_info)

    epsg = None
    for field in frame:
        if field.tag == 'HORIZONTAL_CS_CODE':
            epsg = field.text
            epsg_int = int(epsg.split(':')[1])

    crs = osr.SpatialReference()
    crs.ImportFromEPSG(epsg_int)
    return crs

def _get_orbit_s2_from_root(root):
    gnrl_info = None
    for child in root:
        if child.tag.endswith('General_Info'):
            gnrl_info = child
    assert(gnrl_info is not None), ('metadata not in xml file')

    prod_info = None
    for child in gnrl_info:
        if child.tag == 'Product_Info':
            prod_info = child
    assert(prod_info is not None), ('metadata not in xml file')

    data_take = None
    for child in prod_info:
        if child.tag == 'Datatake':
            data_take = child
    assert(data_take is not None), ('metadata not in xml file')

    row = None
    for field in data_take:
        if field.tag == 'SENSING_ORBIT_NUMBER':
            row = int(field.text)
    return row

def _get_spacecraft_s2_from_root(root):
    gnrl_info = None
    for child in root:
        if child.tag.endswith('General_Info'):
            gnrl_info = child
    assert(gnrl_info is not None), ('metadata not in xml file')

    prod_info = None
    for child in gnrl_info:
        if child.tag == 'Product_Info':
            prod_info = child
    assert(prod_info is not None), ('metadata not in xml file')

    data_take = None
    for child in prod_info:
        if child.tag == 'Datatake':
            data_take = child
    assert(data_take is not None), ('metadata not in xml file')

    version = None
    for field in data_take:
        if field.tag == 'SPACECRAFT_NAME':
            version = field.text[-1].upper()
    return version

def read_geotransform_s2(path, fname='MTD_TL.xml', resolution=10):
    """ get the mapping transformation of the Sentinel-2 image

    Parameters
    ----------
    path : string
        location where the meta data is situated
    fname : string
        file name of the meta-data file
    resolution : {float,integer}, unit=meters, default=10
        resolution of the grid

    Returns
    -------
    geoTransform : tuple, size=(1,6)
        affine transformation coefficients

    See Also
    --------
    s2d2.image_coordinate_tools.map2pix : transform map- to image-coordinates
    s2d2.image_coordinate_tools.map2pix : transform image- to map-coordinates

    Notes
    -----
    The metadata is scattered over the file structure of Sentinel-2, L1C

    .. code-block:: text

        * S2X_MSIL1C_20XX...
        ├ AUX_DATA
        ├ DATASTRIP
        │  └ DS_XXX_XXXX...
        │     └ QI_DATA
        │        └ MTD_DS.xml <- metadata about the data-strip
        ├ GRANULE
        │  └ L1C_TXXXX_XXXX...
        │     ├ AUX_DATA
        │     ├ IMG_DATA
        │     ├ QI_DATA
        │     └ MTD_TL.xml <- metadata about the tile
        ├ HTML
        ├ rep_info
        ├ manifest.safe
        ├ INSPIRE.xml
        └ MTD_MSIL1C.xml <- metadata about the product

    The metadata structure is as follows:
        .. code-block:: text

        * MTD_TL.xml
        └ n1:Level-1C_Tile_ID
           ├ n1:General_Info
           ├ n1:Geometric_Info
           │  ├ Tile_Geocoding
           │  │  ├ HORIZONTAL_CS_NAME
           │  │  ├ HORIZONTAL_CS_CODE
           │  │  ├ Size : resolution={"10","20","60"}
           │  │  │  ├ NROWS
           │  │  │  └ NCOLS
           │  │  └ Geoposition
           │  │     ├ ULX
           │  │     ├ ULY
           │  │     ├ XDIM
           │  │     └ YDIM
           │  └ Tile_Angles
           └ n1:Quality_Indicators_Info

    The following acronyms are used:

    - DS : datastrip
    - TL : tile
    - QI : quality information
    - AUX : auxiliary
    - MTD : metadata
    - MSI : multi spectral instrument
    - L1C : product specification,i.e.: level 1, processing step C
    """
    assert (resolution in (10, 20, 60)), ('please provide correct resolution')
    if not os.path.isfile(os.path.join(path, fname)):
        return None
    root = get_root_of_table(path, fname)

    geoTransform = _get_geotransform_s2_from_root(root, resolution=resolution)
    return geoTransform

def read_orbit_number_s2(path, fname='MTD_MSIL1C.xml'):
    """ read the orbit number of the Sentinel-2 image from the metadata

    """
    root = get_root_of_table(path, fname)
    row = _get_orbit_s2_from_root(root)
    return row

def get_local_bbox_in_s2_tile(fname_1, s2dir):
    _,geoTransform, _, rows, cols, _ = read_geo_info(fname_1)
    bbox_xy = get_bbox(geoTransform, rows, cols)

    s2Transform = read_geotransform_s2(s2dir)
    if s2Transform is None:
        return None

    bbox_i, bbox_j = map2pix(s2Transform, bbox_xy[0:2], bbox_xy[2::])
    bbox_ij = np.concatenate((np.flip(bbox_i), bbox_j)).astype(int)
    return bbox_ij

def get_sun_angles_s2_from_root(root, angle='Zenith'):
    assert angle in ('Zenith','Azimuth',), \
        ('please provide correct angle name')
    geom_info = None
    for child in root:
        if child.tag.endswith('Geometric_Info'):
            geom_info = child
    assert(geom_info is not None), ('metadata not in xml file')

    grids = None
    for tile in geom_info:
        if tile.tag == 'Tile_Angles':
            frame = tile
            for instance in frame:
                if instance.tag=='Sun_Angles_Grid':
                    grids = instance
    assert(grids is not None), ('metadata not in xml file')

    Ang, col_step, row_step = None, None, None
    for box in grids:
        if (box.tag==angle):
            for field in box:
                if field.tag == 'Values_List':
                    Ang = get_array_from_xml(field)
                elif field.tag == 'COL_STEP':
                    col_step = int(field.text)
                elif field.tag == 'ROW_STEP':
                    row_step = int(field.text)

    assert(Ang is not None), ('metadata does not seem to be present')
    return Ang, col_step, row_step

def read_sun_angles_s2(path, fname='MTD_TL.xml'):
    """ This function reads the xml-file of the Sentinel-2 scene and extracts
    an array with sun angles, as these vary along the scene.

    Parameters
    ----------
    path : string
        path where xml-file of Sentinel-2 is situated

    Returns
    -------
    Zn : numpy.array, size=(m,n), dtype=float
        array of the solar zenith angles, in degrees.
    Az : numpy.array, size=(m,n), dtype=float
        array of the solar azimuth angles, in degrees.

    Notes
    -----
    The computation of the solar angles is done in two steps: 1) Computation of
    the solar angles in J2000; 2) Transformation of the vector to the mapping
    frame.

    - 1) The outputs of the first step is the solar direction normalized vector
         with the Earth-Sun distance, considering that the direction of the sun
         is the same at the centre of the Earth and at the centre of the
         Sentinel-2 satellite.

    - 2) Attitude of the satellite platform are used to rotate the solar vector
         to the mapping frame. Also Ground Image Calibration Parameters (GICP
         Diffuser Model) are used to transform from the satellite to the
         diffuser, as Sentinel-2 has a forward and backward looking sensor
         configuration. The diffuser model also has stray-light correction and
         a Bi-Directional Reflection Function model.

    The angle(s) are declared in the following coordinate frame:

        .. code-block:: text

                 ^ North & y
                 |
            - <--|--> +
                 |
                 +----> East & x

    The angles related to the sun are as follows:

        .. code-block:: text

          surface normal              * sun
          ^                     ^    /
          |                     |   /
          |-- zenith angle      |  /
          | /                   | /|
          |/                    |/ | elevation angle
          +----                 +------

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

    The metadata is scattered over the file structure of Sentinel-2, L1C

        .. code-block:: text

            * S2X_MSIL1C_20XX...
            ├ AUX_DATA
            ├ DATASTRIP
            │  └ DS_XXX_XXXX...
            │     └ QI_DATA
            │        └ MTD_DS.xml <- metadata about the data-strip
            ├ GRANULE
            │  └ L1C_TXXXX_XXXX...
            │     ├ AUX_DATA
            │     ├ IMG_DATA
            │     ├ QI_DATA
            │     └ MTD_TL.xml <- metadata about the tile
            ├ HTML
            ├ rep_info
            ├ manifest.safe
            ├ INSPIRE.xml
            └ MTD_MSIL1C.xml <- metadata about the product

    The metadata structure looks like:

        .. code-block:: text

            * MTD_TL.xml
            └ n1:Level-1C_Tile_ID
               ├ n1:General_Info
               ├ n1:Geometric_Info
               │  ├ Tile_Geocoding
               │  └ Tile_Angles
               │     ├ Sun_Angles_Grid
               │     │  ├ Zenith
               │     │  │  ├ COL_STEP
               │     │  │  ├ ROW_STEP
               │     │  │  └ Values_List
               │     │  └ Azimuth
               │     │     ├ COL_STEP
               │     │     ├ ROW_STEP
               │     │     └ Values_List
               │     ├ Mean_Sun_Angle
               │     └ Viewing_Incidence_Angles_Grids
               └ n1:Quality_Indicators_Info

    The following acronyms are used:

    - s2 : Sentinel-2
    - DS : datastrip
    - TL : tile
    - QI : quality information
    - AUX : auxiliary
    - MTD : metadata
    - MSI : multi spectral instrument
    - L1C : product specification,i.e.: level 1, processing step C
    """
    if not os.path.isfile(os.path.join(path, fname)):
        return None, None
    root = get_root_of_table(path, fname)

    mI,nI = _get_image_size_s2_from_root(root)
    Zn = get_sun_angles_s2_from_root(root, angle='Zenith')[0]

    zi = np.linspace(0 - 10, mI + 10, np.size(Zn, axis=0))
    zj = np.linspace(0 - 10, nI + 10, np.size(Zn, axis=1))

    interp = RegularGridInterpolator((zi, zj), Zn)

    iGrd, jGrd = np.arange(0, mI), np.arange(0, nI)
    Igrd, Jgrd = np.meshgrid(iGrd, jGrd)

    Zn = interp((Igrd, Jgrd))
    del zi, zj

    # get Azimuth array
    Az = get_sun_angles_s2_from_root(root, angle='Azimuth')[0]

    ai = np.linspace(0 - 10, mI + 10, np.size(Az, axis=0))
    aj = np.linspace(0 - 10, nI + 10, np.size(Az, axis=1))

    interp = RegularGridInterpolator((ai, aj), Az)
    Az = interp((Igrd, Jgrd))
    del ai, aj
    return Zn, Az

def get_view_angles_s2_from_root(root):
    geom_info = _get_geometric_info_from_root(root)
    frame = _get_frame_from_geometric_info(geom_info)

    Az_grd,Zn_grd = None, None
    bnd, det = None, None
    col_step, row_step = None, None
    for instance in frame:
        if instance.tag == 'Viewing_Incidence_Angles_Grids':
            boi, doi = np.atleast_1d(int(instance.attrib['bandId'])), \
                       np.atleast_1d(int(instance.attrib['detectorId']))
            for box in instance:
                if box.tag in ('Zenith', 'Azimuth'):
                    for field in box:
                        if field.tag == 'Values_List':
                            Grd = np.atleast_3d(get_array_from_xml(field))
                    if box.tag == 'Zenith':
                        if Zn_grd is None:
                            Zn_grd = Grd
                        else:
                            Zn_grd = np.dstack((Zn_grd, Grd))
                    else:
                        if Az_grd is None:
                            Az_grd = Grd
                        else:
                            Az_grd = np.dstack((Az_grd, Grd))
            bnd = boi if bnd is None else np.hstack((bnd, boi))
            det = doi if det is None else np.hstack((det, doi))

    # create geotransform
    for field in box:
        if field.tag == 'COL_STEP':
            dx = float(field.text)
        elif field.tag == 'ROW_STEP':
            dy = float(field.text)
    ul_x,ul_y = _get_ul_coord_s2_from_root(root)

    geoTransform = (ul_x, dx, 0., ul_y, 0., -dy,
                    Zn_grd.shape[0], Zn_grd.shape[1])

    return Az_grd, Zn_grd, bnd, det, geoTransform

def read_view_angles_s2(path, fname='MTD_TL.xml', det_stack=np.array([]),
                        boi_df=None):
    """ This function reads the xml-file of the Sentinel-2 scene and extracts
    an array with viewing angles of the MSI instrument.

    Parameters
    ----------
    path : string
        path where xml-file of Sentinel-2 is situated
    fname : string
        the name of the metadata file, sometimes this is changed
    boi_df : pandas.dataframe, default=4
        each band has a somewhat minute but different view angle

    Returns
    -------
    Zn : numpy.ma.array, size=(m,n), dtype=float
        masked array of the solar zenith angles, in degrees.
    Az : numpy.ma.array, size=(m,n), dtype=float
        masked array of the solar azimuth angles, in degrees.

    See Also
    --------
    list_central_wavelength_s2

    Notes
    -----
    The azimuth angle is declared in the following coordinate frame:

        .. code-block:: text

                 ^ North & y
                 |
            - <--|--> +
                 |
                 └----> East & x

    The angles related to the satellite are as follows:

        .. code-block:: text

               #*#                   #*# satellite
          ^    /                ^    /|
          |   /                 |   / | nadir
          |-- zenith angle      |  /  v
          | /                   | /|
          |/                    |/ | elevation angle
          └----- surface        └------

    The metadata is scattered over the file structure of Sentinel-2, L1C

        .. code-block:: text

            * S2X_MSIL1C_20XX...
            ├ AUX_DATA
            ├ DATASTRIP
            │  └ DS_XXX_XXXX...
            │     └ QI_DATA
            │        └ MTD_DS.xml <- metadata about the data-strip
            ├ GRANULE
            │  └ L1C_TXXXX_XXXX...
            │     ├ AUX_DATA
            │     ├ IMG_DATA
            │     ├ QI_DATA
            │     └ MTD_TL.xml <- metadata about the tile
            ├ HTML
            ├ rep_info
            ├ manifest.safe
            ├ INSPIRE.xml
            └ MTD_MSIL1C.xml <- metadata about the product

    The metadata structure looks like:

        .. code-block:: text

            * MTD_TL.xml
            └ n1:Level-1C_Tile_ID
               ├ n1:General_Info
               ├ n1:Geometric_Info
               │  ├ Tile_Geocoding
               │  └ Tile_Angles
               │     ├ Sun_Angles_Grid
               │     ├ Mean_Sun_Angle
               │     └ Viewing_Incidence_Angles_Grids
               │        ├ Zenith
               │        │  ├ COL_STEP
               │        │  ├ ROW_STEP
               │        │  └ Values_List
               │        └ Azimuth
               │           ├ COL_STEP
               │           ├ ROW_STEP
               │           └ Values_List
               └ n1:Quality_Indicators_Info

    The following acronyms are used:

    - s2 : Sentinel-2
    - DS : datastrip
    - TL : tile
    - QI : quality information
    - AUX : auxiliary
    - MTD : metadata
    - MSI : multi spectral instrument
    - L1C : product specification,i.e.: level 1, processing step C
    """
    if boi_df is None:
        boi_df = list_central_wavelength_msi()
    assert boi_df['gsd'].nunique() == 1, \
        'make sure all bands are the same resolution'
    root = get_root_of_table(path, fname)

    if det_stack.size == 0: # if no stack is given, just import all metadata
        roi = boi_df['gsd'].mode()[0]  # resolution of interest
        # image dimensions
        for meta in root.iter('Size'):
            res = float(meta.get('resolution'))
            if res == roi:
                mI, nI = int(meta[0].text), int(meta[1].text)

        det_stack = read_detector_mask(os.path.join(path,'QI_DATA'),
                                       boi_df,
                                       np.array([0, 0, 0, 0, 0, 0, mI, nI]))

    det_list = np.unique(det_stack)
    bnd_list = np.asarray(boi_df['bandid'])
    Zn_grd, Az_grd = None, None
    # get coarse grids
    for grd in root.iter('Viewing_Incidence_Angles_Grids'):
        bid = float(grd.get('bandId'))
        if not boi_df['bandid'].isin([bid]).any():
            continue
        grd_idx_3 = np.where(bnd_list==bid)[0][0]

        # link to band in detector stack, by finding the integer index
        det_idx = np.flatnonzero(boi_df['bandid'].isin([bid]))[0]

        # link to detector id, within this band
        did = float(grd.get('detectorId'))
        if not np.any(det_list==did):
            continue
        grd_idx_4 = np.where(det_list==did)[0][0]

        col_sep = float(grd[0][0].text) # in meters
        row_sep = float(grd[0][0].text) # in meters
        col_idx = boi_df.columns.get_loc('gsd')
        col_res = float(boi_df.iloc[boi_df.index[det_idx],col_idx])
        row_res = float(boi_df.iloc[boi_df.index[det_idx],col_idx])

        Zarray = get_array_from_xml(grd[0][2])
        if Zn_grd is None:
            Zn_grd = np.zeros((Zarray.shape[0], Zarray.shape[1],
                              len(bnd_list), len(det_list)), dtype=np.float64)
        Zn_grd[:,:,grd_idx_3,grd_idx_4] = Zarray
        Aarray = get_array_from_xml(grd[1][2])
        if Az_grd is None:
            Az_grd = np.zeros((Aarray.shape[0], Aarray.shape[1],
                              len(bnd_list), len(det_list)), dtype=np.float64)
        Az_grd[:,:,grd_idx_3,grd_idx_4] = Aarray

    mI, nI = det_stack.shape[0], det_stack.shape[1]

    # create grid coordinate frame
    I_grd, J_grd = np.mgrid[0:Zn_grd.shape[0], 0:Zn_grd.shape[1]]
    I_grd, J_grd = I_grd.astype('float64'), J_grd.astype('float64')
    I_grd *= (col_sep / col_res)
    J_grd *= (row_sep / row_res)

    return Zn_grd, Az_grd, I_grd, J_grd

def read_mean_sun_angles_s2(path, fname='MTD_TL.xml'):
    """ Read the xml-file of the Sentinel-2 scene and extract the mean sun angles.

    Parameters
    ----------
    path : string
        path where xml-file of Sentinel-2 is situated

    Returns
    -------
    Zn : float
        Mean solar zentih angle of the scene, in degrees.
    Az : float
        Mean solar azimuth angle of the scene, in degrees

    Notes
    -----
    The azimuth angle declared in the following coordinate frame:

        .. code-block:: text

                 ^ North & y
                 |
            - <--|--> +
                 |
                 └----> East & x

    The angles related to the sun are as follows:

        .. code-block:: text

          surface normal              * sun
          ^                     ^    /
          |                     |   /
          |-- zenith angle      |  /
          | /                   | /|
          |/                    |/ | elevation angle
          └----                 └------

    The metadata is scattered over the file structure of Sentinel-2, L1C

        .. code-block:: text

            * S2X_MSIL1C_20XX...
            ├ AUX_DATA
            ├ DATASTRIP
            │  └ DS_XXX_XXXX...
            │     └ QI_DATA
            │        └ MTD_DS.xml <- metadata about the data-strip
            ├ GRANULE
            │  └ L1C_TXXXX_XXXX...
            │     ├ AUX_DATA
            │     ├ IMG_DATA
            │     ├ QI_DATA
            │     └ MTD_TL.xml <- metadata about the tile
            ├ HTML
            ├ rep_info
            ├ manifest.safe
            ├ INSPIRE.xml
            └ MTD_MSIL1C.xml <- metadata about the product

    The following acronyms are used:

    - s2 : Sentinel-2
    - DS : datastrip
    - TL : tile
    - QI : quality information
    - AUX : auxiliary
    - MTD : metadata
    - MSI : multi spectral instrument
    - L1C : product specification,i.e.: level 1, processing step C
    """
    root = get_root_of_table(path, fname)
    Zn, Az = [], []
    for att in root.iter('Mean_Sun_Angle'):
        Zn = float(att[0].text)
        Az = float(att[1].text)
    return Zn, Az

def read_detector_mask(path_meta, boi, geoTransform):
    """ create array of with detector identification

    Sentinel-2 records in a pushbroom fasion, collecting reflectance with a
    ground resolution of more than 270 kilometers. This data is stacked in the
    flight direction, and then cut into granules. However, the sensorstrip
    inside the Multi Spectral Imager (MSI) is composed of several CCD arrays.

    This function collects the geometry of these sensor arrays from the meta-
    data. Since this is stored in a gml-file.

    Parameters
    ----------
    path_meta : string
        path where the meta-data is situated.
    boi : pandas.DataFrame
        list with bands of interest
    geoTransform : tuple, size={(6,),(8,)}
        affine transformation coefficients

    Returns
    -------
    numpy.ma.array, size=(msk_dim[0:1],len(boi)), dtype=int8
        array where each pixel has the ID of the detector, of a specific band

    See Also
    --------
    list_central_wavelength_msi : creates a dataframe for the MSI instrument
    read_view_angles_s2 : read grid of Sentinel-2 observation angles

    Notes
    -----
    The metadata is scattered over the file structure of Sentinel-2, L1C

    .. code-block:: text

        * S2X_MSIL1C_20XX...
        ├ AUX_DATA
        ├ DATASTRIP
        │  └ DS_XXX_XXXX...
        │     └ QI_DATA
        │        └ MTD_DS.xml <- metadata about the data-strip
        ├ GRANULE
        │  └ L1C_TXXXX_XXXX...
        │     ├ AUX_DATA
        │     ├ IMG_DATA
        │     ├ QI_DATA
        │     └ MTD_TL.xml <- metadata about the tile
        ├ HTML
        ├ rep_info
        ├ manifest.safe
        ├ INSPIRE.xml
        └ MTD_MSIL1C.xml <- metadata about the product

    The following acronyms are used:

    - DS : datastrip
    - TL : tile
    - QI : quality information
    - AUX : auxiliary
    - MTD : metadata
    - MSI : multi spectral instrument
    - L1C : product specification,i.e.: level 1, processing step C

    Examples
    --------
    >>> from s2d2.read_sentinel2 import (
    >>>    list_central_wavelength_msi, read_detector_mask
    >>> )

    >>> path_meta = '/GRANULE/L1C_T15MXV_A027450_20200923T163313/QI_DATA'
    >>> boi = ['red', 'green', 'blue', 'near infrared']
    >>> s2_df = list_central_wavelength_msi()
    >>> boi_df = s2_df[s2_df['common_name'].isin(boi)]
    >>> geoTransform = (600000.0, 10.0, 0.0, 10000000.0, 0.0, -10.0)
    >>>
    >>> det_stack = read_detector_mask(path_meta, boi_df, geoTransform)
    """
    assert isinstance(path_meta, str), ('please provide a string')
    assert isinstance(boi, pd.DataFrame), ('please provide a dataframe')
    assert isinstance(geoTransform, tuple)

    det_stack = None
    if len(geoTransform)>6: # also image size is given
        msk_dim = (int(geoTransform[-2]), int(geoTransform[-1]), len(boi))
        det_stack = np.zeros(msk_dim, dtype='int8')  # create stack

    for i in range(len(boi)):
        im_id = boi.index[i] # 'B01' | 'B8A'
        if type(im_id) is int:
            f_meta = os.path.join(path_meta, 'MSK_DETFOO_B'+ \
                                  f'{im_id:02.0f}' + '.gml')
        else:
            f_meta = os.path.join(path_meta, 'MSK_DETFOO_'+ im_id +'.gml')

        if os.path.isfile(f_meta): # old Sentinel-2 files had gml-files
            dom = ElementTree.parse(glob.glob(f_meta)[0])
            root = dom.getroot()

            if det_stack is None:
                msk_dim = get_msk_dim_from_gml(root, geoTransform)
                msk_dim = msk_dim + (len(boi),)
                det_stack = np.zeros(msk_dim, dtype='int8') # create stack

            mask_members = root[2]
            for k in range(len(mask_members)):
                pos_arr, det_num = get_xy_poly_from_gml(mask_members, k)

                # transform to image coordinates
                i_arr, j_arr = map2pix(geoTransform, pos_arr[:,0], pos_arr[:,1])
                ij_arr = np.hstack((j_arr[:,np.newaxis], i_arr[:,np.newaxis]))
                # make mask
                msk = Image.new("L", [np.size(det_stack,1), np.size(det_stack,0)], 0)
                ImageDraw.Draw(msk).polygon(tuple(map(tuple, ij_arr[:,0:2])),
                                            outline=det_num, fill=det_num)
                msk = np.array(msk)
                det_stack[:,:,i] = np.maximum(det_stack[:,:,i], msk)
        else:
            im_meta = f_meta[:-3]+'jp2'
            assert os.path.isfile(im_meta), ('gml and jp2 file not present')
            msk,_,mskTransform,_ = read_geo_image(im_meta)
            if det_stack is None:
                det_stack = np.zeros((*msk.shape[:2], len(boi)), dtype='int8')
            if mskTransform!=geoTransform:
                msk = make_same_size(msk, mskTransform, geoTransform)
            det_stack[..., i] = msk

    det_stack = np.ma.array(det_stack, mask=det_stack==0)
    return det_stack

def read_sensing_time_s2(path, fname='MTD_TL.xml'):
    """

    Parameters
    ----------
    path : string
        path where the meta-data is situated
    fname : string
        file name of the metadata.

    Returns
    -------
    numpy.datetime64
        time of image acquisition

    See Also
    --------
    get_s2_granual_id

    Notes
    -----
    The recording time is the average sensing within the tile, which is a
    combination of forward and backward looking datastrips. Schematically, this
    might look like:

        .. code-block:: text

             x recording time of datastrip, i.e.: the start of the strip
             × recording time of tile, i.e.: the weighted average
                                        _________  datastrip
                              ____x____/        /_________
                             /        /        /        /
                            /        /        /        /
                           /        /        /        /
                          /        /        /        /
                         /        /        /        /
                        /        /        /        /
                       /    ┌---/------┐ /        /
                      /     |  /       |/        /
                     /      | /  ×     /        /
                    /       |/        /|       /
                   /        /--------/-┘      /
                  /        /  tile  /        /
                 /        /        /        /
                /        /        /        /
               /        /        /        /
              /        /        /        /
             /        _________/        /
             _________        /________/

    Examples
    --------
    demonstrate the code when in a Scihub data structure

    >>> import os
    >>> from s2d2.handler.sentinel2 import get_s2_image_locations
    >>> from s2d2.read_sentinel2 import read_sensing_time_s2

    >>> fpath = '/Users/Data/'
    >>> sname = 'S2A_MSIL1C_20200923T163311_N0209_R140_T15MXV_20200923T200821.SAFE'
    >>> fname = 'MTD_MSIL1C.xml'
    >>> s2_path = os.path.join(fpath, sname, fname)
    >>> s2_df,_ = get_s2_image_locations(fname, s2_df)
    >>> s2_df['filepath']
    B02 'GRANULE/L1C_T15MXV_A027450_20200923T163313/IMG_DATA/T115MXV...'
    B03 'GRANULE/L1C_T15MXV_A027450_20200923T163313/IMG_DATA/T15MXV...'
    >>> full_path = '/'.join(s2_df['filepath'][0].split('/')[:-2])
    'GRANULE/L1C_T15MXV_A027450_20200923T163313'
    >>> rec_time = read_sensing_time_s2(path, fname='MTD_TL.xml')
    >>> rec_time

    """
    assert os.path.isfile(os.path.join(path,fname)), \
        ('file does not seem to exist')
    root = get_root_of_table(path, fname)
    for att in root.iter('Sensing_Time'.upper()):
        time_str = att.text
        if time_str.endswith('Z'): time_str = time_str[:-1]
        rec_time = np.datetime64(time_str, 'ns')
    return rec_time

#todo: robustify
def get_xy_poly_from_gml(gml_struct,idx):
    # get detector number from meta-data
    det_id = gml_struct[idx].attrib
    det_id = list(det_id.items())[0][1].split('-')[2]
    det_num = int(det_id)

    # get footprint
    pos_dim = gml_struct[idx][1][0][0][0][0].attrib
    pos_dim = int(list(pos_dim.items())[0][1])
    pos_list = gml_struct[idx][1][0][0][0][0].text
    pos_row = [float(s) for s in pos_list.split(' ')]
    pos_arr = np.array(pos_row).reshape((int(len(pos_row) / pos_dim), pos_dim))
    return pos_arr, det_num

#todo: robustify
def get_msk_dim_from_gml(gml_struct, geoTransform):
    assert isinstance(geoTransform, tuple)
    # find dimensions of array through its map extent in metadata
    lower_corner = np.fromstring(gml_struct[1][0][0].text, sep=' ')
    upper_corner = np.fromstring(gml_struct[1][0][1].text, sep=' ')
    lower_x, lower_y = map2pix(geoTransform, lower_corner[0], lower_corner[1])
    upper_x, upper_y = map2pix(geoTransform, upper_corner[0], upper_corner[1])

    msk_dim = (np.maximum(lower_x, upper_x).astype(int),
               np.maximum(lower_y, upper_y).astype(int))
    return msk_dim

def read_cloud_mask(path_meta, geoTransform):
    """

    Parameters
    ----------
    path_meta : string
        directory where meta data is situated, 'MSK_CLOUDS_B00.gml' is typically
        the file of interest
    geoTransform : tuple
        affine transformation coefficients

    Returns
    -------
    numpy.ndarray, size=(m,n), dtype=int
        array which highlights where clouds could be
    """
    assert isinstance(geoTransform, tuple)

    f_meta = os.path.join(path_meta, 'MSK_CLOUDS_B00.gml')
    root = get_root_of_table(f_meta)

    if len(geoTransform)>6: # also image size is given
        msk_dim = (geoTransform[-2], geoTransform[-1])
        msk_clouds = np.zeros(msk_dim, dtype='int8')
    else:
        msk_dim = get_msk_dim_from_gml(root)
        msk_clouds = np.zeros(msk_dim, dtype='int8')  # create stack

    if len(root)>2: # look into meta-data for cloud polygons
        mask_members = root[2]
        for k in range(len(mask_members)):
            pos_arr = get_xy_poly_from_gml(mask_members, k)[0]

            # transform to image coordinates
            i_arr, j_arr = map2pix(geoTransform, pos_arr[:, 0], pos_arr[:, 1])
            ij_arr = np.hstack((j_arr[:, np.newaxis], i_arr[:, np.newaxis]))

            # make mask
            msk = Image.new("L", [msk_dim[1], msk_dim[0]], 0) # in [width, height] format
            ImageDraw.Draw(msk).polygon(tuple(map(tuple, ij_arr[:,0:2])),
                                        outline=1, fill=1)
            msk = np.array(msk)
            msk_clouds = np.maximum(msk_clouds, msk)
    return msk_clouds

def read_detector_time_s2(path, fname='MTD_DS.xml', s2_df=None):
    """ get the detector metadata of the relative detector timing

    Parameters
    ----------
    path : string
        path where the meta-data is situated
    fname : string
        file name of the metadata.

    Returns
    -------
    det_time : numpy.array, numpy.datetime64
        detector time
    det_name : list, size=(13,)
        list of the different detector codes
    det_meta : numpy.array, size=(13,4)
        metadata of the individual detectors, each detector has the following:
            1) spatial resolution     [m]
            2) minimal spectral range [nm]
            3) maximum spectral range [nm]
            4) mean spectral range    [nm]
    line_period : float, numpy.timedelta64
        temporal sampling distance

    Notes
    -----
    The sensor blocks are arranged as follows, with ~98 pixels overlap:

        .. code-block:: text

             ┌-----┐   ┌-----┐   ┌-----┐   ┌-----┐   ┌-----┐   ┌-----┐
             |DET02|   |DET04|   |DET06|   |SCA08|   |SCA10|   |SCA12|
             └-----┘   └-----┘   └-----┘   └-----┘   └-----┘   └-----┘
         ┌-----┐   ┌-----┐   ┌-----┐   ┌-----┐   ┌-----┐   ┌-----┐
         |DET01|   |DET03|   |SCA05|   |SCA07|   |SCA09|   |SCA11|
         └-----┘   └-----┘   └-----┘   └-----┘   └-----┘   └-----┘

        .. code-block:: text

          a forward looking array, while band order is reverse for aft looking
          <-2592 10m pixels->
          ┌-----------------┐    #*# satellite
          |       B02       |     |
          ├-----------------┤     | flight
          |       B08       |     | direction
          ├-----------------┤     |
          |       B03       |     v
          ├-----------------┤
                 etc.
          the detector order is B02, B08, B03, B10, B04, B05,
          B11, B06, B07, B8A, B12, B01 and B09
    """
# example : line_period / np.timedelta64(1, 's')
    if s2_df is None:
        s2_df = list_central_wavelength_msi()

    det_time = np.zeros((13, 12), dtype='datetime64[ns]')
    det_name = [None] * 13
    det_meta = np.zeros((13, 4), dtype='float')

    root = get_root_of_table(path, fname)

    # image dimensions
    for meta in root.iter('Band_Time_Stamp'):
        bnd = int(meta.get('bandId')) # 0...12
        # <Spectral_Information bandId="8" physicalBand="B8A">
        for stamps in meta:
            det = int(stamps.get('detectorId')) # 1...12
            det_time[bnd,det-1] = np.datetime64(stamps[1].text, 'ns')

    # get line_period : not the ground sampling distance, but the temporal
    # sampling distance
    for elem in root.iter('LINE_PERIOD'):
       if elem.get('unit')=='ms': # datetime can only do integers...
           line_period = np.timedelta64(int(float(elem.text)*1e6), 'ns')

    for meta in root.iter('Spectral_Information'):
        bnd = int(meta.get('bandId'))
        bnd_name = meta.get('physicalBand')
        # at a zero if needed, this seems the correct naming convention
        if len(bnd_name)==2:
            bnd_name = bnd_name[0]+'0'+bnd_name[-1]
        det_name[bnd] = bnd_name
        det_meta[bnd,0] = float(meta[0].text) # resolution [m]
        det_meta[bnd,1] = float(meta[1][0].text) # min
        det_meta[bnd,2] = float(meta[1][1].text) # max
        det_meta[bnd,3] = float(meta[1][2].text) # mean

    # make selection based on the dataframe
    matchers = s2_df.index.values.tolist()
    matching = np.array([i for i,s in enumerate(det_name)
                         if any(xs in s for xs in matchers)])
    det_name = [det_name[i] for i in matching]
    det_time = det_time[matching,:]
    det_meta = det_meta[matching,:]
    return det_time, det_name, det_meta, line_period

def get_flight_bearing_from_detector_mask_s2(Det):
    """

    Parameters
    ----------
    Det : numpy.array, size=(m,n), ndim={2,3}, dtype=integer
        array with numbers of the different detectors in Sentinel-2

    Returns
    -------
    ψ : float, unit=degrees
        general bearing of the satellite
    """
    if Det.ndim==3:
        D = Det[:, :, 0]
    else:
        D = Det

    D_x = convolve2d(D, np.outer([+1., +2., +1.], [-1., 0., +1.]),
                     boundary='fill', mode='same')
    D_x[0, :], D_x[-1, :], D_x[:, 0], D_x[:, -1] = 0, 0, 0, 0
    D_tr = np.abs(D_x) >= 4

    L, Lab_num = label(D_tr, structure=[[1, 1, 1], [1, 1, 1], [1, 1, 1]])

    ψ = np.zeros(Lab_num)
    for i in range(Lab_num):
        [i_x, j_x] = np.where(L == i + 1)
        ψ[i] = np.rad2deg(np.arctan2(np.min(j_x) - np.max(j_x),
                                       np.min(i_x) - np.max(i_x)))
    ψ = np.rad2deg(np.arctan2(np.median(np.sin(np.deg2rad(ψ))),
                                np.median(np.cos(np.deg2rad(ψ)))))
    return ψ

# use the detector start and finish to make a selection for the flight line

def get_raw_imu_s2(ds_path, fname='MTD_DS.xml'):



    return

def get_raw_str_s2(ds_path, fname='MTD_DS.xml'):

    root = get_root_of_table(ds_path, fname)

    anci_info = None
    for child in root:
        if child.tag[-19:] == 'Ancillary_Data_Info':
            anci_info = child
    assert(anci_info is not None), ('metadata not in xml file')

    anci_data = None
    for child in anci_info:
        if child.tag == 'Attitudes':
            anci_data = child
    assert(anci_data is not None), ('metadata not in xml file')

    anci_list = None
    for child in anci_data:
        if child.tag == 'Raw_Attitudes':
            anci_list = child
    assert(anci_list is not None), ('metadata not in xml file')

    str_list = None
    for child in anci_list:
        if child.tag == 'STR_List':
            str_list = child
    assert(str_list is not None), ('metadata not in xml file')

    sat_time, sat_ang, sat_quat, sat_str = None, None, None, None
    for str in str_list:
        id = int(str.attrib['strId'][-1])

        att_list = None
        for entity in str:
            if entity.tag=='Attitude_Data_List':
                att_list = entity
        if att_list is not None:
            str_time = np.empty((len(att_list)), dtype='datetime64[ns]')
            str_ang = np.zeros((len(att_list), 3))
            str_quat = np.zeros((len(att_list), 4))
            for idx, att in enumerate(att_list):
                for field in att:
                    if field.tag == 'QUATERNION_VALUES':
                        str_quat[idx,:] = np.fromstring(field.text,
                                                        dtype=float, sep=' ')
                    elif field.tag == 'ANGULAR_RATE':
                        str_ang[idx] = np.fromstring(field.text,
                                                     dtype=float, sep=' ')
                    elif field.tag == 'GPS_TIME':
                        str_time[idx] = np.datetime64(field.text, 'ns')
            str_id = id*np.ones_like(str_ang[:,0])

        if sat_time is None:
            sat_time, sat_ang = str_time.copy(), str_ang.copy()
            sat_quat, sat_str = str_quat.copy(), str_id.copy()
        else:
            sat_time = np.hstack((sat_time, str_time))
            sat_ang = np.vstack((sat_ang, str_ang))
            sat_quat = np.vstack((sat_quat, str_quat))
            sat_str = np.hstack((sat_str, str_id))

    return sat_time, sat_ang, sat_quat, sat_str

def get_integration_and_sampling_time_s2(ds_path, fname='MTD_DS.xml',
                                         s2_dict=None): #todo: create s2_dict methodology
    """

    Parameters
    ----------
    ds_path : string
        location where metadata of datastrip is situated
    fname : string, default='MTD_DS.xml'
        metadata filename
    s2_dict : dictionary, default=None
        metadata dictionary, can be used instead of the string based method

    Returns
    -------

    Notes
    -----
    The metadata is scattered over the file structure of Sentinel-2, L1C

    .. code-block:: text

        * S2X_MSIL1C_20XX...
        ├ AUX_DATA
        │  └ DS_XXX_XXXX...
        │     └ QI_DATA
        │        └ MTD_DS.xml <- metadata about the data-strip
        ├ GRANULE
        │  └ L1C_TXXXX_XXXX...
        │     ├ AUX_DATA
        │     ├ IMG_DATA
        │     ├ QI_DATA
        │     └ MTD_TL.xml <- metadata about the tile
        ├ HTML
        ├ rep_info
        ├ manifest.safe
        ├ INSPIRE.xml
        └ MTD_MSIL1C.xml <- metadata about the product

    The following acronyms are used:

    - DS : datastrip
    - TL : tile
    - QI : quality information
    - AUX : auxiliary
    - MTD : metadata
    - MSI : multi spectral instrument
    - L1C : product specification,i.e.: level 1, processing step C

    """
    if isinstance(ds_path,dict):
        # only dictionary is given, which already has all metadata within
        s2_dict = ds_path
        assert ('MTD_DS_path' in s2_dict.keys()), 'please run get_s2_dict'

        root = get_root_of_table(s2_dict['MTD_DS_path'], 'MTD_MSIL1C.xml')
    else:
        root = get_root_of_table(ds_path, fname)

    integration_tim = np.zeros(13, dtype='timedelta64[ns]')
    # get line sampling interval, in nanoseconds
    for att in root.iter('Line_Period'.upper()):
        line_tim = np.timedelta64(int(float(att.text)*1e6), 'ns')

    # get line sampling interval, in nanoseconds
    for att in root.iter('Spectral_Band_Information'):
        band_idx = int(att.attrib['bandId'])
        for var in att:
            if var.tag == 'Integration_Time'.upper():
                int_tim = np.timedelta64(int(float(var.text)*1e6), 'ns')
        integration_tim[band_idx] = int_tim

    return line_tim, integration_tim



