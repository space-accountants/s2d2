#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Space Accountants"
__license__ = "MIT License - You must cite this source"
__version__ = "202309"
__maintainer__ = "B. Altena"
__email__ = "info at space hyphen accountants dot eu"

import numpy as np

from s2d2.handler.xml import get_root_of_table
from s2d2.mapping_tools import ecef2llh, ecef2map

def get_flight_bearing_from_gnss_s2(path, spatialRef, rec_tim,
                                    fname='MTD_DS.xml'):
    """ get the direction/argument/heading of the Sentinel-2 acquisition

    Parameters
    ----------
    path : string
        directory where the meta-data is located
    spatialRef : osgeo.osr.SpatialReference
        projection system used
    rec_tim : {integer,float}
        time stamp of interest
    fname : string
        filename of the meta-data file

    Returns
    -------
    az : numpy.array, size=(m,1), float
        array with argument values, based upon the map projection given

    See Also
    --------
    get_s2_image_locations : to get the datastrip id

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
    """

    sat_tim,sat_xyz,_,_ = get_flight_path_s2(path, fname=fname)
    sat_xy = ecef2map(sat_xyz, spatialRef)

    dif_tim = sat_tim - rec_tim
    idx = np.argmin(np.abs(dif_tim))
    dif_xy = sat_xy[idx + 1] - sat_xy[idx]

    az = np.arctan2(dif_xy[0], dif_xy[1]) * 180 / np.pi
    return az

def get_flight_path_s2(ds_path, fname='MTD_DS.xml', s2_dict=None):
    """

    It is also possible to only give the dictionary, as this has such metadata
    within.

    Parameters
    ----------
    ds_path : string
        location of the metadata file
    fname : string, default='MTD_DS.xml'
        name of the xml-file that has the metadata
    s2_dict : dictonary
        metadata of the Sentinel-2 platform

    Returns
    -------
    sat_time : numpy.array, size=(m,1), dtype=np.datetime64, unit=ns
        time stamp of the satellite positions
    sat_xyz : numpy.array, size=(m,3), dtype=float, unit=meter
        3D coordinates of the satellite within an Earth centered Earth fixed
        (ECEF) frame.
    sat_err : numpy.array, size=(m,3), dtype=float, unit=meter
        error estimate of the 3D coordinates given by "sat_xyz"
    sat_uvw : numpy.array, size=(m,3), dtype=float, unit=meter sec-1
        3D velocity vectors of the satellite within an Earth centered Earth
        fixed (ECEF) frame.
    s2_dict : dictonary
        updated with keys: "gps_xyz", "gps_uvw", "gps_tim", "gps_err",
        "altitude", "speed"

    See Also
    --------
    get_flight_orientation_s2 : get the quaternions of the flight path

    Examples
    --------
    Following the file and metadata structure of scihub:

    >>> import os
    >>> import numpy as np
    >>> from s2d2.read_sentinel2 import list_central_wavelength_msi
    >>> from s2d2.handler.sentinel2 import get_s2_image_locations

    >>> s2_dir = '/data-dump/examples/'
    >>> s2_name = 'S2A_MSIL1C_20200923T163311_N0209_R140_T15MXV_20200923T200821.SAFE'
    >>> fname = os.path.join(s2_dir, s2_name, 'MTD_MSIL1C.xml')
    >>> s2_df = list_central_wavelength_s2()

    >>> s2_df, datastrip = get_s2_image_locations(fname, s2_df)
    >>> path_det = os.path.join(s2_dir, s2_name, 'DATASTRIP', datastrip[17:-7])

    >>> sat_tim, sat_xyz, sat_err, sat_uvw = get_flight_path_s2(path_det)

    Notes
    -----
    The metadata structure of MTD_DS looks like:

    .. code-block:: text

        * MTD_DS.xml
        └ n1:Level-1C_DataStrip_ID
           ├ n1:General_Info
           ├ n1:Image_Data_Info
           ├ n1:Satellite_Ancillary_Data_Info
           │  ├ Time_Correlation_Data_List
           │  ├ Ephemeris
           │  │  ├ GPS_Number_List
           │  │  ├ GPS_Points_List
           │  │  │  └ GPS_Point
           │  │  │     ├ POSITION_VALUES
           │  │  │     ├ POSITION_ERRORS
           │  │  │     ├ VELOCITY_VALUES
           │  │  │     ├ VELOCITY_ERRORS
           │  │  │     ├ GPS_TIME
           │  │  │     ├ NSM
           │  │  │     ├ QUALITY_INDEX
           │  │  │     ├ GDOP
           │  │  │     ├ PDOP
           │  │  │     ├ TDOP
           │  │  │     ├ NOF_SV
           │  │  │     └ TIME_ERROR
           │  │  └ AOCS_Ephemeris_List
           │  ├ Attitudes
           │  ├ Thermal_Data
           │  └ ANC_DATA_REF
           │
           ├ n1:Quality_Indicators_Info
           └ n1:Auxiliary_Data_Info

    The following acronyms are used:

    - AOCS : attitude and orbit control system
    - CAMS : Copernicus atmosphere monitoring service
    - DEM : digital elevation model
    - GDOP : geometric dilution of precision
    - GPS : global positioning system
    - GRI : global reference image
    - IERS : international earth rotation and reference systems service
    - IMT : instrument measurement time
    - NSM : navigation solution method
    - NOF SV : number of space vehicles
    - PDOP : position dilution of precision
    - TDOP : time dilution of precision

    """
    if isinstance(ds_path,dict):
        # only dictionary is given, which already has all metadata within
        s2_dict = ds_path
        assert ('MTD_DS_path' in s2_dict.keys()), 'please run get_s2_dict'

        root = get_root_of_table(s2_dict['MTD_DS_path'], 'MTD_MSIL1C.xml')
    else:
        root = get_root_of_table(ds_path, fname)

    for att in root.iter('GPS_Points_List'):
        sat_time = np.empty((len(att)), dtype='datetime64[ns]')
        sat_xyz = np.zeros((len(att), 3), dtype=float)
        sat_err = np.zeros((len(att), 3), dtype=float)
        sat_uvw = np.zeros((len(att), 3), dtype=float)
        counter = 0
        for idx, point in enumerate(att):
            # 'POSITION_VALUES' : tag
            xyz = np.fromstring(point[0].text, dtype=float, sep=' ')
            if point[0].attrib['unit'] == 'mm': xyz *= 1E-3 # convert to meters
            # 'POSITION_ERRORS'
            err = np.fromstring(point[1].text, dtype=float, sep=' ')
            if point[1].attrib['unit'] == 'mm': err *= 1E-3 # convert to meters
            # 'VELOCITY_VALUES'
            uvw = np.fromstring(point[2].text, dtype=float, sep=' ')
            if point[2].attrib['unit'] == 'mm/s': uvw *= 1E-3
            # convert to meters per second

            # 'GPS_TIME'
            gps_tim = np.datetime64(point[4].text, 'ns')

            # fill in the arrays
            sat_time[idx] = gps_tim
            sat_xyz[idx, :], sat_err[idx, :], sat_uvw[idx, :] = xyz, err, uvw
    if s2_dict is None:
        return sat_time, sat_xyz, sat_err, sat_uvw
    else: # include into dictonary
        s2_dict.update({'gps_xyz': sat_xyz, 'gps_uvw': sat_uvw,
                        'gps_tim': sat_time, 'gps_err': sat_err})
        # estimate the altitude above the ellipsoid, and platform speed
        llh = ecef2llh(sat_xyz)
        velo = np.linalg.norm(sat_uvw, axis=1)
        s2_dict.update({'altitude': np.squeeze(llh[:,-1]),
                        'velocity': np.squeeze(velo)})
        return s2_dict

def get_flight_orientation_s2(ds_path, fname='MTD_DS.xml', s2_dict=None):
    """ get the flight path and orientations of the Sentinel-2 satellite during
    acquisition.

    It is also possible to only give the dictionary, as this has such metadata
    within.

    Parameters
    ----------
    ds_path : string
        directory where the meta-data is located
    fname : string, default='MTD_DS.xml'
        filename of the meta-data file
    s2_dict : dictonary, default=None
        metadata of the Sentinel-2 platform

    Returns
    -------
    sat_time : numpy.array, size=(m,1), unit=nanosec
        satellite timestamps
    sat_angles : numpy.array, size=(m,3)
        satellite orientation angles, given in ECEF
    s2_dict : dictonary
        updated with keys: "time", "quat"

    See Also
    --------
    get_flight_path_s2 : get positions of the flight path
    """
    if isinstance(ds_path,dict):
        # only dictionary is given, which already has all metadata within
        s2_dict = ds_path
        assert ('MTD_DS_path' in s2_dict.keys()), 'please run get_s2_dict'

        root = get_root_of_table(s2_dict['MTD_DS_path'], 'MTD_MSIL1C.xml')
    else:
        root = get_root_of_table(ds_path, fname)

    for att in root.iter('Corrected_Attitudes'):
        sat_time = np.empty((len(att)), dtype='datetime64[ns]')
        sat_quat = np.zeros((len(att),4))
        counter = 0
        for idx, val in enumerate(att):

            for field in val:
                if field.tag == 'QUATERNION_VALUES':
                    sat_quat[idx, :] = np.fromstring(field.text,
                                                     dtype=float, sep=' ')
                elif field.tag == 'GPS_TIME':
                    sat_time[idx] = np.datetime64(field.text, 'ns')
    if s2_dict is None:
        return sat_time, sat_quat
    else: # include into dictonary
        s2_dict.update({'imu_quat': sat_quat, 'imu_time': sat_time})
        return s2_dict

#todo: def get_intrinsic_temperatures_s2(ds_path, fname='MTD_DS.xml'):
