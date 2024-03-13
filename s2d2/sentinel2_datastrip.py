
import xml.etree.ElementTree as ElementTree
import numpy as np
import pandas as pd

from .handler.xml import get_root_of_table, get_branch
from .mapping_tools import ecef2llh
from .typing import Path


class Sentinel2Datastrip:
    def __init__(self, path: Path) -> None:
        # add as optional paths of all files used here?
        self.path = path

        # parent specifics
        self.spacecraft = None

        self.datatake_id = None
        self.orbit = None
        self.orbit_counter = None

        #
        self.tile_list = None

        self.gps_flightpath = None
        self.attitudes_corrected = None
        self.attitudes_raw = None

        self.detector_time = None
        self.sat_time = None
        self.sat_ang = None
        self.sat_quat = None
        self.sat_str = None
        self.integration_time = None
        self.sampling_time = None

    def load_metadata(self) -> None:
        root = get_root_of_table(self.path, fname='MTD_DS.xml')

        gen_info = get_branch(root, 'General_Info')
        self._get_ids_from_xmltree(gen_info)
        self._get_abs_orbit_from_xmltree(gen_info)

        img_info = get_branch(root, 'Image_Data_Info')
        tile_info = get_branch(img_info, 'Tiles_Information')
        tile_list = get_branch(tile_info, 'Tile_List')
        self._get_tile_list_from_xmltree(tile_list)

        sen_conf = get_branch(img_info, 'Sensor_Configuration')

        sat_anc = get_branch(root, 'Satellite_Ancillary_Data_Info')
        sat_eph = get_branch(sat_anc, 'Ephemeris')

        gps_pts = get_branch(sat_eph, 'GPS_Points_List')
        self._get_gps_flightpath_from_xmltree(gps_pts)

        aocs_pts = get_branch(sat_eph, 'AOCS_Ephemeris_List')
        self._get_aocs_flightpath_from_xmltree(aocs_pts)

        sat_att = get_branch(sat_anc, 'Attitudes')
        cor_att = get_branch(sat_att, 'Corrected_Attitudes')
        self._get_corrected_attitudes_from_xmltree(cor_att)

        # read_sentinel2.read_detector_time_s2
        self.detector_time = ...
        # read_sentinel2.get_raw_str_s2
        # also in read_sentinel2.get_flight_orientation_s2 ?
        self.sat_time = ...  # full name would make it clearer
        self.sat_ang = ...  # full name would make it clearer
        self.sat_quat = ...  # full name would make it clearer
        self.sat_str = ...  # full name would make it clearer
        # read_sentinel2.get_integration_and_sampling_time_s2
        self.integration_time = ...
        self.sampling_time = ...
        # sensor_readings_sentinel2.get_flight_path_s2
        self.sat_time = ... # full name would make it clearer - is it same as from read_sentinel2.get_raw_str_s2?
        self.sat_xyz = ... # full name would make it clearer
        self.sat_err = ... # full name would make it clearer
        self.sat_uvw = ... # full name would make it clearer


    def _get_ids_from_xmltree(self,
                              general_info: ElementTree.Element) -> None:
        for field in general_info:
            if not (field.tag == 'Datatake_Info'): continue

            self.datatake_id = field.attrib['datatakeIdentifier']
            for instance in field:
                if instance.tag == 'SPACECRAFT_NAME':
                    self.spacecraft = instance.text[-1].upper()
                elif instance.tag == 'SENSING_ORBIT_NUMBER':
                    self.orbit = int(instance.text)

    def _get_abs_orbit_from_xmltree(self,
                                    general_info: ElementTree.Element) -> None:
        for field in general_info:
            if not (field.tag == 'Downlink_Info'): continue
        for instance in field:
            if instance.tag == 'DOWNLINK_ORBIT_NUMBER':
                self.orbit_counter = int(instance.text)

    def _get_tile_list_from_xmltree(self, tile_list: ElementTree.Element) -> None:
        tl = []
        for tile in tile_list:
            tl.append(tile.attrib['tileId'])
        self.tile_list = tl


    def _get_gps_flightpath_from_xmltree(self, gps_list: ElementTree.Element) -> None:
        """
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
               │  │  └ AOCS_Ephemeris_List
               │  │  │  └ AOCS_Ephmeris
               │  │  │     ├ VALID_FLAG
               │  │  │     ├ OPSOL_QUALITY
               │  │  │     ├ POSITION_VALUES
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
               │  ├ Attitudes
               │  ├ Thermal_Data
               │  └ ANC_DATA_REF
               │
               ├ n1:Quality_Indicators_Info
               └ n1:Auxiliary_Data_Info
        """

        df = pd.DataFrame(columns=["pos_x", "pos_y", "pos_z",
                                   "pos_x_err", "pos_y_err", "pos_z_err",
                                   "vel_x", "vel_y", "vel_z",
                                   "vel_x_err", "vel_y_err", "vel_z_err",
                                   "time_gps", "time_err"
                                   "gdop", "pdop", "tdop", "sats"])

        for instance in gps_list:
            entry = {}
            for elem in instance:
                if elem.tag == 'POSITION_VALUES':
                    xyz = np.fromstring(elem.text, dtype=float, sep=' ')
                    if elem.attrib['unit'] == 'mm':
                        xyz *= 1E-3  # convert to meters
                    entry.update({"pos_x": [xyz[0]], "pos_y": [xyz[1]], "pos_z": [xyz[2]]})
                elif elem.tag == 'POSITION_ERRORS':
                    xyz_err = np.fromstring(elem.text, dtype=float, sep=' ')
                    if elem.attrib['unit'] == 'mm':
                        xyz_err *= 1E-3  # convert to meters
                    entry.update({"pos_x_err": [xyz_err[0]], "pos_y_err": [xyz_err[1]], "pos_z_err": [xyz_err[2]]})
                if elem.tag == 'VELOCITY_VALUES':
                    uvw = np.fromstring(elem.text, dtype=float, sep=' ')
                    if elem.attrib['unit'] == 'mm/s':
                        uvw *= 1E-3  # convert to meters per second
                    entry.update({"vel_x": [xyz[0]], "vel_y": [xyz[1]], "vel_z": [xyz[2]]})
                elif elem.tag == 'VELOCITY_ERRORS':
                    uvw_err = np.fromstring(elem.text, dtype=float, sep=' ')
                    if elem.attrib['unit'] == 'mm/s':
                        uvw_err *= 1E-3  # convert to meters
                    entry.update({"vel_x_err": [uvw_err[0]], "vel_y_err": [uvw_err[1]], "vel_z_err": [uvw_err[2]]})
                elif elem.tag == 'GPS_TIME':
                    entry.update({"time_gps": [np.datetime64(elem.text, 'ns')]})
                elif elem.tag == 'TIME_ERROR':
                    entry.update({"time_err": [np.datetime64(elem.text, elem.attrib['unit'])]})
                elif elem.tag == 'GDOP':
                    entry.update({"gdop": [int(elem.text)]})
                elif elem.tag == 'PDOP':
                    entry.update({"pdop": [int(elem.text)]})
                elif elem.tag == 'TDOP':
                    entry.update({"tdop": [int(elem.text)]})
            # combine into dataframe
            df = pd.concat([df, pd.DataFrame.from_dict(entry)], ignore_index=True)
        self.gps_flightpath = df

    def _get_aocs_flightpath_from_xmltree(self, aocs_list: ElementTree.Element) -> None:
        df = pd.DataFrame(columns=["pos_x", "pos_y", "pos_z",
                                   "vel_x", "vel_y", "vel_z",
                                   "time_gps", "orb_ang"])
        for instance in aocs_list:
            entry = {}
            for elem in instance:
                if elem.tag == 'POSITION_VALUES':
                    xyz = np.fromstring(elem.text, dtype=float, sep=' ')
                    if elem.attrib['unit'] == 'mm':
                        xyz *= 1E-3  # convert to meters
                    entry.update({"pos_x": [xyz[0]], "pos_y": [xyz[1]], "pos_z": [xyz[2]]})
                if elem.tag == 'VELOCITY_VALUES':
                    uvw = np.fromstring(elem.text, dtype=float, sep=' ')
                    if elem.attrib['unit'] == 'mm/s':
                        uvw *= 1E-3  # convert to meters per second
                    entry.update({"vel_x": [xyz[0]], "vel_y": [xyz[1]], "vel_z": [xyz[2]]})
                elif elem.tag == 'ORBIT_ANGLE':
                    orb_ang = [float(elem.text)]
                    if elem.attrib['unit'] == 'rad':
                        orb_ang = np.rad2deg(orb_ang)[0]  # convert to degrees
                    entry.update({"orb_ang": [orb_ang]})
                elif elem.tag == 'GPS_TIME':
                    entry.update({"time_gps": [np.datetime64(elem.text, 'ns')]})
            # combine into dataframe
            df = pd.concat([df, pd.DataFrame.from_dict(entry)], ignore_index=True)
        self.aocs_flightpath = df


    def _get_corrected_attitudes_from_xmltree(self, attitudes_list: ElementTree.Element) -> None:
        """
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
               │  ├ Attitudes
               │  │  └ Corrected_Attitudes
               │  │  │  └ Values
               │  │  │     ├ QUATERNION_VALUES
               │  │  │     ├ QUATERNION_VALIDITY
               │  │  │     ├ GPS_TIME
               │  │  │     ├ INUSE_FLAGS
               │  │  │     ├ AOCS_MODE
               │  │  │     ├ AOCS_SUBMODE
               │  │  │     ├ INNOVATION_STR1
               │  │  │     ├ INNOVATION_STR2
               │  │  │     └ ATTITUDE_QUALITY_INDICATOR
               │  ├ Thermal_Data
               │  └ ANC_DATA_REF
               ├ n1:Quality_Indicators_Info
               └ n1:Auxiliary_Data_Info

        These are as follows:
            - QUATERNION_VALUES: float float float float
                space separated list of 4 quaternion values ordered as qv1 qv2 qv3 qs
            - QUATERNION_VALIDITY: bool
                true if quaternion is valid
            - GPS_TIME: yyyy-mm-ddThh:mm:ss.s
            - INUSE_FLAGS: bool bool bool bool bool bool bool ...
                list of 11 boolean flags separated by whitespace for
                STR1 STR2 STR3 GPSR-A GPSR-B VCU-A VCU-B IMU-1 IMU-2 IMU-3 IMU-4
            - AOCS_MODE: int
            - AOCS_SUBMODE: int
            - INNOVATION_STR1: float float float
                difference between GSE filter estimate and second in-use STR measurement
            - INNOVATION_STR2: float float float
                difference between GSE filter estimate and second in-use STR measurement
            - ATTITUDE_QUALITY_INDICATOR: str
        """
        df = pd.DataFrame(columns=["quat", "time_gps"])
        for instance in attitudes_list:
            entry = {}
            for elem in instance:
                if elem.tag == 'QUATERNION_VALUES':
                    quat = np.fromstring(elem.text, dtype=float, sep=' ')
                    entry.update({"quat": [quat.tolist()]})
                elif elem.tag == 'GPS_TIME':
                    entry.update({"time_gps": [np.datetime64(elem.text, 'ns')]})
            # combine into dataframe
            df = pd.concat([df, pd.DataFrame.from_dict(entry)], ignore_index=True)
        self.attitudes_corrected = df



#    def _get_corrected_flightpath(self):

    def get_altitude(self):
        # sensor_readings_sentinel2.get_flight_path_s2
        # estimate the altitude above the ellipsoid
        llh = ecef2llh(self.sat_xyz)
        return np.squeeze(llh[:,-1])

    def get_velocity(self):
        # sensor_readings_sentinel2.get_flight_path_s2
        # estimate platform speed
        velo = np.linalg.norm(self.sat_uvw, axis=1)
        return np.squeeze(velo)

    #TODO def get_intrinsic_temperatures_s2(self):
