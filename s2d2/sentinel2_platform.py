import pandas as pd
import numpy as np

from datetime import date

S2_PLATFORM_SPECS = {
        'COSPAR': ['2015-028A', '2017-013A'],
        'NORAD': [40697, 42063],
        'ESOC': [266, 267],
        'launch_date': [date(2015, 6, 23), date(2017, 3, 7)],
        'J': [ [[558, 30, -30], [30, 819, 30], [-30, 30, 1055]], [[558, 30, -30], [30, 819, 30], [-30, 30, 1055]]]
}

S2_GNSS_SPECS = {
    'X': [232., 232.],
    'Y': [227.5, -72.5],
    'Z': [-810., -810.],
    'R': [[[-0.966, 0., -0.259], [0., 1., 0.], [0.259, 0., -0.966]],
         [[-0.966, 0., -0.259], [0., 1., 0.], [0.259, 0., -0.966]]]
}

class Sentinel2Platform:
    def __init__(self, spacecraft: str) -> None:
        df = pd.DataFrame.from_dict(data=S2_PLATFORM_SPECS, orient='index', columns=['A', 'B'])
        assert spacecraft in df.keys(), f'please provide correct spacecraft id, that is: {list(df.columns.values)}'

        self.id_norad = df[spacecraft]['NORAD']
        self.id_cospar = df[spacecraft]['COSPAR']
        self.id_esoc = df[spacecraft]['ESOC']

        self.mass = None
        self.drag = df[spacecraft]['J']

class solarpanel:
    def __init__(self) -> None:
        self.dimensions = None
        self.orientation = None
        self.angularsuntrackrate = 0.06 # [deg / sec]
        self.angularrewindrate = 0.205 # [deg / sec]
        self.angularrange = 239 # [deg]

class gnss:
    def __init__(self, id: str) -> None:
        df = pd.DataFrame.from_dict(S2_GNSS_SPECS, orient="index", columns=['GPS-A', 'GPS-B'])
        self.dimensions = [df[id]['X'], df[id]['Y'], df[id]['Z']]
        self.orientation = df[id]['R']

def list_jitter_frequencies_s2():
    """

    References
    ----------
    .. [Mi21] Mingione, "Solar array drive mechanism disturbance estimation:
              a robust approach applied to Sentinel-2 mission", Msc thesis of
              Politecnico di Milano, April 2021.
    .. [Al21] Alazard et al., "Characterization of SADM induced disturbances
              and their effects on spacecraft pointing errors" 11th
              international ESA conference on guidance, navigation & control
              systems, 2021.

    """
    d = {'imperfection': [1, 2, 3, 4, 5, 6],
         'angular rates': [13616., 184., 197.3, 14800., 197.3, 185.],
         'gear pair': [(2,3), (2,3), (2,3), (3,4), (3,4), (3,4)],
         'cause': ['gear frequency', 'a tooth on body 2', 'a tooth on body 3',
                   'gear frequency', 'a tooth on body 3', 'a tooth on body 4'],
         'number of sources': [1, 74, 69, 1, 75, 80]}
    df = pd.DataFrame.from_dict(d)
    return df

# solar panel
# angular angular velocity is 0.06°/s during sun-tracking
# and 0.205°/s for the rewind.
# The rotation angle in both directions required for the nominal operation of the solar array is ±119.5°

# mass history file
# https://sentinels.copernicus.eu/documents/d/sentinel/s2b-mhf

# historic outage files
# https://sentinels.copernicus.eu/documents/d/sentinel/s2a-out

# manuever file
#https://sentinels.copernicus.eu/documents/d/sentinel/s2a-man

# The Sentinel-2 spacecraft reference frame is defined as:
#  - The origin is located in the plane of attachment to the launcher and in the centre of the attachment ring.
#  - The X axis is perpendicular to the Satellite/Launcher separation plane, pointing positively from the separation plane towards the Satellite.
#  - The Y axis completes the right-handed orthogonal system.
#  - The Z axis is pointing towards the satellite side which is nominally nadir pointing.
# The Sentinel-2 nadir reference frame is defined as:
#
#
# The attitude of Sentinel-2 expressed in J2000 reference frame is composed of two rotations:
# - Rotation from J2000 to the Nadir reference Frame as described in section 3.7 of [RD.3].
# - Rotation (around Z axis, yaw rotation) from the Nadir reference frame to the Satellite Body Fixed Reference Frame as described in section 3.8 of [RD.3].
# According to [RD.3] the required spacecraft pointing and attitude during extended observation is a
# specified roll angle up to ±20.38 deg. The quaternions are reflecting this attitude mode if extended
# observations are performed

def list_gps_antenna_coords():
    d = {'X': [232., 232.],
         'Y': [227.5, -72.5],
         'Z': [-810., -810.],
         'R': [[[-0.966, 0., -0.259], [0., 1., 0.], [0.259, 0., -0.966]],
               [[-0.966, 0., -0.259], [0., 1., 0.], [0.259, 0., -0.966]]]
         }
    df = pd.DataFrame.from_dict(d, orient="index", columns=['GPS-A', 'GPS-B'])
    return df

# https://documentation.dataspace.copernicus.eu/Data/SentinelMissions/Sentinel2.html#sentinel-2-precise-orbit-determination-pod-products

Manoeuvre history file format description (header)
Key Type Description
Epoch Epoch, a23 File last update epoch as YYYY/MM/DD-HH:MM:SS.SSS
ESOC ID Integer,i3 Satellite ESOC ID (266 for S2A and 267 for S2B)

Manoeuvre history file format description (body)
Epoch Epoch, a23 Burn start or stop time (UTC) as YYYY/MM/DD-HH:MM:SS.SSS
Acceleration Real, f15.8 First component of acceleration in km/s2
Acceleration Real, f15.8 Second component of acceleration in km/s2
Acceleration Real, f15.8 Third component of acceleration in km/s2
Manoeuvre Integer, i1 Record flag. If >0, this is a manoeuvre start record. If 0, it is a manoeuvre end
start-end flag record.
On a start record:
- value 1 indicates that the components are radial, along-track and cross-
track respectively
- value 2 that they are along the J2000.0 X-, Y- and Z- axes respectively

Example:
2018/04/11-10:00:03.026 266
2016/04/29-07:25:31.282 -0.62906997D-08 0.89614348D-10-0.50554106D-06 1
2016/04/29-07:27:12.407 -0.62906997D-08 0.89614348D-10-0.50554106D-06 0
2016/04/29-15:39:26.552 -0.63265823D-08-0.21231811D-08-0.50738033D-06 1
2016/04/29-15:41:07.552 -0.63265823D-08-0.21231811D-08-0.50738033D-06 0

Outages file format description (header)
Key Type Description
Epoch Epoch, a23 File last update epoch as YYYY/MM/DD-HH:MM:SS.SSS
File Title String, a12 OUTAGES FILE
ESOC ID Integer,i3 Satellite ESOC ID (266 for S2A and 267 for S2B)
Table 4-5: Outages file format description (body)
Key Type Description
Epoch Epoch, a23 Outage start epoch as YYYY/MM/DD-HH:MM:SS.SSS
Epoch Epoch, a23 Outage end epoch as YYYY/MM/DD-HH:MM:SS.SSS
Outage
type
String, a3 Type of outage. It may be: input gap (GAP), manoeuvre (MAN) or a combination of
both (MIX)
Example:
2018/04/18-03:10:06.000 OUTAGES FILE 266
2000/01/01-00:00:00.000 2016/04/27-04:43:04.000 GAP
2016/04/28-16:44:34.000 2016/04/29-10:44:14.000 GAP
2016/04/29-15:39:26.000 2016/04/29-15:41:07.000 MAN
2016/05/03-12:51:38.000 2016/05/03-12:53:19.000 MAN
2016/05/05-15:03:47.000 2016/05/06-07:55:35.000 GAP
2016/05/09-14:36:29.000 2016/05/09-23:59:59.000 MIX
