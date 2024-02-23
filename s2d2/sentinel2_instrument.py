import numpy as np
import pandas as pd

from datetime import date

CENTRAL_WAVELENGTH_MSI = pd.DataFrame()

S2A_PLATFORM_SPECS = {
        'COSPAR': '2015-028A',
        'NORAD': 40697,
        'instruments': {'MSI'},
        'launch_date': date(2015, 6, 23),
        'mass': 1129.541,  # [kg]
        'inclination': 98.5621,  # https://www.n2yo.com/satellite/?s=40697
        'revolutions_per_day': 14.30824258387262,
        'J': np.array([[558, 30, -30], [30, 819, 30], [-30, 30, 1055]])
    }

S2B_PLATFORM_SPECS = {
        'COSPAR': '2017-013A',
        'NORAD': 42063,
        'instruments': {'MSI'},
        'launch_date': date(2017, 3, 7),
        'inclination': 98.5664,
        'revolutions_per_day': 14.30818491298178,
        'J': np.array([[558, 30, -30], [30, 819, 30], [-30, 30, 1055]])
    }

def dn_to_toa():
    # read_sentinel2.dn2toa_s2
    pass



