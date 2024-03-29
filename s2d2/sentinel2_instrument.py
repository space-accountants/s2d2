import numpy as np
import pandas as pd

from datetime import date

center_wavelength = {
    "B01": 443,
    "B02": 492,
    "B03": 560,
    "B04": 665,
    "B05": 704,
    "B06": 741,
    "B07": 783,
    "B08": 833,
    "B8A": 865,
    "B09": 945,
    "B10": 1374,
    "B11": 1614,
    "B12": 2202,
}
# convert from nm to µm
center_wavelength = {k: v / 1E3 for k, v in center_wavelength.items()}

full_width_half_max = {
    "B01": 21,
    "B02": 66,
    "B03": 36,
    "B04": 31,
    "B05": 15,
    "B06": 15,
    "B07": 20,
    "B08": 106,
    "B8A": 21,
    "B09": 20,
    "B10": 31,
    "B11": 91,
    "B12": 175,
}
# convert from nm to µm
full_width_half_max = {k: v / 1E3 for k, v in full_width_half_max.items()}

bandid = {
    "B01": 0,
    "B02": 1,
    "B03": 2,
    "B04": 3,
    "B05": 4,
    "B06": 5,
    "B07": 6,
    "B08": 7,
    "B8A": 8,
    "B09": 9,
    "B10": 10,
    "B11": 11,
    "B12": 12,
}
gsd = {
    "B01": 60,
    "B02": 10,
    "B03": 10,
    "B04": 10,
    "B05": 20,
    "B06": 20,
    "B07": 20,
    "B08": 10,
    "B8A": 20,
    "B09": 60,
    "B10": 60,
    "B11": 20,
    "B12": 20,
}

along_pixel_size = {
    "B01": 45.,
    "B02": 7.5,
    "B03": 7.5,
    "B04": 7.5,
    "B05": 15.,
    "B06": 15.,
    "B07": 15.,
    "B08": 7.5,
    "B8A": 15.,
    "B09": 45.,
    "B10": 15.,
    "B11": 15.,
    "B12": 15.,
}

across_pixel_size = {
    "B01": 15.,
    "B02": 7.5,
    "B03": 7.5,
    "B04": 7.5,
    "B05": 15.,
    "B06": 15.,
    "B07": 15.,
    "B08": 7.5,
    "B8A": 15.,
    "B09": 15.,
    "B10": 15.,
    "B11": 15.,
    "B12": 15.,
}

# given as base-to-height
crossdetector_parallax = {
    "B01": 0.055,
    "B02": 0.022,
    "B03": 0.030,
    "B04": 0.034,
    "B05": 0.038,
    "B06": 0.042,
    "B07": 0.046,
    "B08": 0.026,
    "B8A": 0.051,
    "B09": 0.059,
    "B10": 0.030,
    "B11": 0.040,
    "B12": 0.050,
}
# convert to angles in degrees
crossdetector_parallax = \
    {k: 2 * np.tan(v / 2) for k, v in crossdetector_parallax.items()}

common_name = {
    "B01": 'coastal',
    "B02": 'blue',
    "B03": 'green',
    "B04": 'red',
    "B05": 'rededge',
    "B06": 'rededge',
    "B07": 'rededge',
    "B08": 'nir',
    "B8A": 'nir08',
    "B09": 'nir09',
    "B10": 'cirrus',
    "B11": 'swir16',
    "B12": 'swir22'
}
sensor_array = {
    "B01": 'VNIR',
    "B02": 'VNIR',
    "B03": 'VNIR',
    "B04": 'VNIR',
    "B05": 'VNIR',
    "B06": 'VNIR',
    "B07": 'VNIR',
    "B08": 'VNIR',
    "B8A": 'VNIR',
    "B09": 'VNIR',
    "B10": 'SWIR',
    "B11": 'SWIR',
    "B12": 'SWIR',
}

MSI_SPECIFICS = pd.DataFrame({
    "center_wavelength":
        pd.Series(center_wavelength, dtype=np.dtype('float')),
    "full_width_half_max":
        pd.Series(full_width_half_max, dtype=np.dtype('float')),
    "gsd":
        pd.Series(gsd, dtype=np.dtype('float')),
    "across_pixel_size":
        pd.Series(across_pixel_size, dtype=np.dtype('float')),
    "along_pixel_size":
        pd.Series(along_pixel_size, dtype=np.dtype('float')),
    "common_name":
        pd.Series(common_name, dtype=np.dtype('str')),
    "sensor_array":
        pd.Series(sensor_array, dtype=np.dtype('str')),
    "bandid":
        pd.Series(bandid, dtype=np.dtype('int64')),
    "crossdetector_parallax":
        pd.Series(crossdetector_parallax, dtype=np.dtype('float'))
})


def dn_to_toa(img_dn):
    """
    convert the digital numbers of Sentinel-2 to top of atmosphere (TOA), see
    for more details [wwwS2L1C]_.

    Parameters
    ----------
    img_dn : numpy.ndarray, dim={2,3}, size=(m,n)
        grid with intensities

    Returns
    -------
    numpy.ndarray, dim={2,3}, size=(m,n)
        grid with top of atmosphere reflectances

    Notes
    -----
    .. [wwwS2L1C] https://sentinel.esa.int/web/sentinel/technical-guides/sentinel-2-msi/level-1c/algorithm
    """
    img_toa = np.divide(img_dn.astype('float'), 1E4)
    return img_toa
