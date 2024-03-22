from typing import Callable, Optional

from .sentinel2_instrument import MSI_SPECIFICS


def get_bands_of_interest(query: Optional[str] = None,
                          subject: Optional[str] = None,
                          condition: Optional[Callable] = None,
                          selection: Optional[list] = None) -> list[str]:
    """

    Parameters
    ----------
    subject : str
        the characteristic to base the selection on. The following specific
        information about the MSI instrument that is onboard Sentinel-2 can
        be selected:

            * "center_wavelength", unit=µm : central wavelength of the band
            * "full_width_half_max", unit=µm : extent of the spectral
                sensitivity
            * "resolution", unit=m : spatial resolution of a pixel
            * "along_pixel_size", unit=µm : physical size of the sensor
            * "across_pixel_size", unit=µm : size of the photosensative sensor
            * "crossdetector_parallax", unit=degress : in along-track direction
            * "common_name" : general name of the band, if applicable
            * "bandid" : number for identification in the meta data
    """
    band_list = []
    if query is not None:
        band_list = MSI_SPECIFICS.query(query).get('bandid')
        return band_list

    if subject is None:
        band_list = MSI_SPECIFICS.get('bandid')
        return band_list

    assert subject in MSI_SPECIFICS.keys(), \
        (f'subject does not seem to be present, please provide one of the'
         f'following: {list(MSI_SPECIFICS.keys())}')

    coi = MSI_SPECIFICS.get(subject)
    if selection is None:
        idx = condition(coi)
        band_list = MSI_SPECIFICS.iloc[idx].get('bandid')
    else:
        band_list = MSI_SPECIFICS.loc[coi.isin(selection)].get('bandid')

    return band_list
