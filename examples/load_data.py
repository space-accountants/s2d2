import os

from s2d2.sentinel2_product import Sentinel2Product
from s2d2.sentinel2_band import Sentinel2Band

dat_path = 'data/S2A_MSIL1C_20170225T155221_N0204_R054_T18TWL_20170225T155545.SAFE'
gml_path = '/home/basal/parcels/s2d2/data/S2A_MSIL1C_20170225T155221_N0204_R054_T18TWL_20170225T155545.SAFE/GRANULE/L1C_T18TWL_A008774_20170225T155545/QI_DATA'
img_path = '/home/basal/parcels/s2d2/data/S2A_MSIL1C_20170225T155221_N0204_R054_T18TWL_20170225T155545.SAFE/GRANULE/L1C_T18TWL_A008774_20170225T155545/IMG_DATA'

#s2b = Sentinel2Band('B04',32618,(499980.0, 10.0, 0.0, 4600020.0, 0.0, -10.0),10980,10980)
#s2b.read_band(img_path, fname='T18TWL_20170225T155221_B04.jp2')
#s2b.read_detector_gml(gml_path)


s2d2 = Sentinel2Product(dat_path)
s2d2.load_metadata()
s2d2.specify_bands_of_interest(query='gsd==10')
s2d2.update_bands_metadata()

print('.')

smaller_than = lambda x: [i for i in range(len(x)) if x[i]<1.]
equal_to = lambda x: [i for i in range(len(x)) if x[i]==10.]
