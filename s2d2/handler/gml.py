import numpy as np

# todo: robustify
def get_xy_polygon_from_gml(gml_struct, idx: int) -> tuple[np.ndarray, int]:
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
