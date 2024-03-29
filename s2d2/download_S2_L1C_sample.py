import os
import random
import shutil
import zipfile

from urllib.request import urlretrieve

ZEONDO_URL = 'https://zenodo.org/records/5221412/files/'

SAMPLES = ['beijing.zip', 'berlin.zip', 'bogota.zip',
           'cairo.zip', 'dakar.zip', 'dhaka.zip', 'jakarta.zip',
           'johannesburg.zip', 'la_paz.zip', 'las_vegas.zip',
           'mexico_city.zip', 'moscow.zip', 'nairobi.zip', 'new_delhi.zip',
           'new_york.zip', 'ottawa.zip', 'sao_paulo.zip', 'stockholm.zip',
           'sydney.zip', 'tokyo.zip']

def get_random_sample_file() -> str:
    """ provide a random zip-file name of the collection

    Returns
    -------
    str
        for example "beijing.zip"
    """
    rand_int = random.randint(0, len(SAMPLES)-1)
    return SAMPLES[rand_int]

def _make_full_zenodo_url(region_name: str) -> str:
    assert isinstance(region_name, str), 'please provide a string'

    region_name = region_name.lower().strip()
    if not region_name.endswith('.zip'): region_name += '.zip'

    assert region_name in SAMPLES, print(f'{region_name} is not present')
    return ZEONDO_URL + region_name

def download_S2_L1C_scene(region_name: str, out_dir: str = None) -> str:

    if out_dir is None: out_dir = os.getcwd()

    full_url = _make_full_zenodo_url(region_name)
    file_name = os.path.basename(full_url)
    file_path = os.path.join(out_dir, file_name)

    # download data
    print('INFO: start downloading data')
    urlretrieve(full_url, file_path)
    print('INFO: downloading done')
    return file_path

def unpack_S2_L1C_scene(zip_path: str, delete_zip: bool = True) -> str:

    region_name = os.path.splitext(os.path.basename(zip_path))[0]
    print(f'INFO: unpacking data from: {region_name}')

    dat_dir = os.path.dirname(zip_path)
    with zipfile.ZipFile(zip_path, mode="r") as archive:

        for member in archive.namelist():
            rel_path = os.path.join(*(member.split(os.path.sep)[1:]))
            file_name = os.path.basename(member)
            if not file_name: continue

            source = archive.open(member)
            new_path = os.path.join(dat_dir, rel_path)
            new_dir = os.path.dirname(new_path)
            os.makedirs(new_dir, exist_ok = True)
            os.chdir(new_dir)
            target = open(new_path, "wb")
            with source, target:
                shutil.copyfileobj(source, target)
    s2_dir = os.path.join(dat_dir, rel_path.split(os.path.sep)[0])

    # delete zip-file
    if delete_zip:
        print('INFO: deleting zip-file')
        os.remove(zip_path)
    return s2_dir

def main() -> None:
    region_name = get_random_sample_file()

    dat_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    dat_dir = os.path.join(dat_dir, 'data')
    os.makedirs(dat_dir, exist_ok=True)

    zip_path = download_S2_L1C_scene(region_name, out_dir=dat_dir)
    s2_dir = unpack_S2_L1C_scene(zip_path)
    print(f"INFO: data is situated in this folder: {s2_dir}")
    return

if __name__ == '__main__':
    main()
