import os
import logging
import requests

from tqdm.autonotebook import tqdm
from tardis.io.util import get_internal_data_path
from tardis.io.config_internal import get_data_dir
import yaml

logger = logging.getLogger(__name__)

def get_atomic_repo_config():
    atomic_repo_fname = get_internal_data_path('atomic_data_repo.yml')
    return yaml.load(open(atomic_repo_fname))

def download_from_url(url, dst):
    """
    kindly used from https://gist.github.com/wy193777/0e2a4932e81afc6aa4c8f7a2984f34e2
    @param: url to download file
    @param: dst place to put the file
    """

    file_size = int(requests.head(url).headers["Content-Length"])
    if os.path.exists(dst):
        first_byte = os.path.getsize(dst)
    else:
        first_byte = 0
    if first_byte >= file_size:
        return file_size
    header = {"Range": "bytes=%s-%s" % (first_byte, file_size)}
    pbar = tqdm(
        total=file_size, initial=first_byte,
        unit='B', unit_scale=True, desc=url.split('/')[-1])
    req = requests.get(url, headers=header, stream=True)
    with(open(dst, 'ab')) as f:
        for chunk in req.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)
                pbar.update(1024)
    pbar.close()
    return file_size



def download_atomic_data(atomic_data_name=None):
    """
    Download the atomic data from the repository
    Parameters
    ----------
    atomic_data_name: str
        if None

    Returns
    -------
        : None

    """
    atomic_repo = get_atomic_repo_config()

    if atomic_data_name is None:
        atomic_data_name = atomic_repo['default']

    if atomic_data_name not in atomic_repo:
        raise ValueError('Atomic Data name {0} not known'.format(atomic_data_name))
    dst_dir = os.path.join(get_data_dir(), '{0}.h5'.format(atomic_data_name))
    src_url = atomic_repo[atomic_data_name]['url']
    logger.info('Downloading atomic data from {0} to {1}'.format(src_url, dst_dir))
    download_from_url(src_url, dst_dir)