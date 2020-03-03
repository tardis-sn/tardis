from tardis import run_tardis
from tardis.io.atom_data.util import download_atom_data
from dask.distributed import Client

# the data is automatically downloaded
download_atom_data('kurucz_cd23_chianti_H_He')

if __name__ == "__main__":
    client = Client()
    x = client.map(run_tardis, 'run_tardis_with_dask_config.yml')
    import pdb;

    pdb.set_trace()
    x.result()