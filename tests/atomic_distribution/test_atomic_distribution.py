import numpy as np
import xarray as xr
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"

files_hdf5 = [
    "outfile.mean_square_displacement.hdf5",
    "outfile.pair_distribution.hdf5",
    "outfile.powder_diffraction.hdf5",
]


def test_mean_square_displacement(file="outfile.mean_square_displacement"):
    file_ref = folder / file
    file_new = parent / file

    data_ref = np.loadtxt(file_ref)
    data_new = np.loadtxt(file_new)

    np.testing.assert_allclose(data_ref, data_new, err_msg=file_new.absolute())


def test_hdf5(files=files_hdf5, decimals=7):
    for file in files:
        file_ref = folder / file
        file_new = parent / file

        ds_ref = xr.load_dataset(file_ref)
        ds_new = xr.load_dataset(file_new)

        for var in ds_ref.data_vars:
            x = ds_ref[var].round(decimals=decimals)
            y = ds_new[var].round(decimals=decimals)
            np.testing.assert_allclose(x, y, err_msg=var)


if __name__ == "__main__":
    test_mean_square_displacement()
    test_hdf5()
