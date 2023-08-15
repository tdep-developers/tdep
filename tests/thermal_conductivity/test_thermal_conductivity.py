import numpy as np
import xarray as xr
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"

files_hdf5 = [
    "outfile.cumulative_kappa.hdf5",
]


def test_thermal_conductivity(file="outfile.thermal_conductivity"):
    file_ref = folder / file
    file_new = parent / file

    data_ref = np.loadtxt(file_ref)
    data_new = np.loadtxt(file_new)

    np.testing.assert_allclose(data_ref, data_new, err_msg=file_new.absolute())


def test_hdf5(files=files_hdf5):
    for file in files:
        file_ref = folder / file
        file_new = parent / file

        ds_ref = xr.load_dataset(file_ref)
        ds_new = xr.load_dataset(file_new)

        for var in ds_ref.data_vars:
            x = ds_ref[var]
            y = ds_new[var]
            np.testing.assert_allclose(x, y, err_msg=var)


if __name__ == "__main__":
    test_thermal_conductivity()
    test_hdf5()
