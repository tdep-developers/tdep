import numpy as np
import xarray as xr
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"

files_hdf5 = [
    "outfile.grid_thermal_conductivity_4ph.hdf5",
]


def test_thermal_conductivity(file="outfile.thermal_conductivity_4ph", atol=20, rtol=5):
    file_ref = folder / file
    file_new = parent / file

    data_ref = np.loadtxt(file_ref)
    data_new = np.loadtxt(file_new)

    np.testing.assert_allclose(
        data_ref, data_new, atol=atol, rtol=rtol, err_msg=file_new.absolute()
    )


def test_hdf5(files=files_hdf5, atol=1, rtol=0.01):
    for file in files:
        file_ref = folder / file
        file_new = parent / file

        ds_ref = xr.load_dataset(file_ref)
        ds_new = xr.load_dataset(file_new)

        for var in ds_ref.data_vars:
            x = ds_ref[var]
            y = ds_new[var]
            np.testing.assert_allclose(x, y, atol=atol, rtol=rtol, err_msg=var)


if __name__ == "__main__":
    test_thermal_conductivity()
    test_hdf5()
