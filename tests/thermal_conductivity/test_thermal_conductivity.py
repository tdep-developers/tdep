import numpy as np
import xarray as xr
from pathlib import Path
import pytest

parent = Path(__file__).parent
folder = parent / "reference"

files_hdf5 = [
    "outfile.thermal_conductivity_grid.hdf5",
]


def test_thermal_conductivity(file="outfile.thermal_conductivity", atol=20, rtol=5):
    file_ref = folder / file
    file_new = parent / file

    data_ref = np.loadtxt(file_ref)
    data_new = np.loadtxt(file_new)

    np.testing.assert_allclose(
        data_ref, data_new, atol=atol, rtol=rtol, err_msg=file_new.absolute()
    )


def get_variables(file):
    """Helper function to get variables from a dataset"""
    return list(xr.load_dataset(folder / file).data_vars)


def generate_test_cases(files=files_hdf5):
    """Generate test cases combining files and their variables"""
    test_cases = []
    for file in files:
        vars = get_variables(file)
        for var in vars:
            test_cases.append((file, var))
    return test_cases


@pytest.mark.parametrize("file,var", generate_test_cases())
def test_hdf5(file, var, atol=1, rtol=0.01):
    """Test HDF5 files by comparing variables between reference and new datasets"""
    file_ref = folder / file
    file_new = parent / file

    ds_ref = xr.load_dataset(file_ref)
    ds_new = xr.load_dataset(file_new)

    x = ds_ref[var]
    y = ds_new[var]

    np.testing.assert_allclose(x, y, atol=atol, rtol=rtol, err_msg=var)


if __name__ == "__main__":
    test_thermal_conductivity()
    test_hdf5()
