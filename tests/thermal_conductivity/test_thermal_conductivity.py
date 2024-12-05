import numpy as np
import xarray as xr
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"

file = "outfile.thermal_conductivity"
file_hdf5 = "outfile.thermal_conductivity_grid.hdf5"


def test_thermal_conductivity(file=file, atol=1, rtol=0.1):
    file_ref = folder / file
    file_new = parent / file

    data_ref = np.loadtxt(file_ref)
    data_new = np.loadtxt(file_new)

    np.testing.assert_allclose(
        data_ref, data_new, atol=atol, rtol=rtol, err_msg=file_new.absolute()
    )


def test_conductivity_comparison(
    file=parent / file,
    file_grid=parent / file_hdf5,
):
    """Check if the thermal conducivities conincide"""
    data1 = np.loadtxt(file)

    # without off-diagonal contribution
    kappa1 = data1[0][:3] + data1[1][:3]

    ds = xr.load_dataset(file_grid)

    kappa2 = np.diag(ds.thermal_conductivity.sum(axis=(0, 1)) / ds.number_of_qpoints)

    np.testing.assert_allclose(kappa1, kappa2)


if __name__ == "__main__":
    test_thermal_conductivity()
    test_conductivity_comparison()
