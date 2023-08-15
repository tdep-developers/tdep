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

    assert np.allclose(data_ref, data_new), f"CHECK {file_new.absolut()}"


def test_hdf5(files=files_hdf5):
    for file in files:
        file_ref = folder / file
        file_new = parent / file

        ds_ref = xr.load_dataset(file_ref)
        ds_new = xr.load_dataset(file_new)

        for var in ds_ref.data_vars:
            x = ds_ref[var]
            y = ds_new[var]
            assert np.allclose(x, y), file_new


if __name__ == "__main__":
    test_mean_square_displacement()
    test_hdf5()
