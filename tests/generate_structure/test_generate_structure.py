import numpy as np
import xarray as xr
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"


def _read_log(file):
    """Return number of atoms only"""
    natoms = int(open(file).readlines()[9].split()[3])
    return natoms

def test_log(file="generate_structure.log"):
    data_ref = _read_log(folder / file)
    data_new = _read_log(parent / file)

    assert np.allclose(data_ref, data_new), (parent / file).absolute()


if __name__ == "__main__":
    test_log()
