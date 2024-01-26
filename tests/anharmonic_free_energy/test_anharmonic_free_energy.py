import numpy as np
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"


def _read_file(file):
    data = np.loadtxt(file)
    return data


def test_output(file="outfile.anharmonic_free_energy"):
    data_ref = _read_file(folder / file)
    data_new = _read_file(parent / file)

    np.testing.assert_allclose(data_ref, data_new, err_msg=(parent / file).absolute())


if __name__ == "__main__":
    test_output()
