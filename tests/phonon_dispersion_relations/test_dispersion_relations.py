import numpy as np
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"


def test_dispersion(file="outfile.dispersion_relations"):
    file_ref = folder / file
    file_new = parent / file

    data_ref = np.loadtxt(file_ref)
    data_new = np.loadtxt(file_new)

    np.testing.assert_allclose(data_ref, data_new, err_msg=file_new.absolute())


if __name__ == "__main__":
    test_dispersion()
