import numpy as np
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"


def test_dispersion(file="outfile.many_dynamical_matrices_DDB"):
    file_ref = folder / file
    file_new = parent / file

    data_ref = np.loadtxt(file_ref, skiprows=151, comments=["q", "2n"])
    data_new = np.loadtxt(file_new, skiprows=151, comments=["q", "2n"])

    np.testing.assert_allclose(data_ref, data_new, err_msg=file_new.absolute())
    print("all done and all close")


if __name__ == "__main__":
    test_dispersion()
