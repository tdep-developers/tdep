import numpy as np
from pathlib import Path

parent = Path(__file__).parent
folder = f"{parent}/reference"


def test_dispersion(file="outfile.many_dynamical_matrices_DDB"):
    file_ref = f"{folder}/{file}"
    file_new = f"{parent}/{file}"

    data_ref = np.loadtxt(file_ref)
    data_new = np.loadtxt(file_new)

    np.testing.assert_allclose(data_ref, data_new, err_msg=file_new.absolute())
    print ("all done and all close")


if __name__ == "__main__":
    test_dispersion()
