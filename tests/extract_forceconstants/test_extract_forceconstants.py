import numpy as np
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"


def _read_fc(file, offset1, offset2):
    rows = []
    with open(file) as f:
        # read header
        next(f)
        next(f)
        nneighbors = int(next(f).split()[0])

        lines = f.readlines()

    for ii in range(nneighbors):
        row = lines[ii * offset1 + offset2 : (ii + 1) * offset1]
        rows.append(np.fromstring("".join(row), sep=" "))

    return np.array(rows)


def _read_fc2(file):
    return _read_fc(file, offset1=5, offset2=2)


def _read_fc3(file):
    return _read_fc(file, offset1=15, offset2=3)


def test_fc2(file="outfile.forceconstant"):
    data_ref = _read_fc2(folder / file)
    data_new = _read_fc2(parent / file)

    np.testing.assert_allclose(data_ref, data_new, err_msg=(parent / file).absolute())


def test_fc3(file="outfile.forceconstant_thirdorder"):
    data_ref = _read_fc3(folder / file)
    data_new = _read_fc3(parent / file)

    np.testing.assert_allclose(data_ref, data_new, err_msg=(parent / file).absolute())


if __name__ == "__main__":
    test_fc2()
    test_fc3()
