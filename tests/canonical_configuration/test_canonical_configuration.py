import numpy as np
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"


def _check_numeric(character):
    try:
        float(character)
        return True
    except ValueError:
        return False


def _read_log(file):
    rows = []
    with open(file) as f:
        for line in f:
            rows.append([s for s in line.split() if _check_numeric(s)])

    return np.nan_to_num(np.array(rows, dtype=float))


def test_log(file="canonical_configuration.dat"):
    data_ref = _read_log(folder / file)
    data_new = _read_log(parent / file)

    np.testing.assert_allclose(data_ref, data_new, err_msg=(parent / file).absolute())


if __name__ == "__main__":
    test_log()
