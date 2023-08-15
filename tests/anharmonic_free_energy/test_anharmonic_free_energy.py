import numpy as np
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"


def _read_file(file):
    with open(file) as f:
        next(f)
        next(f)
        F = float(next(f).split()[2])
        next(f)
        next(f)
        U0 = float(next(f).split()[2])
        Fp = float(next(f).split()[2])

    return np.array([F, U0, Fp])


def test_output(file="outfile.anharmonic_free_energy"):
    data_ref = _read_file(folder / file)
    data_new = _read_file(parent / file)

    assert np.allclose(data_ref, data_new), (parent / file).absolute()


if __name__ == "__main__":
    test_output()
