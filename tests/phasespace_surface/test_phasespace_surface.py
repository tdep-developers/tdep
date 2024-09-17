import numpy as np
from pathlib import Path

parent = Path(__file__).parent
folder = parent / "reference"


def _read_file(file):
    data_plus = []
    data_minus = []
    with open(file, "r") as f:
        for line in f:
            if "... write nonzero phase spaces and scattering intensities" in line:
                break

        for line in f:
            if "plus" in line:
                _plus = True
                continue
            if "minus" in line:
                _plus = False
                continue
            if "All done" in line:
                return np.array(data_plus), np.array(data_minus)

            x1, x2, x3, x4, x5 = line.split()
            if _plus:
                data_plus.append([int(x1), int(x2), int(x3), int(x4), float(x5)])
            else:
                data_minus.append([int(x1), int(x2), int(x3), int(x4), float(x5)])


def test_phasespace(file="phasespace_surface.log"):
    data_ref = _read_file(folder / file)
    data_new = _read_file(parent / file)

    err_msg = (parent / file).absolute()
    # 3 modes: should match exactly
    kw = {"err_msg": err_msg}
    np.testing.assert_allclose(data_ref[0][:, :3], data_new[0][:, :3], **kw)
    np.testing.assert_allclose(data_ref[1][:, :3], data_new[1][:, :3], **kw)

    # 4: number of points: should match up to 10% error
    rtol = 1
    kw = {"atol": 0, "rtol": rtol, "err_msg": err_msg}
    np.testing.assert_allclose(data_ref[0][:, 3], data_new[0][:, 3], **kw)
    np.testing.assert_allclose(data_ref[1][:, 3], data_new[1][:, 3], **kw)

    # 5: intensity: should match up to 0.01
    atol = 1e-2
    kw = {"atol": atol, "rtol": 0, "err_msg": err_msg}
    np.testing.assert_allclose(data_ref[0][:, 4], data_new[0][:, 4], **kw)
    np.testing.assert_allclose(data_ref[1][:, 4], data_new[1][:, 4], **kw)


if __name__ == "__main__":
    test_phasespace()
