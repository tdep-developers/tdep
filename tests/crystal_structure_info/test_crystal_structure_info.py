from pathlib import Path
import filecmp

parent = Path(__file__).parent
folder = parent / "reference"


def test_output(file="outfile.qpoints_dispersion"):
    file_ref = folder / file
    file_new = parent / file

    assert filecmp.cmp(file_ref, file_new), file_new.absolute()


if __name__ == "__main__":
    test_output()
