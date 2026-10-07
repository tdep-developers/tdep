import re

import numpy as np
from pathlib import Path

parent = Path(__file__).parent


def _exitcode(case):
    return int((parent / case / "outfile.exitcode").read_text())


def _n_internal_dof(case):
    log = (parent / case / "outfile.log").read_text()
    match = re.search(r"No\. of internal degrees of freedom:\s*(\d+)", log)
    assert match, f"no DOF count in {parent / case / 'outfile.log'}"
    return int(match.group(1))


def _read_poscar(file):
    lines = Path(file).read_text().splitlines()
    lattice = np.array([line.split() for line in lines[2:5]], dtype=float)
    lattice *= float(lines[1])
    species = lines[5].split()
    counts = [int(n) for n in lines[6].split()]
    symbols = [s for s, n in zip(species, counts) for _ in range(n)]
    positions = np.array([line.split()[:3] for line in lines[8 : 8 + sum(counts)]], dtype=float)
    return lattice, symbols, positions


def test_average_is_symmetry_projected():
    """
    Scenario: the average structure is the symmetry-projected average of the trajectory
      Given a tetragonal AB crystal (P4mm, A = Ga, B = N) with A at (0,0,0) and B at (1/2,1/2,0.600)
      And an MD trajectory of two frames:
        frame 1 = B at z = 0.580
        frame 2 = B at z = 0.640, and A additionally shifted by (+0.010, 0, 0)
      When the average structure is computed
      Then 1 internal degree of freedom is reported
      And outfile.ucposcar has A at (0,0,-0.005) and B at (1/2,1/2,0.605)
      And the lattice is unchanged
    """
    case = "tetragonal_ab_average"
    assert _exitcode(case) == 0
    assert _n_internal_dof(case) == 1

    lattice_ref, _, _ = _read_poscar(parent / case / "infile.ucposcar")
    lattice, _, positions = _read_poscar(parent / case / "outfile.ucposcar")
    np.testing.assert_allclose(lattice, lattice_ref, atol=1e-8)

    positions_expected = np.array([[0.0, 0.0, -0.005], [0.5, 0.5, 0.605]])
    # compare modulo lattice translations
    diff = positions - positions_expected
    diff -= np.round(diff)
    np.testing.assert_allclose(diff, 0.0, atol=1e-8)


if __name__ == "__main__":
    test_average_is_symmetry_projected()
