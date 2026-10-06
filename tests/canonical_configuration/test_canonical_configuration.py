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


def _read_configurations(file, nconf=2):
    """Read the concatenated contcar files of one run, one array per configuration"""
    with open(file) as f:
        lines = f.readlines()

    assert len(lines) % nconf == 0, file
    nlines = len(lines) // nconf
    configurations = []
    for ii in range(nconf):
        chunk = lines[ii * nlines : (ii + 1) * nlines]
        numbers = [float(s) for line in chunk for s in line.split() if _check_numeric(s)]
        configurations.append(np.array(numbers))

    return configurations


def _differ(a, b):
    return a.shape != b.shape or not np.allclose(a, b)


def test_given_same_seed_when_run_twice_then_configurations_are_identical():
    """
    Scenario: The same seed reproduces the configurations
      Given the seed 1
      When canonical_configuration is run twice with 2 configurations each
      Then both runs give identical configurations
    """
    confs_a = _read_configurations(parent / "seed_1_a.dat")
    confs_b = _read_configurations(parent / "seed_1_b.dat")

    for conf_a, conf_b in zip(confs_a, confs_b):
        np.testing.assert_allclose(conf_a, conf_b)


def test_given_seed_when_generating_two_configurations_then_they_differ():
    """
    Scenario: Configurations within one seeded run differ
      Given the seed 1
      When canonical_configuration generates 2 configurations in one run
      Then the two configurations differ
    """
    conf_1, conf_2 = _read_configurations(parent / "seed_1_a.dat")

    assert _differ(conf_1, conf_2), "configurations 1 and 2 are identical"


def test_given_different_seeds_when_run_then_configurations_differ():
    """
    Scenario: Different seeds give different configurations
      Given the seeds 0 and 1
      When canonical_configuration is run once per seed
      Then the two seeds give different configurations
    """
    confs_0 = _read_configurations(parent / "seed_0.dat")
    confs_1 = _read_configurations(parent / "seed_1_a.dat")

    for conf_0, conf_1 in zip(confs_0, confs_1):
        assert _differ(conf_0, conf_1), "seeds 0 and 1 give the same configuration"


def test_given_no_seed_when_run_twice_then_configurations_differ():
    """
    Scenario: Runs without a seed differ
      Given no seed
      When canonical_configuration is run twice, 1 s apart
      Then the two runs give different configurations
    """
    confs_a = _read_configurations(parent / "noseed_a.dat")
    confs_b = _read_configurations(parent / "noseed_b.dat")

    for conf_a, conf_b in zip(confs_a, confs_b):
        assert _differ(conf_a, conf_b), "runs without seed give the same configuration"


if __name__ == "__main__":
    test_log()
    test_given_same_seed_when_run_twice_then_configurations_are_identical()
    test_given_seed_when_generating_two_configurations_then_they_differ()
    test_given_different_seeds_when_run_then_configurations_differ()
    test_given_no_seed_when_run_twice_then_configurations_differ()
