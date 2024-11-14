Test suite
===

## Create test files for all binaries

```bash
make testfiles
```

will run each TDEP binary from the most recent TDEP build in `../bin` in the specific folder, and create the corresponding output. Your `PATH` is modified by `00-set_path.sh` to make sure that you are not accidentally running binaries from a previous TDEP build.

## Check all outputs

To check the files that were created against reference data, please run

```bash
make test
```

This will run python scripts in each folder to check the produced output. **Please note that this requires some python dependencies which you find in `requirements.txt` and install via `pip install -r requirements.txt`**.

Run

```bash
make all
```

to run everything.

## Clean up

Use

```bash
make clean
```

to clean test files.

## Reference files

Were created with TDEP commit `8c01e0343e4098f1d160efc141f8af6ae7f54941` on a macbook with M1 chip.

## Known issues

- Errors like
  ```
  atomic_distribution/test_atomic_distribution.py::test_hdf5
    <frozen importlib._bootstrap>:241: RuntimeWarning: numpy.ndarray size changed, may indicate binary incompatibility. Expected 16 from C header, got 96 from PyObject
  ```

  hint towards compatibility issues in your python environment. **They are unrelated to TDEP.**

## Optional tests

There are optinal tests in `./optional`. The test mechanism to test these optional binaries is the same.
