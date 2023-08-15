Test suite
===

## Run binaries and create test files

```bash
make testfiles
```

will run the TDEP binaries **from your current PATH** and create their output.

## Run tests

Each folder has a `test_binary_name.py` file, either run this directly with `python test_binary_name.py` or use `pytest` to run them. This also works from the root folder. `pytest` has to be installed in order to use it.

## Clean up

Use

```bash
make clean
```

to clean test files.

## Reference files

Were created with TDEP commit `8c01e0343e4098f1d160efc141f8af6ae7f54941` on a macbook with M1 chip.