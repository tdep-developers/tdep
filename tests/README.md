Test suite
===

## Run binaries and create test files

```bash
make testfiles
```

will run the TDEP binaries **from your current PATH** and create their output.

## Run tests

To check the files that were created against reference data, please run

```bash
make test
```

This will run python scripts in each folder to check the produced output. **Please note that this requires 3 python dependencies which you find in `requirements.txt`**

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