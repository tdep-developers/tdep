#!/usr/bin/env bash
set -ex

# Create the lib directory to store the Fortran library in
LIB_DIR="${PREFIX}/lib"
mkdir -p $LIB_DIR

# Create the bin directory to store the Fortran binaries in
echo "FKDEV: echo BIN_DIR"
echo $BIN_DIR
BIN_DIR="${PREFIX}/bin"
mkdir -p $BIN_DIR
echo "FKDEV: echo BIN_DIR after creation"
echo $BIN_DIR

cd $SRC_DIR

# get important settings
echo "FKDEV: cp settings"
cp .conda/important_settings .

# echo "FKDEV: ls PREFIX lib"
# ls $PREFIX/lib 
# echo "FKDEV: ls PREFIX include"
# ls $PREFIX/include
# echo "FKDEV: ls BUILD_PREFIX/bin"
# ls $BUILD_PREFIX/bin


# let's build this thing
echo "FKDEV: build"
# bash build_things.sh --nthreads_make 2
bash build_things.sh

# Install binaries explicitly

for file in bin/*
do
	echo $file
	cp $file $BIN_DIR
done

# Optional: print installed binaries to verify
echo "Installed binaries:"
ls -l $BIN_DIR
