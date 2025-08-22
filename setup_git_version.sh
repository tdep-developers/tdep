#! /usr/bin/env sh

# It is to be used with the MEson build system. There is no need to run it when
# using the ./build_things.sh script.

cd src/libolle
cat>gitinformation<<EOF
#define gitbranch '${gitbranch}'
#define gitrevision '${gitrevision}'
EOF
cd ../..
