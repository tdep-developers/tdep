#! /usr/bin/env sh

cd src/libolle
cat>gitinformation<<EOF
#define gitbranch '${gitbranch}'
#define gitrevision '${gitrevision}'
EOF
cd ../..
