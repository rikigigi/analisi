#!/bin/bash
set -e

if [ -z "$SP_DIR" ] || [ -z "$BUILD_DIR" ] || [ -z "$SOURCE_DIR" ] 
then
   echo "please export the env variables SP_DIR, BUILD_DIR, SOURCE_DIR"
   exit -1
fi

CP="cp -v"

if [ -z "$1" ] 
then
   echo "copying the files"
else
   CP="$1"
fi

mkdir -p "$SP_DIR/pyanalisi"
$CP "$BUILD_DIR"/pyanalisi*.so "$SP_DIR/pyanalisi/"
for f in common.py __init__.py
do
$CP "$SOURCE_DIR/pyanalisi/$f" "$SP_DIR/pyanalisi/"
done

