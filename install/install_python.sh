#!/bin/bash
set -e

if [ -z "$SP_DIR" ] || [ -z "$BUILD_DIR" ] || [ -z "$SOURCE_DIR" ] 
then
   echo "please export the env variables SP_DIR, BUILD_DIR, SOURCE_DIR"
   echo "you can set SP_DIR=__DETECT__ to use the output of the command"
   echo "    python -c 'import sysconfig; print(sysconfig.get_paths()[\"purelib\"])'"
   echo "as SP_DIR. Usually this is the path where you want to install stuff."
   exit -1
fi

if [ "$SP_DIR" = "__DETECT__" ]
then
   export SP_DIR="$(python -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
   echo "Using detected SP_DIR='$SP_DIR'"
fi

CP="cp -v"

if [ -z "$1" ] 
then
   echo "copying the files"
else
   echo "using the command '$1' to install the files"
   CP="$1"
fi

mkdir -p "$SP_DIR/pyanalisi"
$CP "$BUILD_DIR"/pyanalisi*.so "$SP_DIR/pyanalisi/"
for f in common.py __init__.py
do
$CP "$SOURCE_DIR/pyanalisi/$f" "$SP_DIR/pyanalisi/"
done

