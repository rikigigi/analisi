#!/bin/bash
set -e
ANALISI="$BUILD_DIR/analisi"
RES_DIR="$SOURCE_DIR/tests/data/cli"
INPUT="$SOURCE_DIR/tests/data/lammps2020.bin"
INPUT_L="$SOURCE_DIR/tests/data/gk_integral.dat"


if [ -e "$RES_DIR" ]
then
echo RES_DIR $RES_DIR ok
else
echo RES_DIR CREATING $RES_DIR
mkdir -p "$RES_DIR"
fi

TESTS=( "MSD_normal_full"   "-i $INPUT -Q"
        "MSD_normal"        "-i $INPUT -Q -s 10 -S 50" 
        "MSD_cm"            "-i $INPUT -q -s 10 -S 50"
        "MSD_cm_reference"  "-i $INPUT -Q --mean-square-displacement-self"
        "vibrational"       "-i $INPUT -V"
        "vel_histogram"     "-i $INPUT -v 100 -M 6.0"
        "vel_histogram_def" "-i $INPUT -v 1000"
        "GK_heat_single"    "-l $INPUT_L -H -a c_flux[1]"
        "GK_heat_double"    "-l $INPUT_L -H -a c_flux[1] c_vcm[1][1]"
        "pair_corr_no_t"    "-i $INPUT -g 100 -F 0.0 4.0 -S 1 -s 8"
        "pair_corr_t"       "-i $INPUT -g 100 -F 0.0 4.0 -S 10 -s 8"
        "neighbours"        "-i $INPUT --neighbour 10"
        "spherical_harm"    "-i $INPUT -Y 1 -S 2 -s 14 -F 0.0 2.0"
)
N_arr=${#TESTS[@]}
N=$(( N_arr - 1 ))
RETURN=0
FLIST=""
for i in $(seq 0 2 $N )
do
    ip=$(( i + 1 ))
    NAME=${TESTS[i]}
    ARGS=${TESTS[ip]}
    OUTPUT="`$ANALISI $ARGS 2> /dev/null`"
    FTEST="$RES_DIR/$NAME"
    if [ -f "$FTEST" ]
    then
        OUTPUT_F="`cat $FTEST`"
        DIFF="`diff <(echo "$OUTPUT") <(echo "$OUTPUT_F")`" 
        if [ -z "$DIFF" ]
        then
            echo SUCCESS $NAME
        else
            echo FAILED $NAME
            echo "$DIFF"
            FLIST="$FLIST
FAILED $NAME"
            RETURN=$(( RETURN + 1 ))
        fi
    else
        echo NOT FOUND: creating $FTEST
        echo "$OUTPUT" > "$FTEST"
    fi

done
echo "$FLIST"
exit $RETURN
