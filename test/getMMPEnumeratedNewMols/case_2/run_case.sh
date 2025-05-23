#!/bin/bash

if [ -z "$LILLYMOL_HOME" ] || [ -z "$BUILD_DIR" ]; then 
    echo "System variables LILLYMOL_HOME and BUILD_DIR are required for running the test"
    echo "Please export LILLYMOL_HOME (local path to LillyMol code)"
    echo "Please export BUILD_DIR (the folder name under the bin folder after build)"
    echo "Example: export LILLYMOL_HOME=/home/user/LillyMol"
    echo "Example: export BUILD_DIR=Linux-gcc-7.2.1" 
    exit 1
else
    BIN_DIR="$LILLYMOL_HOME/contrib/python/mmp"
fi

test_command=getMMPEnumeratedNewMols
case=case_2
case_id="Case 2"

test_top="$LILLYMOL_HOME/test"
test_cmd_top="$test_top/$test_command"

diff_tool=../../fileDiff.sh

command="$BIN_DIR/$test_command.py"

if [ ! -x "$command" ]; then
    echo "$command is not executable or does not exist"
    exit 1
fi

in1="$test_cmd_top/$case/in/test_data_05.smi"
in2="$test_cmd_top/$case/in/test_data_06.pairs"
out=test_data_06.csv
gold_out="$test_cmd_top/$case/out/test_data_06.csv"

echo "Testing: $command"

$command -i "$in1" -p "$in2" -o "$out" --frag_left_col FRAG_L --frag_right_col FRAG_R -H -b 0.3 2>>err.txt

$diff_tool "$out" "$gold_out"
ret=$?

if [ $ret -eq 1 ]
then
    echo "$case_id : TEST PASS"
else
    echo "$case_id : TEST FAIL"
fi

rm -f "$out"
rm -f err.txt
