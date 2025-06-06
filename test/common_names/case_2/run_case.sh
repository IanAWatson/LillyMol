#! /bin/bash
if [ -z "$LILLYMOL_HOME" ] || [ -z "$BUILD_DIR" ]
then 
    # undefined BIN_DIR
    echo "System variables LILLYMOL_HOME and BUILD_DIR are required for running the test"
    echo "Please export LILLYMOL_HOME(local path to LillyMol code) "
    echo "Please export BUILD_DIR(the folder name under the bin folder after build)"
    echo "Example: export LILLYMOL_HOME=/home/user/LillyMol"
    echo "Example: export BUILD_DIR=Linux-gcc-7.2.1" 
    exit 1
else
    BIN_DIR=$LILLYMOL_HOME/bin/$BUILD_DIR
fi

command=$BIN_DIR/common_names
case_id="Case 2"
echo "Testing:  $command"

if [ ! -e "$command" ]
then
    echo "Executable is not found"
    exit 1
fi

name1='output.smi'
golden=out/output.smi
stdout='stdout'
stderr='stderr'
diff_tool=../../fileDiff.sh
$command -S output -s 10000 -D + -v in/input1.smi in/input2.smi >${stdout} 2>${stderr}

$diff_tool ${name1} ${golden}
ret=$?
if [ $ret == 1 ]
then
  echo "$case_id : TEST PASS"
else
  echo "$case_id : TEST FAIL"
fi
rm $name1
rm ${stdout}
rm ${stderr}
