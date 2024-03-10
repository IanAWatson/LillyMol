#! /bin/bash

in_name=$1
out_name=$2

echo "Files ${in_name} ${out_name}" >&2
if [[ -z "${out_name}" ]] ; then
  echo "One or both files missing ${in_name} ${out_name}" >&2
  exit 2
fi

ret=0
#echo name1
#echo name1_out

result=$(diff -y -W 72 $in_name $out_name)

if [ $? -eq 0 ]
then
    #return 1 for file match
    ret=1
else
    s1=$(stat -c %s ${in_name})
    s2=$(stat -c %s ${out_name})
    echo "Sizes ${s1} and ${s2}" >&2
    line_count_1=$(wc -l < $in_name)
    #echo $line_count_1
    line_count_2=$(wc -l < $out_name)
    #echo $line_count_2
    if [ $line_count_1 = $line_count_2 ]
    then
        #return 2 for line count match
        #echo "line count match"
        ret=2
    fi
fi
#echo $ret
exit $ret
