#!/bin/bash
#
# From "A table of genus two handlebody-knots up to six crossings"
# by Ishii-Kishimoto-Moriuchi
#
list="0 0 0 1 4 16"
maxcrossings=`echo "$list" | wc -w`

if [ -n "$1" ]
then
  if [ "$1" -le "$maxcrossings" ]
  then
    maxcrossings=$1
  else
    echo "usage: $0 [maxcrossings]"
    echo "  maxcrossings: integer not larger then $maxcrossings"
    exit 1
  fi
fi

echo "HK0_1"
crossings=0
for total in $list
do
  crossings=$[ $crossings + 1 ]
  if [ "$crossings" -gt "$maxcrossings" ]; then break; fi
  for n in `seq 1 $total`
  do
    echo "HK${crossings}_${n}"
  done
done

