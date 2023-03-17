#!/bin/bash
#
# From "A table of n-component handlebody links of genus n plus 1 up to six crossings"
# by Bellettini-Paolini-Paolini-Wang
#
list="0 0 0 1 1 15"
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

crossings=0
for total in $list
do
  crossings=$[ $crossings + 1 ]
  if [ "$crossings" -gt "$maxcrossings" ]; then break; fi
  for n in `seq 1 $total`
  do
    echo "HL${crossings}_${n}"
  done
done

