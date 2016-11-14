#!/bin/bash
#
# https://oeis.org/A002863
# https://oeis.org/A002864
#
#echo "K0a1"
list="0:0 0:0 1:1 1:1 2:2 3:3 7:7 21:18 49:41 165:123 552:367 2176:1288 9988:4878 46972:19536"
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
for pair in $list
do
  crossings=$[ $crossings + 1 ]
  if [ "$crossings" -gt "$maxcrossings" ]; then break; fi
  total=`echo "$pair" | cut -f1 -d':'`
  alternating=`echo "$pair" | cut -f2 -d':'`
  nonalternating=$[ $total - $alternating ]
  for n in `seq 1 $alternating`
  do
    echo "K${crossings}a${n}"
  done
  for n in `seq 1 $nonalternating`
  do
    echo "K${crossings}n${n}"
  done
done

