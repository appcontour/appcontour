#!/bin/bash
#
echo "L2a1"
echo "L4a1"
echo "L5a1"
nums=`seq 1 5`
for n in $nums
do
  echo "L6a$n"
done
echo "L6n1"
list="7:2 21:8 55:28 174:113 548:459"
crossings="6"

for pair in $list
do
  crossings=$[ $crossings + 1 ]
  aa=`echo $pair | cut -f1 -d:`
  nn=`echo $pair | cut -f2 -d:`
  for n in `seq 1 $aa`
  do
    echo "L${crossings}a${n}"
  done
  for n in `seq 1 $nn`
  do
    echo "L${crossings}n${n}"
  done
done
