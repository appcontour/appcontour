#!/bin/bash

ex=$1

if [ -z "$ex" ]
then
  echo "usage: $0 example"
  exit 1
fi

if [ ! -d "${ex}.transformations" ]
then
  echo "cannot find transformations directory"
  exit 1
fi
cd ${ex}.transformations

echo 'I will move all duplicated contours into the "duplicated"'
echo "subdirectory of ${ex}.transformations"

mkdir duplicated

count=`ls *_${ex}* | cut -f1 -d'_' | tail -n 1`
count=${count#0}
count=${count#0}
count=${count#0}
count=${count#0}
count=${count#0}

nums=`seq $count`

for nn in $nums
do
  n=`printf "%05d" $nn`
  list=`ls ${n}* | tail -n +2`
  for f in $list
  do
    mv $f duplicated
  done
done
