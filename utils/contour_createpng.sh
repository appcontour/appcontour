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

list=`ls *_${ex}*`

for file in $list
do
  num=`echo "$file" | cut -f1 -d_`
  contour printmorse $file | showcontour --ge xfig --xfigspecial --skiprtime 3.0 --title $num
  fig2dev -L png $num.fig >$num.png
  echo contour $num converted...
done
