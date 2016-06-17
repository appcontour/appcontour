#!/bin/bash
#
# the resulting file, named links_gausscodes.txt must be placed
# in /usr/local/share/data/
#
file="Links.rdf"

if [ -n "$1" ]
then
  file=$1
fi

if [ ! -f "$file" ]
then
  echo "Fatal: cannot find file $file"
  echo "usage: $0 [file] >links_gausscodes.txt"
  exit 1
fi

function parse ()
{
  while read f1 f2 f3
  do
    knot=`echo "$f1" | cut -f2 -d: | cut -f1 -d'>'`
    gc=`echo "$f3" | cut -f2 -d'"' | tr -d ','`
    echo "$knot:$gc"
  done
}

grep "<invariant:Gauss_Code>" $file | parse

