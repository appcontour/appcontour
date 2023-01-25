#!/bin/bash
#

p=$1
q=$2

if [ -z "$q" ]
then
  echo "Usage: $0 p q"
  exit 999
fi

p=$[ $p ]
q=$[ $q ]

if [ "$p" -lt "2" -o "$q" -lt "2" ]
then
  echo "Both p and q must be larger than 1"
  exit 3
fi

function rawbraid ()
{
  p=$1
  q=$2

  for i in `seq $p`
  do
    moves=$[ $q - 1 ]
    for k in `seq $moves`
    do
      vertsbefore=$[ $moves - $k ]
      vertsafter=$[ $k - 1 ]
      for s in `seq $vertsbefore`
      do
        echo -n "|"
      done
      echo -n "X"
      for s in `seq $vertsafter`
      do
        echo -n "|"
      done
      echo ""
    done
  done
}

#
# q is the number of strands
#
function braidtoknot ()
{
  q=$1
  for n in `seq $q`
  do
    verts=$[ $n - 1 ]
    spaces=$[ $q - $verts - 1 ]

    for m in `seq $verts`
    do
      echo -n "|"
    done
    for m in `seq $spaces`
    do
      echo -n " "
    done
    echo -n "^ "
    for m in `seq $spaces`
    do
      echo -n " "
    done
    for m in `seq $verts`
    do
      echo -n "|"
    done
    echo ";"
  done

  while read braidline
  do
    echo -n "$braidline "
    for n in `seq $q`
    do
      echo -n "|"
    done
    echo ";"
  done

  for n in `seq $q -1 1`
  do
    verts=$[ $n - 1 ]
    spaces=$[ $q - $verts - 1 ]
    for m in `seq $verts`
    do
      echo -n "|"
    done
    for m in `seq $spaces`
    do
      echo -n " "
    done
    echo -n "U "
    for m in `seq $spaces`
    do
      echo -n " "
    done
    for m in `seq $verts`
    do
      echo -n "|"
    done
    echo ";"
  done
}

echo "knot {"
rawbraid $p $q | braidtoknot $q
echo "}"
