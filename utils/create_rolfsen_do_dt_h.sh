#!/bin/bash
#
# use file "Rolfsen_DT.html" from knotscape package as input file
#
function display ()
{
  echo "  \"$1\", \"$2\","
}

crossings=""
lastknots=""
while read a b c
do
  #echo "$a-$b-$c"
  if [ "$c" = "crossings" ]
  then
    a=$b
    b=$c
  fi
  if [ "$b" = "crossings" ]
  then
    crossings="$a"
    alternating="a"
    if [ "$c" = "nonalternating" ]
    then
      alternating="n"
    fi
  fi
  if echo "$b" | grep -E -q "^[0-9]+$"
  then
    if echo "$a" | grep -q =
    then
      lastknots="1"
      a1=`echo "$a" | cut -f1 -d=`
      a=`echo "$a" | cut -f2 -d=`
      display ${crossings}_$a1 ${crossings}${alternating}_$b
    fi
    if [ -n "$lastknots" ]
    then
      display ${crossings}_$a "#${crossings}${alternating}_$b"
      aprev=$[ $a - 1 ]
      display ${crossings}_${a}R "${crossings}${alternating}_$b"
      display ${crossings}_${aprev}KA "${crossings}${alternating}_$b"
    else
      display ${crossings}_$a ${crossings}${alternating}_$b
    fi
    #echo "{\"${crossings}_$a\", \"${crossings}${alternating}_$b\"},"
  fi
done

echo "  0"
