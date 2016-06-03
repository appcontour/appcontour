#!/bin/bash

function display ()
{
  echo "  \"$1\", \"$2\","
}

crossings=""
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
      a1=`echo "$a" | cut -f1 -d=`
      a=`echo "$a" | cut -f2 -d=`
      display ${crossings}_$a1 ${crossings}${alternating}_$b
    fi
    display ${crossings}_$a ${crossings}${alternating}_$b
    #echo "{\"${crossings}_$a\", \"${crossings}${alternating}_$b\"},"
  fi
done

echo "  0"
