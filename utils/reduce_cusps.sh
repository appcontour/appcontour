#!/bin/bash
#
example="$1"
contour="contour"

if [ -z "$example" ]; then echo "Usage: $0 example"; exit 1; fi
if [ ! -f "$example" ]; then echo "Cannot find $example"; exit 2; fi

#
# check if this is a Huffman contour
#
if ! $contour -q ishuffman $example 2>/dev/null
then
  echo "This is not a Huffman contour"
  exit 3
fi

#
# count number of cusps
#

cuspnum=$( $contour -q info $example 2>/dev/null| grep -i "^cusps: *[0-9]" | cut -f2 -d: )
cuspnum=$( echo $cuspnum )
if [ "$cuspnum" = "0" ]
then
  echo "There are no cusps!" >&2
  echo ""
  exit 0
fi

#
# are there applicable swallowtails?
#
aprules=$( $contour -q testallrules $example 2>/dev/null)
if echo "$aprules" | grep -q -w -i cn1
then
  echo "Can remove a swallowtail" >&2
  echo "CN1"
  exit 0
fi

#
# first look for adjacent cusps
#

regex="^arc *[0-9]*:.*[[(][0-9]* *[0-9]* "
lines=$( $contour print $example 2>/dev/null| grep -i "$regex" )

if [ -n "$lines" ]
then
  firstline=$( echo "$lines" | head -n 1 )
  echo "There are adjacent cusps" >&2
  arc=$( echo $firstline | cut -f1 -d: | cut -f2 -d' ' )
  firstline=$( echo $firstline | cut -f2 -d: | tr "[" "(" | tr "])" " " )
  dvals=$( echo $firstline | cut -f2 -d'(' | tr -s ' ' )
  dvals=$( echo $dvals )
  n1=$( echo $dvals | cut -f1 -d' ' )
  n2=$( echo $dvals | cut -f2 -d' ' )
  n3=$( echo $dvals | cut -f3 -d' ' )
  delta1="+"
  delta2="+"
  if [ "$n2" -lt "$n1" ]; then delta1="-"; fi
  if [ "$n3" -lt "$n2" ]; then delta2="-"; fi
  if [ "$delta1" != "$delta2" ]
  then
    echo "C2"
    exit 0
  fi
  s1="0";s2="2"
  if [ "$delta1" = "-" ]; then s1="2"; s2="0"; fi
  regex2="[-]a *${arc}:$s1 *[-]a *${arc}:$s2 "
  lines=$( $contour listmergearcs $example 2>/dev/null | grep -i invn4 | grep -i "$regex2" )
  firstline=$( echo "$lines" | head -n 1 )
  # echo "arc $arc, deltas=$delta1 $delta2"
  rule=$( echo $firstline | cut -f1 -d')' | cut -f2 -d'(' )
  echo "$rule CN1"
  exit 0
fi

echo "Cannot find adjacent cusps"
exit 4
