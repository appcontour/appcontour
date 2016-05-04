#!/bin/bash
#
# generate kanenobu knot for given p and q
#
# see http://www.ams.org/journals/proc/1986-097-01/S0002-9939-1986-0831406-7/S0002-9939-1986-0831406-7.pdf
#
# entries in knot atlas with same Alexander polynomial are (up to 10 crossings):
#
# 8_8 10_129 10_137
#

p=$1
q=$2

if [ -z "$p" ]
then
  echo "usage: $0 p q"
  exit 1
fi

xp='x'
xq='x'

if [ "$p" -lt 0 ]
then
  xp='X'
  p=$[ - $p ]
fi
if [ "$q" -lt 0 ]
then
  xq='X'
  q=$[ - $q ]
fi

echo "knot {"
echo "        ^       ;"
echo "  |   ^   ^   | ;"

for i in `seq 1 $p`
do
  echo "  |  |  $xp  |  | ;"
done

echo " | ^ || ^ || ^ |;"
echo "  X X X   x x x ;"
echo " | U |X   x| U |;"
echo " |   U | | U   |;"

for i in `seq 1 $q`
do
  echo "  |     $xq     | ;"
done

echo "    U       U   ;"
echo "}"

