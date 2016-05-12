#!/bin/bash
#
# generate kanenobu knot for given p and q
#
# see http://www.ams.org/journals/proc/1986-097-01/S0002-9939-1986-0831406-7/S0002-9939-1986-0831406-7.pdf
#
# entries in knot atlas with same Alexander polynomial are (up to 10 crossings):
#
# p = q = 0 -> 4_1 # 4_1
# p = 0, q = 1 -> 10_137 (or it's mirror image)
#
# in caso di "half" twists ci sono anche i nodi 8_8(?) 8_9 10_129 10_155
#

p=$1
q=$2
halftwists="$3"

if [ -z "$p" ]
then
  echo "usage: $0 p q [halftwists]"
  echo ""
  echo "$0 1 1"
  echo " is the same as"
  echo "$0 2 2 halftwists"
  exit 1
fi

if [ -z "$halftwists" ]
then
  p=$[ 2 * $p ]
  q=$[ 2 * $q ]
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

