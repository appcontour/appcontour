#!/bin/bash
#
# generate kanenobu knot K(p,q) for given p and q
#
# [ref: Taizo Kanenobu, Exampes on Polynoomial Invariants of Knots and Links,
# Math. Ann. 275, 555-572 (1986)]
#
# see also http://www.ams.org/journals/proc/1986-097-01/S0002-9939-1986-0831406-7/S0002-9939-1986-0831406-7.pdf
# where K_{p,q} is K(-2p,-2q)
#
# entries in knot atlas with same Alexander polynomial are (up to 10 crossings):
#
# p = q = 0 -> 4_1 # 4_1
# p = 0, q = 2 -> 10_137 (or it's mirror image)
#
# in caso di "half" twists ci sono anche i nodi 8_8(?) 8_9 10_129 10_155
#

p=$1
q=$2
completetwists="$3"

if [ -z "$p" ]
then
  echo "usage: $0 p q [completetwists]"
  echo ""
  echo "$0 1 1 completetwists"
  echo " is the same as"
  echo "$0 2 2"
  exit 1
fi

if [ -n "$completetwists" ]
then
#
# should also change sign
#
  p=$[ 2 * $p ]
  q=$[ 2 * $q ]
fi
xp='X'
xq='X'

if [ "$p" -lt 0 ]
then
  xp='x'
  p=$[ - $p ]
fi
if [ "$q" -lt 0 ]
then
  xq='x'
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

