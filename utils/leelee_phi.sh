#!/bin/bash
#
# generate the \Phi_n handlebody knot described in
#
# Jung Hoon Lee, Sangyop Lee, Inequivalent handlebody-knots with homeomorphic complements
# Algebraic and Geometric Topology 12 (2012) 1059-1079
#

p=$1

if [ -z "$p" ]
then
  echo "usage: $0 p"
  echo ""
  echo "where p is the number of complete twists"
  exit 1
fi

xp='X'

if [ "$p" -lt 0 ]
then
  xp='x'
  p=$[ - $p ]
fi

echo "knot {"
echo "        ^      ;"
echo "  |    ^     | ;"

for i in `seq 1 $p`
do
  echo "   $xp       | | ;"
  echo "   $xp       | | ;"
done

echo " | ^ |     | | ;"
echo "  x X      | | ;"
echo "  x X      | | ;"
echo " | y |     | | ;"
echo " |  X      | | ;"
echo "  y     U    | ;"
echo "        U      ;"
echo "}"

