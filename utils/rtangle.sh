#!/bin/bash
#

fraction=$1

num=$( echo $fraction | cut -f1 -d/ )
den=$( echo $fraction | cut -s -f2 -d/ )

function usage ()
{
  echo "usage: $0 fraction"
  echo "  where fraction is of the form p/q or -p/q"
}

isdecimal() {
  # filter octal/hex/ord()
  numb=$(printf '%s' "$1" | tr -d '-' | sed "s/^0*\([1-9]\)/\1/; s/'/^/")

  test "$numb" && printf '%f' "$numb" >/dev/null 2>&1
}

if [ -z "$num" -o -z "$den" ]
then
  usage
  exit 2
fi

if ! isdecimal $num
then
  usage
  exit 3
fi

if ! isdecimal $den
then
  usage
  exit 4
fi

cross="X"
if [ "$num" -lt 0 ]
then
  cross="x"
  num=$[ 0 - $num ]
fi

if [ $den -le 0 ]
then
  echo "denominator must be a positive number"
  exit 5
fi

if [ $num -eq 0 ]
then
  echo "numerator must be nonzero"
  exit 6
fi

echo "knot {"
echo " ^;"

numm1=$[ $num - 1 ]
echo -n "|"
for j in $( seq $numm1 )
do
  echo -n " ^"
done
echo " |;"

for i in $( seq $den )
do
  if [ $i -ne 1 ]
  then
    echo -n "|"
    for j in $( seq $numm1 )
    do
      echo -n " $cross"
    done
    echo " |;"
  fi

  for j in $( seq $num )
  do
    echo -n " $cross"
  done
  echo " ;"

done

echo -n "|"
for j in $( seq $numm1 )
do
  echo -n " U"
done
echo " |;"

echo " U;"
echo "}"
