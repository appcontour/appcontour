#
# perform a random move among direct and inverse moves
#
appcont=$1
repeat=$2
tmpfile=/tmp/contour_randommove_$$.sketch
tmpfile2=/tmp/contour_randommove_$$_2.sketch

if [ -z "$appcont" ]
then
  echo "Usage: $1 contour [iterations]"
  exit 1
fi

if [ -z "$repeat" ]; then repeat=1; fi

if ! contour ishuffman $appcont -q 2>/dev/null
then
  echo "Contour $appcont is not a labelled apparent contour"
  exit 2
fi

contour print $appcont >$tmpfile 2>/dev/null

for n in `seq $repeat`
do
  list1=`contour rules $tmpfile -q 2>/dev/null`
  list2=`contour listma $tmpfile -q 2>/dev/null`
  list3=`contour listinvl $tmpfile -q 2>/dev/null`
  list4=`contour listinvs $tmpfile -q 2>/dev/null`

  list="$list1 $list2 $list3 $list4"
  list=`echo $list`
  list=`echo "$list" | tr ' ' '\n'`

  count=`echo $list | wc -w`

  rand=$[ $RANDOM % $count ]
  rand=$[ $rand + 1 ]

  rule=`echo "$list" | head -n $rand | tail -n 1`

  echo Randomly selected rule: $rule >&2
  contour rule $rule $tmpfile >$tmpfile2 2>/dev/null
  mv $tmpfile2 $tmpfile
done

cat $tmpfile
rm $tmpfile
