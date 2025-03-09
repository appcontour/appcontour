#!/bin/bash
#

getgenlist=""
onlyinvariant=""
onlyhash=""
full="yes"

while [ "${1:0:1}" = "-" -a -n "${1:1:1}" ]
do
  case $1 in
    -l)
      getgenlist="yes"
      shift
      ;;
    -i)
      onlyinvariant="yes"
      full=""
      shift
      ;;
    -h | -H | --hash)
      onlyhash="yes"
      full=""
      shift
      ;;
    --ccemb)
      ccemb=$2
      echo "Required ccemb is $ccemb" >&2
      shift 2
      ;;
    *)
      echo "Invalid option: $1"
      exit 2
      ;;
  esac
done

fpgroup=$1
tmpinput=""
if [ "$fpgroup" = "-" ]
then
  tmpinput="/tmp/g_image_input.$$"
  cat >$tmpinput
  fpgroup=$tmpinput
fi

group=$2

if [ -z "$2" ]
then
  echo "usage: $0 [-l][-i][-H|--hash] some_wirtinger_presentation.fpgroup group"
  echo "where group can be e.g. A5"
  echo "use '-' as filename if input comes from stdin"
  echo "an embedding can be indicated in place of a wirtinger presentation"
  echo ""
  echo "  -i      Output canonical text (sorted and without id of homomorphism)"
  echo "  --hash  Output the result of md5sum on the canonical text"
  exit 2
fi

igroup="$group"
if [ -n "$3" ]
then
  igroup=$3
  igroup=`echo "$igroup" | tr -d ' '`
  echo "Selecting homomorphisms with image $igroup"
fi

#
#
# option -x 10000 tells the window size. A large value is set to prevent
# from splitting output lines! (2024/05/06)
#

tmpfile=/tmp/g_image.$$.homos
gap1="gap -r -b -q -x 10000"
gap2="tr -d '\r'"
gap3="grep -v '^#I'"

function feedgap ()
{
#  echo -n "Testing: $which"
  genimage=`echo "Order(Group($genlist));" | $gap1 | $gap2 | $gap3`

#  printf "\r"

#  if [ "$genimage" = "$gorder" ]
#  then
    subimage=`echo "Order(Group($mlist,$llist));" | $gap1 | $gap2 | $gap3`
#    if [ "$subimage" -lt "$genimage" ]
#    then
      #normalcl=`echo "Order(NormalClosure( Group($mlist,$llist), Group($mlist)));" | $gap1 | $gap2 | $gap3`
      gapcmd1="nc:=NormalClosure( Group($mlist,$llist), Group($mlist))"
      gapcmd2="Order(nc)"
      gapcmd3="StructureDescription(nc)"
      gapoutput=`echo "$gapcmd1;$gapcmd2;$gapcmd3;" | $gap1 | $gap2 | $gap3`
      normalclosureorder=`echo "$gapoutput" | head -n 2 | tail -n 1`
      normalclosure=`echo "$gapoutput" | head -n 3 | tail -n 1`
      echo -n "Proper homomorphism #$which: surface group: $subimage - normal closure: $normalclosureorder"
      echo " ($normalclosure)"
#    fi
#  fi
}

function parseone ()
{
  genlist=""
  while read line
  do
    w1=`echo "$line" | cut -f1 -d' '`
    w2=`echo "$line" | cut -f2 -d' '`
    w3=`echo "$line" | cut -f3 -d' '`
    w5=`echo "$line" | cut -f5- -d' '`
    sel=`echo $w5 | tr ' ' ','`
    if [ "$w2" = "End" ]
    then
      mlist=`echo ${mvec[*]} | tr ' ' ','`
      llist=`echo ${lvec[*]} | tr ' ' ','`
      genlist=`echo "$genlist" | cut -c2-`
      if echo "$ontolist" | grep -q "^$which$"
      then
        if echo "$properlist" | grep -q "^$which$"
        then
          feedgap
        fi
      fi
      break
    fi
    if [ "$w1" = "Selected" ]
    then
      selectedpos=${w3:1:1}
      selectedpos=$[ $selectedpos - 1 ]  # number from 0
      pairpos=$[ $selectedpos / 2 ]
      pairwhich=$[ $selectedpos - 2 * $pairpos ]
      if [ $pairwhich = "0" ]; then mvec[$pairpos]=$sel; fi
      if [ $pairwhich = "1" ]; then lvec[$pairpos]=$sel; fi
    else
      line=`echo $line | tr ' ' ','`
      genlist="$genlist,$line"
    fi
  done
}

function parselist ()
{
  while read dum hom which rest
  do
    if [ "$hom" = "Homomorphism" ]
    then
      which=`echo "$which" | cut -c2-`
      # echo "Found homomorphism number $which"
      parseone
    fi
  done
}

function homoimageone ()
{
  genlist=""
  while read line
  do
    w1=`echo "$line" | cut -f1 -d' '`
    w2=`echo "$line" | cut -f2 -d' '`
    w3=`echo "$line" | cut -f3 -d' '`
    w5=`echo "$line" | cut -f5- -d' '`
    sel=`echo $w5 | tr ' ' ','`
    if [ "$w2" = "End" ]
    then
      echo "$genlist" | cut -c2-
      break
    fi
    if [ "$w1" = "Selected" ]
    then
      continue
    fi
    line=`echo $line | tr ' ' ','`
    genlist="$genlist,$line"
  done
}

function homoimagelist ()
{
  while read dum hom which rest
  do
    if [ "$hom" = "Homomorphism" ]
    then
      which=`echo "$which" | cut -c2-`
      # echo "Found homomorphism number $which"
      homoimageone
    fi
  done
}

function properparseone ()
{
  while read line
  do
    w1=`echo "$line" | cut -f1 -d' '`
    w2=`echo "$line" | cut -f2 -d' '`
    w3=`echo "$line" | cut -f3 -d' '`
    w5=`echo "$line" | cut -f5- -d' '`
    sel=`echo $w5 | tr ' ' ','`
    if [ "$w2" = "End" ]
    then
      mlist=`echo ${mvec[*]} | tr ' ' ','`
      llist=`echo ${lvec[*]} | tr ' ' ','`
      echo "Order(Group($mlist,$llist));"
      break
    fi
    if [ "$w1" = "Selected" ]
    then
      selectedpos=${w3:1:1}
      selectedpos=$[ $selectedpos - 1 ]  # number from 0
      pairpos=$[ $selectedpos / 2 ]
      pairwhich=$[ $selectedpos - 2 * $pairpos ]
      if [ $pairwhich = "0" ]; then mvec[$pairpos]=$sel; fi
      if [ $pairwhich = "1" ]; then lvec[$pairpos]=$sel; fi
    fi
  done
}

function properparselist ()
{
  while read dum hom which rest
  do
    if [ "$hom" = "Homomorphism" ]
    then
      which=`echo "$which" | cut -c2-`
      # echo "Found homomorphism number $which"
      properparseone
    fi
  done
}

function ordergroup ()
{
  while read line
  do
    echo "Order(Group($line));"
  done
}

function describegroup ()
{
  while read line
  do
    echo "StructureDescription(Group($line));"
  done
}

case $group in
  S4)
    gorder=24
    ;;
  A4)
    gorder=12
    ;;
  S5)
    gorder=120
    ;;
  A5)
    gorder=60
    ;;
  S6)
    gorder=720
    ;;
  A6)
    gorder=360
    ;;
  S7)
    gorder=5040
    ;;
  A7)
    gorder=2520
    ;;
  *)
    echo "Cannot deal with group $group"
    exit 3
    ;;
esac

if [ -z "$onlyhash" ]
then
  echo "Numbers correspond to the order of the group"
fi

echo -n "Computing homomorphisms (in $tmpfile)..." >&2

#if grep -q "^embedding" $fpgroup
ccembdef=1
if [ -n "$ccemb" ]; then ccembdef=$ccemb; fi
#
# I don't know the reason for option '-ccemb $ccembdef' below
# possibly we can remove it
#
if contour wirtinger $fpgroup --ccemb $ccembdef >/dev/null 2>/dev/null
then
  components=`contour countcc $fpgroup -q`
  command="contour wirtinger $fpgroup -Q | contour ks_$group --list"
  if [ "$components" -gt 1 ]
  then
    if [ -z "$ccemb" ]
    then
      echo "Cannot compute g_image invariant for links ($components components). Use option --ccemb <component>"
      exit 2
    fi
    command="contour wirtinger $fpgroup --ccemb $ccemb -Q | contour ks_$group --list"
    contour wirtinger $fpgroup --ccemb $ccemb -Q | contour ks_$group --list 2>/dev/null >$tmpfile
  else
    contour wirtinger $fpgroup -Q | contour ks_$group --list 2>/dev/null >$tmpfile
  fi
else
  contour ks_$group $fpgroup --list 2>/dev/null >$tmpfile
  command="contour ks_$group $fpgroup --list"
fi

number=`grep "Homomorphism #" $tmpfile | wc -l`
echo " $number found" >&2
cat <<EOT >&2
They can be listed with the command:
 
  $command

EOT

numselected=`grep "^Selected" $tmpfile | cut -f1-3 -d' ' | sort -u | wc -l`
genus=$[ $numselected / 2 ]
if [ "$genus" -lt "1" ]
then
  echo "Cannot continue: genus = $genus"
  exit 2
fi

echo "There are $numselected selected elements" >&2
echo "They are assumed to be listed as 'meridian1,longitude1,meridian2,longitude2,...'" >&2
echo "Handlebody of genus $genus" >&2

if [ -n "$getgenlist" ]
then
  cat $tmpfile | homoimagelist | describegroup | $gap1 | $gap2 | $gap3 | sort -u
  rm $tmpfile
  if [ -f "$tmpinput" ]; then rm $tmpinput; fi
  exit
fi

echo -n "Computing ontolist..." >&2
ontolist=`cat $tmpfile | homoimagelist | describegroup | $gap1 | $gap2 | $gap3 | cat -n | tr -d ' ' | tr -d '"' | tr '\t' ':' | tr -d '"' | grep ":$igroup\$" | cut -f1 -d:`
#ontolist=`cat $tmpfile | homoimagelist | ordergroup | $gap1 | $gap2 | $gap3 | cat -n | tr -d ' ' | tr '\t' ':' | grep ":$gorder\$" | cut -f1 -d:`
numonto=`echo "$ontolist" | wc -l`
echo " $numonto found" >&2

echo -n "Computing properlist..." >&2
properlist=`cat $tmpfile | properparselist | $gap1 | $gap2 | $gap3 | cat -n | tr -d ' ' | tr '\t' ':'  | tr -d '"' | grep -v ":$gorder\$" | cut -f1 -d:`
numproper=`echo "$properlist" | wc -l`
echo " $numproper found" >&2

if [ -z "$onlyhash" ]
then
  echo "num:$number onto:$numonto proper:$numproper"
fi

#for n in $ontolist
#do
#  echo "Onto homomorphism number $n"
#done

echo -n "Computing normal closures..." >&2

if [ -n "$onlyhash" ]
then
  cat $tmpfile | parselist | grep '^Proper ' | sed -e 's/ #[0-9]*: /: /' | sort | md5sum | cut -f1 -d' '
fi

if [ -n "$onlyinvariant" ]
then
  echo "===== cut here ===== hash computed on the section below ====="
  cat $tmpfile | parselist | grep '^Proper ' | sed -e 's/ #[0-9]*: /: /' | sort
fi

if [ -n "$full" ]
then
  echo "A real invariant can be obtained by piping stdout through" >&2
  echo "   | grep '^Proper ' | sed -e 's/ #[0-9]*: /: /' | sort | md5sum | cut -f1 -d' '" >&2
  cat $tmpfile | parselist
fi

echo " done" >&2

#
# clean up temporary file
#
rm $tmpfile
if [ -f "$tmpinput" ]; then rm $tmpinput; fi
