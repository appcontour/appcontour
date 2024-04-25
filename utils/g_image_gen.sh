#!/bin/bash
#

getgenlist=""

case $1 in
  -l)
    getgenlist="yes"
    shift
    ;;
esac

fpgroup=$1
group=$2

echo "Numbers correspond to the order of the group"

if [ -z "$2" ]
then
  echo "usage: $0 some_wirtinger_presentation.fpgroup group"
  echo "where group can be e.g. A5"
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
#

tmpfile=/tmp/g_image.$$.homos
gap1="gap -r -b -q"
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
  A7)
    gorder=2520
    ;;
  *)
    echo "Cannot deal with group $group"
    exit 3
    ;;
esac

if grep -q "^embedding" $fpgroup
then
  components=`contour countcc $fpgroup -q`
  if [ "$components" -gt 1 ]
  then
    echo "Cannot compute g_image invariant for links ($components components)"
    exit 2
  fi
  contour wirtinger $fpgroup -Q | contour ks_$group -v 2>/dev/null >$tmpfile
  command="contour wirtinger $fpgroup -Q | contour ks_$group -v"
else
  contour ks_$group $fpgroup -v 2>/dev/null >$tmpfile
  command="contour ks_$group $fpgroup -v"
fi

number=`grep "Homomorphism #" $tmpfile | wc -l`
numselected=`grep "^Selected" $tmpfile | cut -f1-3 -d' ' | sort -u | wc -l`
genus=$[ $numselected / 2 ]
if [ "$genus" -lt "1" ]
then
  echo "Cannot continue: genus = $genus"
  exit 2
fi

echo "There are $numselected selected elements"
echo "They are assumed to be listed as 'meridian1,longitude1,meridian2,longitude2,...'"
echo "Handlebody of genus $genus"

if [ -n "$getgenlist" ]
then
  cat $tmpfile | homoimagelist | describegroup | $gap1 | $gap2 | $gap3 | sort -u
  rm $tmpfile
  exit
fi

ontolist=`cat $tmpfile | homoimagelist | describegroup | $gap1 | $gap2 | $gap3 | cat -n | tr -d ' ' | tr -d '"' | tr '\t' ':' | tr -d '"' | grep ":$igroup\$" | cut -f1 -d:`
#ontolist=`cat $tmpfile | homoimagelist | ordergroup | $gap1 | $gap2 | $gap3 | cat -n | tr -d ' ' | tr '\t' ':' | grep ":$gorder\$" | cut -f1 -d:`
properlist=`cat $tmpfile | properparselist | $gap1 | $gap2 | $gap3 | cat -n | tr -d ' ' | tr '\t' ':'  | tr -d '"' | grep -v ":$gorder\$" | cut -f1 -d:`
numonto=`echo "$ontolist" | wc -l`
numproper=`echo "$properlist" | wc -l`

cat <<EOT >&2
There are $number Homomorphisms of which $numonto are onto and $numproper are proper.
They can be listed with the command:
 
  $command

EOT


#for n in $ontolist
#do
#  echo "Onto homomorphism number $n"
#done

cat $tmpfile | parselist

#
# clean up temporary file
#
rm $tmpfile
