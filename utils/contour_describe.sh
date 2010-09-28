#!/bin/bash
#
#

function usage ()
{
  echo "usage: $0 [-z] <contour-description-file>"
}

function mynumber ()
{
  lnum=$1
  case $lnum in
    0)
      echo "no"
    ;;
    1)
      echo "one"
    ;;
    2)
      echo "two"
    ;;
    3)
      echo "three"
    ;;
    *)
      if [ -n "$nnumber" ]
      then
        $nnumber -l $lnum
      else
        echo $lnum
      fi
    ;;
  esac
}

function myplural()
{
  if [ "$1" = 1 ]
  then
    echo "$2"
  else
    echo "$3"
  fi
}

function describe ()
{
  numoutput="0"
  while read num holes
  do
    numoutput=$[ $numoutput + 1 ]
    separator=""
    if [ "$numoutput" -gt 1 ]
    then
      if [ "$numoutput" = "$totoutput" ]
      then
        separator=" and "
      else
        separator=", "
      fi
    fi
    echo -n "$separator"
    plural="1"
    if [ "$num" = "1" ]; then plural=""; fi
    thenum=`mynumber $num`

    case $holes in
      0)
        object="sphere"
        if [ -n "$plural" ]; then object="spheres"; fi
        ;;
      1)
        object="torus"
        if [ -n "$plural" ]; then object="tori"; fi
        ;;
      *)
        object="sphere with $holes handles"
        if [ -n "$plural" ]; then object="spheres with $holes handles"; fi
        ;;
    esac
    echo -n $thenum $object
    negatives=`echo "$lista" | grep "^${holes}:-" | wc -l`
    negatives=`echo $negatives`
    if [ "$negatives" -gt 0 ]
    then
      if [ "$num" = 1 ]
      then
        echo -n " (negatively oriented)"
      else if [ "$num" = "$negatives" ]
        then
          echo -n " (all negatively oriented)"
        else
          negis=`myplural $negatives is are`
          theneg=`mynumber $negatives`
          echo -n " ($theneg of which $negis negatively oriented)"
        fi
      fi
    fi
  done
}

function zorkdescribe ()
{
  ccid=$1
  count=0

  list=`$ccontour -q ccchilds $ccid $file 2>/dev/null`
  if [ -n "$list" ]
  then
    if [ "$ccid" = "0" ]
    then
      echo -n "You can see"
    else
      echo -n "The ${description[$ccid]} contains"
    fi
    objnum=`echo $list | wc -w`
    for o in $list
    do
      count=$[ $count + 1 ]
      if [ "$count" -gt "1" ]
      then
        if [ "$count" -lt "$objnum" ]
        then
          echo -n ","
        else
          echo -n " and"
        fi
      fi
      echo -n " a ${description[$o]}"
    done
    echo "."
  else
    if [ -n "$list" ]
    then
      echo "You can see nothing."
#    else
#      echo "The ${description[$ccid]} contains nothing."
    fi
  fi
  for o in $list
  do
    zorkdescribe $o
  done
}

declare -a holes
declare -a description
declare -a objwithgivenholes
declare -a digits
declare -a counter
declare -a adjective

adjective[0]="white"
adjective[1]="black"
adjective[2]="red"
adjective[3]="blue"
adjective[4]="green"
adjective[5]="purple"
adjective[6]="brown"
adjective[7]="yellow"
adjective[8]="pink"
adjective[9]="violet"

material[0]="metal"
material[1]="stone"
material[2]="marble"
material[3]="wooden"
material[4]="granite"
material[5]="gold"
material[6]="silver"
material[7]="bronze"
material[8]="plastic"
material[9]="sandstone"

zork=""

while getopts "z" option
do
  case $option in
  z)
    zork=1
    ;;
  *)
    echo "Invalid option"
    usage
    exit 1
    ;;
  esac
done

shift $(( $OPTIND - 1 ))

file=$1
if [ -z "$file" ]
then
  usage
  exit 0
fi

ccontour=`which contour`
if [ "$?" != "0" ]
then
  echo "Cannot find contour software! Check your installation"
  exit 1
fi

nnumber=`which number`
if [ "$?" != "0" ]
then
  nnumber=""
fi

if [ ! -f $file ]
then
  echo "Cannot find file $1"
fi

if ! $ccontour ishuffman $file -q 2>/dev/null
then
  if [ -n "$zork" ]
  then
    echo "It is pitch black, you are likely to be eaten by a grue."
  else
    echo "This is not a contour with Huffman labelling"
  fi
  exit 2
fi

cc=`$ccontour countcc $file -q 2>/dev/null`
lista=""

for comp in `seq $cc`
do
  euler=`$ccontour extractcc $comp $file 2>/dev/null | $ccontour characteristic -q 2>/dev/null`
  nholes=$[ ( 2 - $euler ) / 2 ]
  holes[$comp]=$nholes
  count=${objwithgivenholes[$nholes]}
  digits[$nholes]=${#count}
  counter[$nholes]="0"
  if [ -z "$count" ]; then count="0"; fi
  count=$[ $count + 1 ]
  objwithgivenholes[$nholes]=$count
  orientation=`$ccontour -q ccorientation $comp $file 2>/dev/null`
  # echo "component $comp, ${holes[$comp]} holes"
  if [ -z "$lista" ]
  then
    lista="${holes[$comp]}"
  else
    lista="${lista}#${holes[$comp]}:$orientation"
  fi
done

for comp in `seq $cc`
do
  nholes=${holes[$comp]}
  ldigits=${digits[$nholes]}
  count=${counter[$nholes]}
  lcount=$count
  adj=""
  if [ "$ldigits" -gt 0 ]
  then
    idx=$[ $lcount % 10 ]
    lcount=$[ $lcount / 10 ]
    adj=${adjective[$idx]}
    ldigits=$[ $ldigits - 1 ]
  fi
  if [ "$ldigits" -gt 0 ]
  then
    idx=$[ $lcount % 10 ]
    lcount=$[ $lcount / 10 ]
    adj="$adj ${material[$idx]}"
    ldigits=$[ $ldigits - 1 ]
  fi
  for iter in `seq $ldigits`
  do
    idx=$[ $lcount % 10 ]
    lcount=$[ $lcount / 10 ]
    adj="${adjective[$idx]} $adj"
    ldigits=$[ $ldigits - 1 ]
  done
  case $nholes in
    0)
      descr="sphere"
      ;;
    1)
      descr="torus"
      ;;
    *)
      theholes=`mynumber $nholes`
      descr="sphere with $theholes handles"
      ;;
  esac
  if [ -n "$adj" ]
  then
    descr="$adj $descr"
  fi
  description[$comp]="$descr"
  count=$[ $count + 1 ]
  counter[$nholes]=$count
done

lista=`echo "$lista" | tr '#' '\n' | sort -n`
listanew=`echo "$lista" | cut -f1 -d':' | uniq -c | tr -s ' '`
totoutput=`echo "$listanew" | wc -l`
totoutput=`echo $totoutput`

if [ -n "$zork" ]
then
  echo "Clearing"
  echo "You are in a clearing, with a forest surrounding you on all sides."
  case $cc in
    0)
      echo "There are no objects here."
      ;;
    1)
      echo "There is one object here."
      echo ""
      ;;
    *)
      thecc=`mynumber $cc`
      echo "There is a total of $thecc objects here."
      echo ""
      ;;
  esac
fi

if [ -n "$zork" ]
then
  zorkdescribe 0
else
  if [ "$cc" -gt "1" ]
  then
    echo -n "There are "
  else
    echo -n "There is "
  fi

  echo "$listanew" | describe

  echo "."
fi
