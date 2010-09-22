#!/bin/bash
#
#

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
    case $num in
      0)
        thenum="no"
      ;;
      1)
        thenum="one"
      ;;
      2)
        thenum="two"
      ;;
      3)
        thenum="three"
      ;;
      *)
        thenum=$num
      ;;
    esac

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
  done
}

declare -a holes

if [ -z "$1" ]
then
  echo "usage: $0 <contour-description-file>"
  exit 0
fi

ccontour=`which contour`

if [ "$?" != "0" ]
then
  echo "Cannot find contour software! Check your installation"
  exit 1
fi

file=$1

if [ ! -f $file ]
then
  echo "Cannot find file $1"
fi

if ! $ccontour ishuffman $file -q 2>/dev/null
then
  echo "This is not a contour with Huffman labelling"
  exit 2
fi

cc=`$ccontour countcc $file -q 2>/dev/null`
lista=""

for comp in `seq $cc`
do
  euler=`$ccontour extractcc $comp $file 2>/dev/null | $ccontour characteristic -q 2>/dev/null`
  holes[$comp]=$[ ( 2 - $euler ) / 2 ]
  # echo "component $comp, ${holes[$comp]} holes"
  if [ -z "$lista" ]
  then
    lista="${holes[$comp]}"
  else
    lista="${lista}#${holes[$comp]}"
  fi
done

listanew=`echo "$lista" | tr '#' '\n' | sort -n | uniq -c | tr -s ' '`
totoutput=`echo "$listanew" | wc -l`
totoutput=`echo $totoutput`

if [ "$cc" -gt "1" ]
then
  echo -n "There are "
else
  echo -n "There is "
fi

echo "$listanew" | describe

echo "."
