#!/bin/bash
#
#

function usage ()
{
  echo "usage: $0 [-v] <contour-description-file>"
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
#    case $num in
#      0)
#        thenum="no"
#      ;;
#      1)
#        thenum="one"
#      ;;
#      2)
#        thenum="two"
#      ;;
#      3)
#        thenum="three"
#      ;;
#      *)
#        thenum=$num
#        if [ -n "$nnumber" ]
#        then
#          thenum=`$nnumber -l $num`
#        fi
#      ;;
#    esac
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
          echo -n " (of which $theneg $negis negatively oriented)"
        fi
      fi
    fi
  done
}

declare -a holes
verbose=""

while getopts "v" option
do
  case $option in
  v)
    verbose=1
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
  echo "This is not a contour with Huffman labelling"
  exit 2
fi

cc=`$ccontour countcc $file -q 2>/dev/null`
lista=""

for comp in `seq $cc`
do
  euler=`$ccontour extractcc $comp $file 2>/dev/null | $ccontour characteristic -q 2>/dev/null`
  holes[$comp]=$[ ( 2 - $euler ) / 2 ]
  orientation=`$ccontour -q ccorientation $comp $file 2>/dev/null`
  # echo "component $comp, ${holes[$comp]} holes"
  if [ -z "$lista" ]
  then
    lista="${holes[$comp]}"
  else
    lista="${lista}#${holes[$comp]}:$orientation"
  fi
done

lista=`echo "$lista" | tr '#' '\n' | sort -n`
listanew=`echo "$lista" | cut -f1 -d':' | uniq -c | tr -s ' '`
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
