#!/bin/bash
#
#

function usage ()
{
  echo "usage: $0 [-o][-z] <contour-description-file>"
  echo "  option -z gives zork-like description with orientation information"
  echo "  option -o gives a description with orientation information"
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
      if [ -n "$zork" ]
      then
        echo -n "You can see"
      else
        echo "External components:"
      fi
    else
      if [ -n "$zork" ]
      then
        echo -n "The ${description[$ccid]} contains"
      else
        echo "Component #$ccid contains:"
      fi
    fi
    objnum=`echo $list | wc -w`
    for o in $list
    do
      count=$[ $count + 1 ]
      if [ "$count" -gt "1" -a -n "$zork" ]
      then
        if [ "$count" -lt "$objnum" ]
        then
          echo -n ","
        else
          echo -n " and"
        fi
      fi
      if [ -n "$zork" ]
      then
        echo -n " a ${longdescription[$o]}"
      else
        echo " ${description[$o]}"
      fi
    done
    if [ -n "$zork" ]
    then
      echo "."
    fi
  else
    if [ -n "$list" ]
    then
      if [ -n "$zork" ]
      then
        echo "You can see nothing."
      else
        echo "Empty apparent contour!"
      fi
#    else
#      echo "The ${description[$ccid]} contains nothing."
    fi
  fi
  for o in $list
  do
    zorkdescribe $o
  done
}

function zorkdescribepairs ()
{
  for comp1 in `seq $cc`
  do
    if [ "${holes[$comp1]}" != 1 ]; then continue; fi
    parent1=`$ccontour -q ccparent $comp1 $file 2>/dev/null`
    for comp2 in `seq $comp1 $cc`
    do
      if [ "$comp1" = "$comp2" ]; then continue; fi
      if [ "${holes[$comp2]}" != 1 ]; then continue; fi
      # sono fratelli? sono padre/figlio?
      parent2=`$ccontour -q ccparent $comp2 $file 2>/dev/null`
      fgoption=""
      parent=""
      if [ "$parent1" = "$parent2" ]
      then
        fgoption="--outside"
      fi
      if [ "$comp1" = "$parent2" ]
      then
        fgoption="--inside"
        parent=$comp1
        child=$comp2
      fi
      if [ "$comp2" = "$parent1" ]
      then
        fgoption="--inside"
        parent=$comp2
        child=$comp1
      fi
      if [ -z "$fgoption" ]; then continue; fi
      fg=`$ccontour -q extractcc $comp1,$comp2 $file 2>/dev/null | $ccontour -q fg $fgoption 2>/dev/null`
      relators=`echo $fg | cut -f2 -d';' | cut -f1 -d'>'`
      relators=`echo $relators`
      if [ -z "$relators" ]
      then
        echo "The ${description[$comp1]} and the ${description[$comp2]} are not linked."
        continue
      fi
      linking=`$ccontour -q extractcc $comp1,$comp2 $file 2>/dev/null | $ccontour -q linkingnumber $fgoption 2>/dev/null`
      if [ "$linking" = "0" ]
      then
        alexander=`getalexanderlink2 $comp1 $comp2 $file $fgoption`
        apparently="linked"
        if [ -z "$alexander" ]; then apparently="apparently linked"; fi
        if [ -n "$parent" ]
        then
          echo "The ${description[$child]} is $apparently inside the ${description[$parent]}."
        else
          echo "The ${description[$comp1]} is $apparently with the ${description[$comp2]}."
        fi
        continue
      fi
      case $linking in
        1)
        times="once"
        ;;

        2)
        times="twice"
        ;;

        3)
        times="three times"
        ;;

        *)
        times="many times"
        ;;
      esac
      if [ -n "$parent" ]
      then
        echo "The ${description[$child]} is linked $times into the ${description[$parent]}."
      else
        echo "The ${description[$comp1]} is linked $times to the ${description[$comp2]}."
      fi
    done
  done
}

function getalexander ()
{
  comp=$1
  file=$2
  foxd=$3
  sideopt=$4

  if [ -z "$foxd" ]; then foxd=1; fi
  alexander=`$ccontour extractcc $comp $file 2>/dev/null | $ccontour alexander --foxd $foxd $sideopt -q 2>/dev/null | grep -iv warning`
  alexander=`echo $alexander | tr -d ' '`
  alexander=${alexander#+}
  if [ "$alexander" = "1" ]; then alexander=""; fi
  if [ "$alexander" = "-1" ]; then alexander=""; fi
  if [ "$alexander" = "[+1]" ]; then alexander=""; fi
  if [ "$alexander" = "[1]" ]; then alexander=""; fi
  if [ "$alexander" = "[-1]" ]; then alexander=""; fi
  if [ "$alexander" = "-[+1]" ]; then alexander=""; fi
  echo $alexander
}

function getalexanderlink2 ()
{
  comp1=$1
  comp2=$2
  file=$3
  sideopt=$4

  alexander=`$ccontour extractcc $comp1,$comp2 $file 2>/dev/null | $ccontour alexander --foxd 1 $sideopt -q 2>/dev/null | grep -iv warning`
  alexander=`echo $alexander | tr -d ' '`
  alexander=${alexander#+}
  if [ "$alexander" = "0" ]; then alexander=""; fi
  if [ "$alexander" = "-0" ]; then alexander=""; fi
  if [ "$alexander" = "[+0]" ]; then alexander=""; fi
  if [ "$alexander" = "[0]" ]; then alexander=""; fi
  if [ "$alexander" = "[-0]" ]; then alexander=""; fi
  if [ "$alexander" = "-[+0]" ]; then alexander=""; fi
  echo $alexander
}

declare -a holes
declare -a description
declare -a longdescription
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

zork="1"
orient=""

while getopts "zno" option
do
  case $option in
  z)
    zork=1
    ;;
  n)
    zork=""
    ;;
  o)
    orient="1"
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

#
# number command is contained in the "bsd-games" package
# and converts a number to english words, like one, two, three,...
#
nnumber=`which number 2>/dev/null`
if [ "$?" != "0" ]
then
  nnumber=""
fi

$ccontour ishuffman $file -q 2>/dev/null
status=$?
if [ "$status" = "10" ]
then
  echo "Cannot find file $file"
  exit $status
fi

if [ "$status" != "0" ]
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
    lista="${holes[$comp]}:$orientation"
  else
    lista="${lista}#${holes[$comp]}:$orientation"
  fi
done

for comp in `seq $cc`
do
  nholes=${holes[$comp]}
  if [ -n "$zork" ]
  then
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
    longdescr=""
    ialexander=""
    oalexander=""
    case $nholes in
      0)
        descr="sphere"
        ;;
      1)
        descr="torus"
        ifg=`$ccontour extractcc $comp $file 2>/dev/null | $ccontour fg --inside -q 2>/dev/null | tr -d ' ' | cut -f2 -d';' | cut -f1 -d'>'`
        ofg=`$ccontour extractcc $comp $file 2>/dev/null | $ccontour fg --outside -q 2>/dev/null | tr -d ' ' | cut -f2 -d';' | cut -f1 -d'>'`
        if [ -n "$ofg" ]; then oalexander=`getalexander $comp $file $nholes --outside`; fi
        if [ -n "$ifg" ]; then ialexander=`getalexander $comp $file $nholes --inside`; fi
        if [ -n "$ofg" -a -z "$ifg" ]
        then
          longdescr="torus (apparently knotted)"
          if [ -n "$oalexander" ]
          then
            longdescr="knotted torus"
          fi
        fi
        if [ -n "$ifg" -a -z "$ofg" ]
        then
          longdescr="torus (with an apparently knotted hole)"
          if [ -n "$ialexander" ]
          then
            longdescr="torus with a knotted hole"
          fi
        fi
        if [ -n "$ifg" -a -n "$ofg" ]
        then
          longdescr="mysteriously tangled torus"
        fi
        ;;
      2)
        descr="double torus"
        ifg=`$ccontour extractcc $comp $file 2>/dev/null | $ccontour fg --inside -q 2>/dev/null | tr -d ' ' | cut -f2 -d';' | cut -f1 -d'>'`
        ofg=`$ccontour extractcc $comp $file 2>/dev/null | $ccontour fg --outside -q 2>/dev/null | tr -d ' ' | cut -f2 -d';' | cut -f1 -d'>'`
        if [ -n "$ofg" ]; then oalexander=`getalexander $comp $file $nholes --outside`; fi
        if [ -n "$ifg" ]; then ialexander=`getalexander $comp $file $nholes --inside`; fi
        if [ -n "$ofg" -a -z "$ifg" ]
        then
          longdescr="double torus (apparently knotted)"
          if [ -n "$oalexander" ]
          then
            longdescr="knotted double torus"
          fi
        fi
        if [ -n "$ifg" -a -z "$ofg" ]
        then
          longdescr="doble torus (with an apparently knotted hole)"
          if [ -n "$ialexander" ]
          then
            longdescr="double torus with knotted holes"
          fi
        fi
        if [ -n "$ifg" -a -n "$ofg" ]
        then
          longdescr="mysteriously tangled double torus"
        fi
        ;;
      *)
        theholes=`mynumber $nholes`
        descr="sphere with $theholes handles"
        ;;
    esac
    if [ -n "$adj" ]
    then
      descr="$adj $descr"
      if [ -n "$longdescr" ]; then longdescr="$adj $longdescr"; fi
    fi
  else
    descr="component #$comp, genus $nholes"
  fi
  description[$comp]="$descr"
  longdescription[$comp]="$descr"
  if [ -n "$longdescr" ]
  then
    if [ -n "$oalexander" ]
    then
      longdescr="${longdescr}, you can read \"Alexander: $oalexander\" written on it"
    fi
    if [ -n "$ialexander" ]
    then
      longdescr="${longdescr}, you can read \"Alexander: $ialexander\" written in the inside"
    fi
  fi
  if [ -n "$longdescr" ]; then longdescription[$comp]="$longdescr"; fi
  count=$[ $count + 1 ]
  counter[$nholes]=$count
done

lista=`echo "$lista" | tr '#' '\n' | sort -n`
listanew=`echo "$lista" | cut -f1 -d':' | uniq -c | tr -s ' '`
totoutput=`echo "$listanew" | wc -l`
totoutput=`echo $totoutput`

if [ -z "$orient" ]
then
  if [ -n "$zork" ]
  then
    echo "You are in a clearing, with a forest surrounding you on all sides."
  fi
  case $cc in
    0)
      if [ -n "$zork" ]
      then
        echo "There are no objects here."
      else
        echo "Empty apparent contour!"
      fi
      ;;
    1)
      if [ -n "$zork" ]
      then
        echo "There is one object here."
        echo ""
      else
        echo "Total number of connected components: $cc"
      fi
      ;;
    *)
      if [ -n "$zork" ]
      then
        thecc=`mynumber $cc`
        echo "There is a total of $thecc objects here."
        echo ""
      else
        echo "Total number of connected components: $cc"
      fi
      ;;
  esac
fi

if [ -z "$orient" ]
then
  zorkdescribe 0
  zorkdescribepairs
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
