#!/bin/bash
#
origexample=$1
example=${origexample}

searchdirs="/usr/local/share/appcontour/examples /usr/share/appcontour/examples"
extensions=".morse .sketch .knot"

if [ ! -f "$example" ]
then
  # first search in current directory
  for ext in $extensions
  do
    if [ -f ${example}${ext} ]
    then
      example=${example}${ext}
      break
    fi
  done
fi
if [ ! -f "$example" ]
then
  for d in $searchdirs
  do
    if [ -f ${d}/${example} ]
    then
      example=${d}/${example}
      break
    fi
    for ext in $extensions
    do
      if [ -f ${d}/${example}${ext} ]
      then
        example=${d}/${example}${ext}
        break
      fi
    done
    if [ -f "$example" ]
    then
      break
    fi
  done
fi

if [ "$example" != "$origexample" ]
then
  echo "Found matching file: $example"
fi

#ns="--oldnames"
ns="--newnames"

showcontouroptions="--steps 4"

function listma ()
{
  eval $commandchain | contour $ns listma 2>/dev/null
}

function listwr ()
{
  eval $commandchain | contour $ns listwr 2>/dev/null
}

function listst ()
{
  eval $commandchain | contour $ns listinvcn1 2>/dev/null
}

function listpu ()
{
  eval $commandchain | contour $ns listinvcn3 2>/dev/null
}

function displayinfo ()
{
  eval $commandchain | contour info 2>/dev/null
}

function printexample ()
{
  eval $commandchain | contour print 2>/dev/null
}

function morse ()
{
  eval $commandchain | contour printmorse 2>/dev/null
}

function show ()
{
  eval $commandchain | contour printmorse 2>/dev/null | $showcontour --title "$examplename $rules" 2>/dev/null &
}

function back ()
{
  prevrules="$rules"
  prevrule=""
  rules=""
  for r in $prevrules
  do
    if [ -n "$prevrule" ]
    then
      rules="$rules $prevrule"
    fi
    prevrule="$r"
  done
}

function buildcommandchain ()
{
  commandchain="$catcommand"
  for rule in $rules
  do
    commandchain="$commandchain | contour $ns applyrule $rule 2>/dev/null"
  done
}

if [ -z "$example" ]
then
  echo "usage: $0 example"
  exit 1
fi

if [ ! -f "$example" ]
then
  echo "cannot find example file $example"
  exit 1
fi

examplefname=`basename $example`
examplename=${examplefname%.morse}
examplename=${examplename%.sketch}
examplename=${examplename%.knot}

if [ "$examplefname" == "$examplename".knot ]
then
  echo "This is a knot format, converting to morse..."
  catcommand="contour knot2morse $example"
else
  catcommand="cat $example"
fi

#
# try see if we have showcontour...
#
paths="/usr/local/bin . ../show ../appcontour/show"
showcontour=""
for p in $paths
do
  if [ -x "$p/showcontour" ]
  then
    grident=`$p/showcontour --grident`
    if [ "$grident" != "null" ]
    then
      showcontour="$p/showcontour $showcontouroptions"
      echo "found showcontour [$grident] in $p"
    fi
    break
  fi
done

echo "examplename: $examplename"

rules=""

buildcommandchain
displayinfo

while true
do
  if [ -n "$rules" ]
  then
    echo "Applied rules: $rules"
  fi
  buildcommandchain
  applicablesimple=`eval $commandchain | contour $ns testallrules 2>/dev/null | tail -1`
  applicablema=`eval $commandchain | contour $ns listma -q 2>/dev/null | tail -1`
  applicablewr=`eval $commandchain | contour $ns listwr -q 2>/dev/null | tail -1`
  applicablest=`eval $commandchain | contour $ns listinvcn1 -q 2>/dev/null | tail -1`
  applicablepu=`eval $commandchain | contour $ns listinvcn3 -q 2>/dev/null | tail -1`
  echo "Applicable rules: $applicablesimple"
  #echo "Applicable mergearcs rules: $applicablema"
  applicable=" $applicablesimple $applicablema $applicablewr $applicablest $applicablepu "
  set -o history
  read -e -p "Contour> " command arg
  if [ "$?" != "0" ]
  then
    exit
  fi
  if [ -z "$command" ]
  then
    continue
  fi

  history -s $command $arg

  case $command in
    listst|swallowtail)
      listst
      ;;
    listpu|puncture)
      listpu
      ;;
    listwr|wrinkle)
      listwr
      ;;
    listma|mergearcs)
      listma
      ;;
    info)
      displayinfo
      ;;
    print)
      printexample
      ;;
    morse)
      morse
      ;;
    show)
      if [ -n "$arg" ]
      then
        c1=`echo "$arg" | cut -f1 -d':'`
        c2=`echo "$arg" | cut -s -f2 -d':'`
        if [ -n "$c2" ]
        then
          c1="${c1}:${c2}"
        fi
        if echo "$applicable" | grep -qi " $c1 "
        then
          echo "OK, applying rule $arg would result in this"
          savedrules="$rules"
          rules="$rules $arg"
          buildcommandchain
          show
          rules="$savedrules"
          buildcommandchain
        else
          echo "$arg is not an applicable rule"
        fi
      else
        show
      fi
      ;;
    back)
      back
      ;;
    exit|quit)
      exit
      ;;
    help)
      echo "Valid commands are:"
      echo "help, quit, print, info, morse, back"
      echo -n "wrinkle, mergearcs, swallowtail, puncture"
      if [ -n "$showcontour" ]
      then
        echo ", show [rule]"
      else
        echo ""
      fi
      ;;
    *)
    c1=`echo "$command" | cut -f1 -d':'`
    c2=`echo "$command" | cut -s -f2 -d':'`
    if [ -n "$c2" ]
    then
      c1="${c1}:${c2}"
    fi
    if echo "$applicable" | grep -qi " $c1 "
    then
      echo "OK, applying rule $command"
      rules="$rules $command"
      buildcommandchain
#      echo "new command chain: $commandchain"
#      eval $commandchain
    else
      echo "Invalid command $command, type help for help"
    fi
    ;;
  esac
done
done
