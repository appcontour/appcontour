#!/bin/bash
#
example=$1

showcontouroptions="--steps 4"

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
    commandchain="$commandchain | contour applyrule $rule 2>/dev/null"
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
  applicable=`eval $commandchain | contour testallrules 2>/dev/null | tail -1`
  echo "Applicable rules: $applicable"
  applicable=" $applicable "
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
      show
      ;;
    back)
      back
      ;;
    exit|quit)
      exit
      ;;
    help)
      echo "Valid commands are:"
      echo -n "help, quit, print, info, morse, back"
      if [ -n "$showcontour" ]
      then
        echo ", show"
      else
        echo ""
      fi
      ;;
    *)
    if echo "$applicable" | grep -qi " $command "
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
