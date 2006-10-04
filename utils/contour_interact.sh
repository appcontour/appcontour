#!/bin/bash
#
example=$1

function displayinfo ()
{
  eval $commandchain | contour info 2>/dev/null
}

function printexample ()
{
  eval $commandchain | contour print 2>/dev/null
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

catcommand="cat $example"

examplename=`basename $example ".morse"`
examplename=${examplename%.sketch}

echo "examplename: $examplename"

rules=""

buildcommandchain
displayinfo

while true
do
  echo "applied rules: $rules"
  buildcommandchain
  applicable=`eval $commandchain | contour testallrules 2>/dev/null | tail -1`
  echo "Applicable rules: $applicable"
  applicable=" $applicable "
  echo -n "command: "
  read command arg
  if [ "$?" != "0" ]
  then
    exit
  fi
  if [ -z "$command" ]
  then
    continue
  fi

  case $command in
    info)
      displayinfo
      ;;
    print)
      printexample
      ;;
    back)
      back
      ;;
    exit|quit)
      exit
      ;;
    help)
      echo "Valid commands are:"
      echo "help, exit, print, back"
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
