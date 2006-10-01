#!/bin/bash
#
# analizza una scena descrivendo il tipo topologico 
# delle superfici

if [ -z "$1" ]
then
  echo "usage: $0 sketch"
  exit 1
fi

esempio=$1
if [ ! -f "$esempio" ]
then
  echo "non trovo il file $esempio"
  exit 2
fi

removemorse=""
if grep -q "^knot " $esempio
then
  echo "This seems a knot description, creating a morse tubular description"
  nomeesempio=`basename $esempio ".knot"`
  morsefile="${nomeesempio}.morse"
  if [ -f "$morsefile" ]
  then
    echo "file $morsefile already exists, cannot proceed"
    exit 2
  fi
  contour knot2morse $esempio >$morsefile
  esempio="$morsefile"
  removemorse="1"
fi

ccnum=`contour -q countcc $esempio`

e="e"
a="a"
if [ "$ccnum" -gt 1 ]
then
  e="i"
  a="e"
fi

if [ "$ccnum" -gt 1 ]
then
  echo "La superficie e' composta da $ccnum componenti connesse."
  for cc in `seq 1 $ccnum`
  do
    charac=`contour extractcc $cc $esempio | contour -q characteristic 2>/dev/null`
    echo "La numero $cc ha una caratteristica di Eulero pari a $charac"
  done
else
  charac=`contour -q characteristic $esempio`
  echo "La superficie e' connessa, e la sua caratteristica di Eulero e' $charac"
fi


