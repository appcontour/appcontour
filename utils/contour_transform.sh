#!/bin/bash
#

#
####  funzioni  #####
#
function incasella
{
  nuovofile=$1
  pattern=$2
  dir=$3

#  echo "Incasello $nuovofile in $dir..."

  cmpfile=`ls $dir/${pattern}* 2>/dev/null | tail -1`
  if [ -n "$cmpfile" ]
  then
    cat $nuovofile $cmpfile | contour compare 2>/dev/null >/dev/null
    res=$?
    case $res in
      0) mv $nuovofile $dir
        return 1
        ;;
      1)
        if [ ! -d $dir/.worse ]
        then
          mkdir $dir/.worse
        fi
        incasella $nuovofile $pattern $dir/.worse
        return $?
        ;;
      2)
        if [ ! -d $dir/.better ]
        then 
          mkdir $dir/.better
        fi
        incasella $nuovofile $pattern $dir/.better
        return $?
        ;;
      *) echo "ERRORE"
        exit 4
        ;;
    esac
  else
    cp $nuovofile $dir
    return 0
  fi
}

counter=1
function appiattisci
{
  local dir dest
  dir=$1
  dest=$2

  if [ -d $dir/.better ]
  then
    appiattisci $dir/.better $dest
    rmdir $dir/.better
  fi

  for file in `ls $dir`
  do
    bfile=`basename $file`
#    echo "Moving $dir/$file as number $counter"
    cntr=`printf "%05d" $counter`
    mv $dir/$file $dest/${cntr}_$bfile
  done
  counter=$[ $counter + 1 ]

  if [ -d $dir/.worse ]
  then
    appiattisci $dir/.worse $dest
    rmdir $dir/.worse
  fi
}

#
#######   fine funzioni  ########
#

if [ -z "$1" ]
then
  echo "usage: $0 sketch [numero-mosse]"
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

mosse=100
if [ -n "$2" ]
then
  mosse=$2
fi

mosse=$[ mosse + 1 ]
nomeesempio=`basename $esempio ".morse"`

fatti=$nomeesempio.transformations

myname=`basename $0`
mydir=`dirname $0`

pre="${mydir}/"
if [ "$0" = "$myname" ]
then
  pre=""
fi

basta=""
if contour print $esempio 2>/dev/null | grep -q "\[no info"
then
  basta="1"
  echo "Some informations on arcs are missing for $esempio:"
  contour print $esempio
fi

if ! contour isappcon $esempio 2>/dev/null >/dev/null 
then
  basta="1"
  contour isappcon $esempio
fi

if [ -e deposito ]
then
  basta="1"
  echo "you must first remove directory 'deposito'"
fi

if [ -e $fatti ]
then
  echo "you must firse remove directory $fatti"
  basta="1"
fi

if [ -n "$basta" ]
then
  if [ -n "$removemorse" ]
  then
    rm $esempio
  fi
  exit 1
fi

mkdir deposito
mkdir $fatti

#cp $esempio deposito/$nomeesempio
contour print $esempio >deposito/$nomeesempio

#
# ci deve essere un deposito da cui pescare gli sketch da
# trasformare

cd deposito
mkdir .albero
incasella $nomeesempio $nomeesempio .albero
goon=1
while [ -n "$goon" ]
do
  goon=""
  for file in `ls`
  do
    goon=1
    tail=${file#$nomeesempio}
    rr=`echo "$tail" | cut -f${mosse} -d.`
    if [ -z "$rr" ]
    then
#      echo "Trasformo $file..."
      rules=`contour testallrules <$file 2>/dev/null | tail -1`

#      echo rules=$rules
      for r in $rules
      do
        rmod=`echo $r | tr ":" "-"`
        echo "Transforming $file with rule $r"
        contour applyrule $r <$file 2>/dev/null >$file.$rmod
        incasella $file.$rmod $nomeesempio .albero
      done
    fi
    rm $file
  done
done

appiattisci .albero ../$fatti
rmdir .albero
cd ..
rmdir deposito

if [ -n "$removemorse" ]
then
  echo "removing $esempio"
  rm $esempio
fi
