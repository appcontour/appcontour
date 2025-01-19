#!/bin/bash
#

num=$1
choice=$2

function usage ()
{
  echo "usage: $0 <num> [<choice>] <embeddings-file"
}

if [ -z "$num" ]
then
  usage
fi

if [ -z "$choice" ]
then
  choice=0
  echo "# using default value of 0 for choice"
fi

emb=`grep "^#$num " | cut -f2- -d' '`

echo "embedding:$choice $emb"
