#!/bin/bash
#
# nella directory contenente l'esito di "trasforma" ci sono
# tutti i risultati delle trasformazioni con un nome iniziante
# con "nnnnn_" dove nnnnn e' un numero a partire da 00001 in
# avanti.
# controlla la corrispondenza della conversione in morse

list=`ls | grep "^[0-9]"`

for s1 in $list
do
  echo -n "checking morse and region equivalency for $s1:"
  if ( cat $s1 ; contour printmorse $s1 2>/dev/null ) | contour compare 2>/dev/null >/dev/null
  then
    echo OK
  else
    echo FAIL
    exit 1
  fi
done
