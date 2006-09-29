#!/bin/bash
#
# nella directory contenente l'esito di "trasforma" ci sono
# tutti i risultati delle trasformazioni con un nome iniziante
# con "nnnnn_" dove nnnnn e' un numero a partire da 00001 in
# avanti.  In teoria due sketch s1 e s2 dovrebbero essere
# s1 < s2; s1 = s2 o s1 > s2 se rispettivamente i numeri del
# nome (n1 e n2) sono nello stesso ordine.
#
# controllo che questo sia vero!

list=`ls [0-9]*`

for s1 in $list
do
  echo "comparing $s1 with all others"
  for s2 in $list
  do
    cat $s1 $s2 | contour compare >/dev/null 2>/dev/null
    res=$?
    n1=`echo "$s1" | cut -f1 -d_`
    n2=`echo "$s2" | cut -f1 -d_`
    case $res in
      0)
        if [ "$n1" != "$n2" ]
        then
          echo "TROVATA DIFFERENZA: $s1, $s2"
        fi
        ;;
      1)
        if [ "$n1" -le "$n2" ]
        then
          echo "TROVATA DIFFERENZA: $s1, $s2"
        fi
        ;;
      2)
        if [ "$n1" -ge "$n2" ]
        then
          echo "TROVATA DIFFERENZA: $s1, $s2"
        fi
        ;;
      *)
        echo "exitcode invalido: $res, s1 = $s1, s2 = $s2"
        exit 1
        ;;
    esac
  done
done
