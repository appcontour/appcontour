#!/bin/bash
#
# a quick way to obtain a complex set of generators that end up
# with an ideal generate by just two polynomials

giove="
ideal(t) {
+4t-3;
+2t+3;
}"

echo "$giove" | ./contour --experimental alexander

#
# this is the ideal (no simplification) from kanenobu -2 8
#
idealkanm2_8="
#
# --foxd 2
#
# *** Warning: result can be noncanonical ***
ideal(t) {
+4-7t+4t^2;
+4-6t+t^2+t^3;
+4-6t+t^2+t^3;
-4+5t+2t^2-2t^3;
}"

echo "$idealkanm2_8" | ./contour --experimental alexander

#kanenobu.sh 8 -2 | ./contour --out alexander --noidealsimplify -Q --foxd 2 | ./contour alexander --experimental

#
# this should produce an integer overflow
#
idealkan8m2=`kanenobu.sh 10 -2 | ./contour --out alexander -Q --foxd 2`

idealkan8m2bis="
#
# --foxd 2
#
# WARNING: above max allowed int size, cannot complete simplification
# WARNING: inhibiting gcd computation due to integer size
# *** Warning: result can be noncanonical ***
ideal(t) {
+332177328;
+12+68240532t;
+1353876-1607472t;
+1+16057617t+55111t^2;
+69022+140468496t+544144t^2;
-477317+19746927t-12953t^2;
+5177-82979183t+709211t^2+616t^3;
-4513+652924t+26830052t^2+t^3;
+2589-58686413t+787863t^2+17104t^3;
}"

#echo "$idealkan8m2" | ./contour --experimental alexander

#
# questi due sono lo stesso ideale, la base di Groebner deve coincidere
#
ideala1="
ideal(t) {
+4;
+2t^2+t+1;
}"

ideala2="
ideal(t) {
+4;
+2t^2-t-1;
}"

#echo "$ideala1" | ./contour --experimental alexander
#echo "$ideala2" | ./contour --experimental alexander

ideal1="
ideal(t) {
+9;
+3t-3;
t^2+t+1;
}"

ideal2="
ideal(t) {
+6;
+3t-3;
t^2+t;
}"

ideal3="
ideal(t) {
+9;
3t^2+3t+3;
}"

ideal4="
ideal(t) {
+9t^3-9;
+3t^4-3t^3-3t+3;
+t^5+t^4+t^3-t^2-t-1;
}"

ideal4bis="
ideal(t) {
+9t^3-9;
0;
+3t^4+6t^3-3t-6;
+t^5+t^4+t^3-t^2-t-1;
+3t^4+6t^3-3t-6;
}"

echo "$ideal1" | ./contour alexander --experimental
echo "========="

echo "$ideal2" | ./contour alexander --experimental
echo "========="

echo "$ideal3" | ./contour alexander --experimental
echo "========="

echo "$ideal4" | ./contour alexander --experimental
echo "========="

echo "$ideal4bis" | ./contour alexander --experimental

