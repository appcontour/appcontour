#!/bin/bash
#
#testname:groebner basis for Laurent in one indeterminate
#

giove="ideal(t) {
+4t-3;
+2t+3;
}"

echo "Canonifying ideal:"
echo "$giove"
echo "gives"
echo "$giove" | ./contour --experimental alexander
echo "=================="
#
# this is the ideal (no simplification) from kanenobu -2 8
#
idealkanm2_8="ideal(t) {
+4-7t+4t^2;
+4-6t+t^2+t^3;
+4-6t+t^2+t^3;
-4+5t+2t^2-2t^3;
}"

echo "Canonifying ideal:"
echo "$idealkanm2_8"
echo "gives"
echo "$idealkanm2_8" | ./contour --experimental alexander
echo "=================="

#kanenobu.sh 8 -2 | ./contour --out alexander --noidealsimplify -Q --foxd 2 | ./contour alexander --experimental

#
# this should NOT produce an integer overflow
#
idealkan8m2=`kanenobu.sh 10 -2 | ./contour --out alexander -Q --foxd 2`
echo "Canonifying ideal obtained from kanenoby 10 -2"
echo "gives"
echo "=================="

idealkan8m2bis="ideal(t) {
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

echo "Canonify ideal:"
echo "$idealkan8m2bis"
echo "gives"
echo "$idealkan8m2" | ./contour --experimental alexander
echo "=================="

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

echo "The next two should coincide"
echo "$ideala1" | ./contour --experimental alexander
echo "=================="
echo "$ideala2" | ./contour --experimental alexander
echo "=================="

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

echo "Result of some ideals, see source for details:"
echo "$ideal1" | ./contour alexander --experimental
echo "=================="

echo "$ideal2" | ./contour alexander --experimental
echo "=================="

echo "$ideal3" | ./contour alexander --experimental
echo "=================="

echo "$ideal4" | ./contour alexander --experimental
echo "=================="

echo "$ideal4bis" | ./contour alexander --experimental
echo "=================="

