#!/bin/sh

#script to transform a 3DGB pdb file from the classical pipeline to structure with color and bonds

pdb_in=$1
pdb_out=$2


echo "generate modified structure with bonds"

cat -n "$pdb_in" | grep "C0" | awk '{ if ($4=="C01") print $2,$1,$4,$5,$6,"1",$8,$9,$10,$11,$12;else print $2,$1,$4,$5,$6,"2",$8,$9,$10,$11,$12}' > temp1.txt


awk '{printf("ATOM  %5d  %-4s%-3s %1s%4d    %8.3f%8.3f%8.3f %5.2f%6.2f\n", \
             $2, $3, $4, $5, $6, $7, $8, $9, $10, $11)}' temp1.txt | sed 's/,/./g' > temp2.txt

line=$(cat temp1.txt | grep "C01" | wc -l | awk '{print $1-1}')
echo "number of ATOMs"
echo "$line"

cat -n "$pdb_in" | grep "C01" | awk '{print "CONNECT",$1,$1+1}' |  head -"$line" > temp3.txt

line=$(cat temp1.txt | grep "C02" | wc -l | awk '{print $1-1}')
echo "number of ATOMs"
echo "$line"
cat -n "$pdb_in" | grep "C02" | awk '{print "CONNECT",$1,$1+1}' |  head -"$line" >> temp3.txt
awk '{printf("CONECT %4d %4d\n",$2,$3)}' temp3.txt >> temp2.txt
echo "END" >> temp2.txt


cat temp2.txt > "$pdb_out"
