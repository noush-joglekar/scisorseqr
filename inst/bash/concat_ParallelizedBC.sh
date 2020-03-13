#!/bin/bash

f_name=$1;
outFolder=$2;

cat $outFolder/tmp* | tr -d "[']" > $outFolder/"DeconvBC_"$f_name".csv"
rm $outFolder/tmp*

cat $outFolder/"DeconvBC_"$f_name".csv" | awk '{tot++1;polyA[$2]++;totLen+=$11;\
if($5!~/no_bc|,/){bc++1;len+=$11;};if(!(see[$5])) {see[$5]++1;count++1;}; tso[$9]++1;}END \
{print "Total reads >= 200 bp:",tot; \
for(p in polyA) {print p":",polyA[p],"\nPercent",p":",polyA[p]/tot}; \
print "Unique Barcodes:", count-1;
print "Reads with barcode:",bc; \
print "Reads with barcode given polyA:", bc/polyA["poly_T_found"]; \
print "Percent barcoded reads:", bc/tot; \
print "Percent TSO found:",tso["TSO_found"]/tot; \
print "Percent doubleTSO:",tso["DoubleTSO"]/tot; \
print "Average barcoded read length:",len/bc; \
print "Average read length of file:", totLen/tot;}' > $outFolder/"DeconvBC_"$f_name"_summary"
