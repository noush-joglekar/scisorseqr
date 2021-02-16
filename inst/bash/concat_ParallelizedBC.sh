#!/bin/bash
## By Anoushka Joglekar 04/2020

f_name=$1;
outFolder=$2;

cat $outFolder/tmp* | tr -d "[ ']" > $outFolder/"DeconvBC_"$f_name".csv"
rm $outFolder/tmp*

cat $outFolder/"DeconvBC_"$f_name".csv" | awk '{tot++1;polyA[$2]++;totLen+=$11;\
if($5!~/no_bc|,/){bc++1;len+=$11;};if(!(see[$5])) {see[$5]++1;count++1;}; tso[$9]++1;}END \
{print "Total reads >= 200 bp:",tot; \
for(p in polyA) {print p":",polyA[p],"\nPercent",p":",polyA[p]*100/tot}; \
print "Unique Barcodes:", count-2;
print "Reads with barcode:",bc; \
print "Reads with barcode given polyA:", bc*100/polyA["poly_T_found"]; \
print "Percent barcoded reads:", bc*100/tot; \
print "Percent TSO found:",tso["TSO_found"]*100/tot; \
print "Percent doubleTSO:",tso["DoubleTSO"]*100/tot; \
print "Average barcoded read length:",len/bc; \
print "Average read length of file:", totLen/tot;}' > $outFolder/"DeconvBC_"$f_name"_summary"


num=$(ls $outFolder/"DeconvBC_"*"_summary" | wc -l);
if [ "$num" -gt 1 ]
then
awk -v num=$num '{split($0,t,": "); a[FNR]+=$NF; title[FNR]=t[1];}END \
{for(i=1;i<=FNR;i++)print title[i],":",a[i]/num;}' \
$outFolder/"DeconvBC_"*"_summary" > $outFolder/AllFiles_Summary
fi
