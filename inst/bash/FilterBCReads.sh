#!/bin/bash

rawOutput=$1;
filteredOutput=$2;
concat=$3;

if [[ $concat = "TRUE" ]]
then
cat $rawOutput/DeconvBC*.csv | grep -v "no_clust" | grep -v "no_bc" | grep -v "," >> $filteredOutput/"FilteredDeconvBC_AllFiles.csv"
sed -i $'1i\\\nName\tT9_Status\tStrand\tT9_position\tBarcodes\tBarcode_position\tCluster\tUMI\tTSO_Status\tTSO_position\tLength' $filteredOutput/"FilteredDeconvBC_AllFiles.csv"
else
for file in $rawOutput/DeconvBC*.csv
do
	f_name=`echo $file | awk '{n=split($1,a,/\/|DeconvBC_|.csv/); print a[n-1]}'`
	cat $file | grep -v "no_clust" | grep -v "no_bc" | grep -v "," > $filteredOutput/"FilteredDeconvBC_"$f_name".csv"
	sed -i $'1i\\\nName\tT9_Status\tStrand\tT9_position\tBarcodes\tBarcode_position\tCluster\tUMI\tTSO_Status\tTSO_position\tLength' $filteredOutput/"FilteredDeconvBC_"$f_name".csv"
done
fi

