#!/bin/bash
## By Anoushka Joglekar 03/2020
## Edited 07/2020

rawOutput=$1;
filteredOutput=$2;
concat=$3;
filterReads=$4;
fqFolder=$5

if [[ $concat = "TRUE" ]]
then
cat $rawOutput/DeconvBC*.csv | grep -v "no_clust" | grep -v "no_bc" | grep -v "," \
>> $filteredOutput/"FilteredDeconvBC_AllFiles.csv"
sed -i $'1i\\\nName\tT9_Status\tStrand\tT9_position\tBarcodes\tBarcode_position\tCluster\tUMI\tTSO_Status\tTSO_position\tLength' \
$filteredOutput/"FilteredDeconvBC_AllFiles.csv"
if [[ $filterReads = "TRUE" ]]
then
for file in $rawOutput/DeconvBC*.csv
do
       	f_name=`echo $file | awk '{n=split($1,a,/\/|DeconvBC_|.csv/); print a[n-1]}'`
       	cat $file | grep -v "no_clust" | grep -v "no_bc" | grep -v "," | cut -f 1 > x;
       	zcat $fqFolder/$f_name".fastq.gz" | awk 'NR==FNR {a[$1]; next} $1 in a \
{print; getline; print; getline; print; getline; print}' x - | gzip >> $filteredOutput/"Barcoded_AllFiles.fastq.gz";
       	rm x;
done
fi

else
for file in $rawOutput/DeconvBC*.csv
do
	f_name=`echo $file | awk '{n=split($1,a,/\/|DeconvBC_|.csv/); print a[n-1]}'`
	cat $file | grep -v "no_clust" | grep -v "no_bc" | grep -v "," > \
	$filteredOutput/"FilteredDeconvBC_"$f_name".csv"
	sed -i $'1i\\\nName\tT9_Status\tStrand\tT9_position\tBarcodes\tBarcode_position\tCluster\tUMI\tTSO_Status\tTSO_position\tLength' \
	$filteredOutput/"FilteredDeconvBC_"$f_name".csv"
	if [[ $filterReads = "TRUE" ]]
  then
       	cat $file | grep -v "no_clust" | grep -v "no_bc" | grep -v "," | cut -f 1 > x;
       	zcat $fqFolder/$f_name".fastq.gz" | awk 'NR==FNR {a[$1]; next} $1 in a \
{print; getline; print; getline; print; getline; print}' x - | gzip > $filteredOutput/"Barcoded_"$f_name".fastq.gz";
       	rm x;
  fi
done
fi

