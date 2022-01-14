#!/bin/bash
## By Anoushka Joglekar. Edited 04/2020


outDir=$1;
tabFile=$2;
isIso=$3;
isTSS=$4;
isPolyA=$5;

if [[ $isIso = "TRUE" ]]
then
	cat $outDir/$tabFile | awk 'BEGIN{OFS="\t"} NR>1 {a[$2"_"$4"_"$5]+=1;}END \
{for(i in a) {split(i,c,"_"); print c[1]"_"c[3],c[2],a[i]}}' \
| grep -v "none" > $outDir/"NumIsoPerCell"

	cat $outDir/$tabFile | awk '{a[$2"_"$5"~"$4]++;}END {for (isoBC in a) \
{split(isoBC,i,"~"); OFS="\t"; print i[1],i[2],a[isoBC]}}' | \
awk '!(a[$1]) {a[$1]=$0; next} \
a[$1] {w=$1; $1=""; a[w]=a[w]"\t"$2"\t"$3}END \
 {for(i in a) {print a[i]}}' FS="\t" OFS="\t" > $outDir/"IsoXNumInCell"

	cat $outDir/$tabFile | awk 'BEGIN{OFS="\t"} NR>1 {count[$2"_"$5"~"$3]++}END \
{for (gc in count) {split(gc,g,"~"); \
print g[1],g[2],count[gc];}}' > $outDir/"NumIsoPerCluster"

fi

if [[ $isTSS = "TRUE" ]]
then
	cat $outDir/$tabFile | awk 'BEGIN{OFS="\t"} NR>1 {if($6!="NA"){count[$2"_"$6"~"$3]++}}END \
{for (gc in count) {split(gc,g,"~"); \
print g[1],g[2],count[gc];}}' > $outDir/"NumTSSPerCluster"
fi

if [[ $isPolyA = "TRUE" ]]
then
	cat $outDir/$tabFile | awk 'BEGIN{OFS="\t"} NR>1 {if($7!="NA"){count[$2"_"$7"~"$3]++}}END \
{for (gc in count) {split(gc,g,"~"); \
print g[1],g[2],count[gc];}}' > $outDir/"NumPolyAPerCluster"
fi

