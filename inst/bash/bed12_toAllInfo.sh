#!/bin/sh
# By Anoushka Joglekar 06.2022
# File to convert from IsoQuant corrected BED file to
# a corrected all-info file
# Column 6: Corrected Intron Chain
# Column 9: Corrected Exon Chain
# Column 10: Old intron chain
# Column 11: Old exon chain

corrBed=$1;
allInfo=$2;
allInfoCorrected=$3;

cat $corrBed | awk 'NR>1 {exonChain=""; rn=$4; chr=$1; str=$6; nE=$10; nz=split($11,bsz,","); nt=split($12,bst,","); \
if(nE>=2) {for(i=1;i<=nE;i++) {exon=chr"_"1+$2+bst[i]"_"$2+bst[i]+bsz[i]"_"str; \
exonChain=exonChain";%;"exon} print rn"\t"exonChain} }' > corrBed_ai

awk -v allInfo=$allInfo 'BEGIN{comm="cat corrBed_ai"; while(comm|getline) {ec[$1]=$2;} comm="zcat "allInfo; \
while(comm|getline) {if($1 in ec) {$13=$11; $12=$10; $11=$9; $10=$6; $9=ec[$1]; split($11,o,/_|;%;/); \
split($9,n,/_|;%;/); if(n[5]!=o[5]) {gsub(n[5],o[5],$9);} {if($9==$11) {$6=$10;} else {intronChain=""; \
nE=split($9,exChain,";%;"); for(i=2;i<=(nE-1);i++){split(exChain[i],start,"_"); split(exChain[i+1],end,"_"); \
intron=start[1]"_"start[3]+1"_"end[2]-1"_"end[4]; intronChain=intronChain";%;"intron; $6=intronChain}} } OFS="\t"; \
print }}}' | gzip -c > $allInfoCorrected

rm corrBed_ai

echo "Done!"
