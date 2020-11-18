#!/bin/bash
# by Hagen Tilgner 5/2014
# Edited 07/2020 AJ

function checkForFile {
    file=$1;
    if [ ! -f $file ]
    then
       echo "ERROR:"$file" is not a regular file ... exiting "
       exit;
    fi
}

function checkForDirectory {
    dir=$1;
    if [ ! -d $dir ]
    then
       echo "ERROR:"$dir" is not a regular directory ... exiting "
       exit;
    fi
}

###############
# 0. the arguments
echo "+++++++++++++++++++++++++ 1. arguments";

echo "++++++++++++++++++ 1a. flexible";
unzipCommand=${1}
echo "unzipCommand="$unzipCommand

inputDir=${2}
outputDir=$inputDir
echo "inputDir=outputDir="$inputDir
checkForDirectory $inputDir

CAGEBED=${3}
echo "CAGEBED="$CAGEBED
checkForFile $CAGEBED

POLYABED=${4}
echo "POLYABED="$POLYABED
checkForFile $POLYABED

scriptDir=${5}
echo "scriptDir="$scriptDir
checkForDirectory $scriptDir

distance=${6}
echo "cpDistance="$distance

echo "+++++++++++++++++++++++++ 2. full-length analysis";
echo "++++++++++++++++++++++ 2a. CAGE";
echo "+++++++++++++++++++ 2a.1 execution";

awk -v exonMapping=$inputDir/mapping.bestperRead.RNAdirection.withConsensIntrons.gff.gz \
-v readsVsGenes=$inputDir/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz \
-v anno=$CAGEBED -v unzipCommand=$unzipCommand -v cageCutoff=1 -v verbose=2 \
-f ${scriptDir}/v2.1a.compareTranscriptStarts_vs_CageTSS.awk | gzip -c > \
$outputDir/mapping.bestperRead.RNAdirection.withConsensIntrons.cageRecordTSS.tab.gz

echo "+++++++++++++++++++ 2a.2 stats";
for i in 50 45 40 35 30 25 20 15 10 5; do n=`$unzipCommand \
$outputDir/mapping.bestperRead.RNAdirection.withConsensIntrons.cageRecordTSS.tab.gz \
| awk -v c=$i '$2<=c || $4<=c' | wc -l`; m=`$unzipCommand \
$outputDir/mapping.bestperRead.RNAdirection.withConsensIntrons.cageRecordTSS.tab.gz \
| wc -l`; o=`awk -v m=$m -v n=$n 'BEGIN{print n/m;}'`; \
echo -e $i"\t"$n"\t"$o; done

echo "++++++++++++++++++++++ 2b. polyA";
echo "+++++++++++++++++++ 2b.1 execution";
awk -v exonMapping=$inputDir/mapping.bestperRead.RNAdirection.withConsensIntrons.gff.gz \
-v readsVsGenes=$inputDir/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz \
-v anno=$POLYABED -v unzipCommand=$unzipCommand -v cageCutoff=0.0001 -v verbose=2 \
-f ${scriptDir}/v2.1a.compareTranscriptEnds_vs_PolyAs.awk | gzip -c > \
$outputDir/mapping.bestperRead.RNAdirection.withConsensIntrons.polyAsiteRecordPA.tab.gz

echo "+++++++++++++++++++ 2b.2 stats";
for i in 50 45 40 35 30 25 20 15 10 5; do n=`$unzipCommand \
$outputDir/mapping.bestperRead.RNAdirection.withConsensIntrons.polyAsiteRecordPA.tab.gz | \
awk -v c=$i '$2<=c || $4<=c' | wc -l`; m=`$unzipCommand \
$outputDir/mapping.bestperRead.RNAdirection.withConsensIntrons.polyAsiteRecordPA.tab.gz | \
wc -l`; o=`awk -v m=$m -v n=$n 'BEGIN{print n/m;}'`; echo -e $i"\t"$n"\t"$o; done


echo "+++++++++++++++++++++++++ 3. getting 5' and 3' complete files";
echo "++++++++++++++++++++++ 3.a stretches";

awk -v unzipCommand=$unzipCommand -v inputDir=$inputDir -v distance=$distance \
-v outputDir=$outputDir 'BEGIN{comm=unzipCommand" "\
outputDir"/mapping.bestperRead.RNAdirection.withConsensIntrons.cageRecordTSS.tab.gz"; \
while(comm|getline){d=1000; if(substr($1,1,1)=="#"){continue;} \
if($2<=distance){d=$2; TSS[$1]=$3; }  if($4<=distance && $4<d){TSS[$1]=$5;}}
comm=unzipCommand" "outputDir"/mapping.bestperRead.RNAdirection.withConsensIntrons.polyAsiteRecordPA.tab.gz"; \
while(comm|getline){d=1000; if(substr($1,1,1)=="#"){continue;} \
if($4<=distance){d=$4; polyA[$1]=$5;}  if($2<=distance && $2<d){polyA[$1]=$3;}} \
for(k in TSS){if(k in polyA){r2r[k]=k;}} \
comm=unzipCommand" "inputDir"/newIsoforms_vs_Anno_ignoreAnno/stretches.gz"; OFS="\t"; \
while(comm|getline){split($3,a,/\=|\@/); \
if(a[2] in r2r){print $1,$2,a[1]"="r2r[a[2]]"@"a[3],$4,TSS[r2r[a[2]]], polyA[r2r[a[2]]]; }}}' | \
gzip -c > $outputDir/newIsoforms_vs_Anno_ignoreAnno/CagePolyA.complete.stretches.gz

awk -v unzipCommand=$unzipCommand -v inputDir=$inputDir -v distance=$distance \
-v outputDir=$outputDir 'BEGIN{comm=unzipCommand" "\
outputDir"mapping.bestperRead.RNAdirection.withConsensIntrons.cageRecordTSS.tab.gz"; \
while(comm|getline) {d=1000; if(substr($1,1,1)=="#"){continue;} \
if($2<=distance){d=$2; TSS[$1]=$3; }  if($4<=distance && $4<d){TSS[$1]=$5;} \
if($2>distance && $4>distance){TSS[$1]="NoTSS"}} comm=unzipCommand" "\
outputDir"mapping.bestperRead.RNAdirection.withConsensIntrons.polyAsiteRecordPA.tab.gz"; \
while(comm|getline){d=1000; if(substr($1,1,1)=="#"){continue;} if($4<=distance){d=$4; polyA[$1]=$5;} \
if($2<=distance && $2<d){polyA[$1]=$3;} if($2>distance && $4>distance){polyA[$1]="NoPolyA"} } \
comm=unzipCommand" "inputDir"/newIsoforms_vs_Anno_ignoreAnno/stretches.gz"; OFS="\t"; \
while(comm|getline){split($3,a,/\=|\@/); {print $1,$2,$3,$4,TSS[a[2]], polyA[a[2]]; }} }' | \
gzip -c > $outputDir/newIsoforms_vs_Anno_ignoreAnno/incompleteStretches.gz

echo "++++++++++++++++++++++ 3.b transcriptWiseGenes";
$unzipCommand $inputDir/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz | \
awk -v unzipCommand=$unzipCommand -v outputDir=$outputDir \
'BEGIN{comm=unzipCommand" "outputDir"/newIsoforms_vs_Anno_ignoreAnno/CagePolyA.complete.stretches.gz"; \
while(comm|getline){split($3,a,/\=|\@/); use[a[2]]=1;}}{if($1 in use){print $0;}}' | gzip -c > \
$outputDir/CagePolyA.complete.mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz

echo "++++++++++++++++++++++ 3.c mapping file";
$unzipCommand $inputDir/mapping.bestperRead.RNAdirection.withConsensIntrons.gff.gz | \
awk -v unzipCommand=$unzipCommand -v outputDir=$outputDir \
'BEGIN{comm=unzipCommand" "outputDir"/newIsoforms_vs_Anno_ignoreAnno/CagePolyA.complete.stretches.gz"; \
while(comm|getline){split($3,a,/\=|\@/); use[a[2]]=1;}}{split($10,a,"\""); if(a[2] in use){print $0;}}' | gzip -c > \
$outputDir/CagePolyA.complete.mapping.bestperRead.RNAdirection.withConsensIntrons.gff.gz

echo "+++++++++++++++++++++++++ 3. done ";
exit;
