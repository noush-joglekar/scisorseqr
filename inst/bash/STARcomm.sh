#!/bin/bash
# Anoushka Joglekar
# April 2020

fastqDir=$1;
misc=$2;
starOut=$3;
progPath=$4;
refGenome=$5
numThreads=$6;

### For testing:
# /athena/tilgnerlab/store/hut2006/soft/src/star-mapper/2.5.2b/STAR-2.5.2b/bin/Linux_x86_64/STARlong
# /athena/tilgnerlab/store/hut2006/data/seq/genomes/M.musculus/mm10/wholeGenomeUnzipped/starIndex_gencode10_sequins/


## Create a fastq guide
for i in $(ls $fastqDir/*.fastq.gz) ; do echo $i | awk -v path=$(pwd) '{n=split($1,a,/\/|.fastq/); print a[n-1]"\t"path"/"$1}' ; done > $misc/fastqGuide

## Running STAR program over all files listed in the guide
n=`cat Misc/fastqGuide | wc -l`; \
for i in `seq 1 $n` ; do \
name=`cat Misc/fastqGuide | head -$i | tail -1 | awk '{print $1;}'` ; \
echo "### treating "$name >> $starOut/REPORT.star; \
file=`cat Misc/fastqGuide | head -$i | tail -1 | awk '{print $2}'` ; \
$progPath --readFilesCommand zcat --runMode alignReads \
--outSAMattributes NH HI NM MD --readNameSeparator space --outFilterMultimapScoreRange 1 \
--outFilterMismatchNmax 2000 --scoreGapNoncan -20 --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen -1 \
--scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 --alignEndsType Local \
--seedSearchStartLmax 50 --seedPerReadNmax 100000 --seedPerWindowNmax 1000 \
--alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 --runThreadN $numThreads \
--genomeDir $refGenome \
--readFilesIn $file >> $starOut/REPORT.star; \
mv Aligned.out.sam $starOut/${name}_Aligned.out.sam; \
mv Log.final.out $starOut/Log.final.out_${name}; \
mv Log.out $starOut/Log.out_${name}; \
mv Log.progress.out $starOut/Log.progress.out_${name}; \
mv SJ.out.tab $starOut/SJ.out.tab_${name}; \
done


for i in $(ls $starOut/SJ.out.tab*) ; do cat $i; done >> $starOut/SJ.out.all.tab

for i in $(ls $starOut/*.sam) ; \
do name=`echo $i | awk '{n=split($1,a,/\/|.sam/); print a[n-1]}'`; \
samtools view -bh $i -o $starOut/$name.bam ; \
rm $starOut/$name.sam; \
done


## Create a bam guide for further processing
for i in $(ls $starOut/*.bam) ; do echo $i | awk -v path=$(pwd) '{print path"/"$1}' ; done > $misc/bamGuide
