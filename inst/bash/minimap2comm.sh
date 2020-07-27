#!/bin/bash
# Anoushka Joglekar
# April 2020

fastqDir=$1;
mmOut=$2;
progPath=$3;
annotation=$4;
numThreads=$5;

### For testing:
# /athena/tilgnerlab/store/hut2006/data/annotations/M.Musculus/mm10/gencode.vM21.annotation.gtf.gz

## Create a fastq guide
for i in $(ls $fastqDir/*.fastq.gz) ; do echo $i | awk -v path=$(pwd) \
'{n=split($1,a,/\/|.fastq/); print a[n-1]"\t"path"/"$1}' ; done > Misc/fastqGuide

## Running minimap2 over all files listed in the guide
n=`cat Misc/fastqGuide | wc -l` ; for i in `seq 1 $n` ; do name=`cat Misc/fastqGuide \
| head -$i | tail -1 | awk '{print $1;}'` ; \
echo "### treating "$name >> $mmOut/REPORT.minimap ; \
file=`cat Misc/fastqGuide | head -$i | tail -1 | awk '{print $2}'` ; \
$progPath -t $numThreads -ax splice --secondary=no $annotation \
$file > $mmOut/$name.sam ; done &>> $mmOut/REPORT.minimap2

for i in $(ls $mmOut/*.sam) ; \
do name=`echo $i | awk '{n=split($1,a,/\/|.sam/); print a[n-1]}'`; \
samtools view -bh $i -o $mmOut/$name.bam ; \
rm $mmOut/$name.sam; \
done


## Create a bam guide for further processing
for i in $(ls $mmOut/*.bam) ; do echo $i | awk -v path=$(pwd) '{print path"/"$1}' ; done > Misc/bamGuide
