This directory contains the Dec. 2011 (GRCm38/mm10) assembly of the mouse genome
(mm10, Genome Reference Consortium Mouse Build 38 (GCA_000001635.2)) in one gzip-compressed FASTA file per chromosome.

Repeats from RepeatMasker and Tandem Repeats Finder (with period
of 12 or less) are shown in lower case; non-repeating sequence is
shown in upper case.

This assembly was produced by the Mouse Genome Sequencing Consortium,
and the National Center for Biotechnology Information (NCBI).
For more information on the mouse genome, see the project website:

See also: http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/mouse/
          http://www.ncbi.nlm.nih.gov/genome/52

Files included in this directory:

  - chr*.fa.gz: compressed FASTA sequence of each chromosome.
    Each chromosome is in a separate file in a gzip Fasta format.
    Repeats -- which are shown in lower case -- are annotated by
    RepeatMasker run at the sensitive setting and Tandem Repeats Finder
    (repeats of period 12 or less).

  - md5sum.txt - MD5 checksum of these files to verify correct transmission

The main assembly is contained in the chrN.fa.gz files, where N is the name
of the chromosome.  The chrN_random.fa.gz files contain clones that are not
yet finished or cannot be placed with certainty at a specific place on
the chromosome.  The chrUn_random.fa.gz sequence are unfinished clones,
or clones that can not be tenatively placed on any chromosome.

------------------------------------------------------------------
If you plan to download a large file or multiple files from this 
directory, we recommend that you use ftp rather than downloading the 
files via our website. To do so, ftp to hgdownload.cse.ucsc.edu, then 
go to the directory goldenPath/mm10/chromosomes. To download multiple 
files, use the "mget" command:

    mget <filename1> <filename2> ...
    - or -
    mget -a (to download all the files in the directory)

Alternate methods to ftp access.
    
Using an rsync command to download the entire directory:
    rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/ .
For a single file, e.g. chrM.fa.gz
    rsync -avzP 
        rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chrM.fa.gz .
    
Or with wget, all files:
    wget --timestamping 
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/*'
With wget, a single file:
    wget --timestamping 
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chrM.fa.gz' 
        -O chrM.fa.gz
    
To uncompress the fa.gz files:
    gunzip <file>.fa.gz


All the files in this directory are freely available for public use.

------------------------------------------------------------------
This file last updated: 2012-02-09 - 09 February 2012
------------------------------------------------------------------
