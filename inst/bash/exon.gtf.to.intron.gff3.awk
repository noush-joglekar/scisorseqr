BEGIN{
    OFS="\t";
}

{
    ##########
    # 0. skipping unwanted data: 
    
    # 0a. must have the correct source
    if(source && toupper($2)!=toupper(source) ){next;}
       
    # 0b. must have the correct feature
    if(feature && toupper($3)!=toupper(feature)){next;}

    # 0c. the chromosome must be one of the legal chromosomes
    if(chromosomeList && !(toupper($1) in legalChromosomes)){next;}



    #####
    # 1. if we get here, we are interested in this exon

    # 1a. getting the transcript ID of this exon
    split($transcol,b,"\"");    
    trID=b[2];
    trIDwithChrom=b[2]"@"$1; 
    
    # 1b. exonID of this exon
    exIDwithChrom=$1"_"$4"_"$5"_"$7;

    # 1c. the ID of this chromosome-strand pair
    chromStrandID=$1"\t"$7;


    #####
    # 2. is this line consistent with previously found data ?

    # 2a. the transcript strand
    if(trIDwithChrom in transcriptStrand){	
	if($7 != transcriptStrand[trIDwithChrom]){
	    print "ERROR: unequal transcript strands:"$7" and "transcriptStrand[trIDwithChrom] > "/dev/stderr";
	    print $0 > "/dev/stderr";	    
	    exit(0);
	
	}
    }

    # 2b. remembering the strand of this transcript:
    transcriptStrand[trIDwithChrom]=$7;
    


    #####
    # 3. checking the ordering of the file:
    if(trIDwithChrom in lastLine){
	if($7=="+" && $4<=lastEnd[trIDwithChrom]){
	    print "sorting not in transcript direction:"> "/dev/stderr";
	    print lastLine[trIDwithChrom]> "/dev/stderr";
	    print "appeared before"> "/dev/stderr";
	    print $0 > "/dev/stderr";	    
	    exit(0);
	    
	}
	if($7=="-" && $5>=lastStart[trIDwithChrom]){
	    print "sorting not in transcript direction:"> "/dev/stderr";
	    print lastLine[trIDwithChrom]> "/dev/stderr";
	    print "appeared before"> "/dev/stderr";
	    print $0 > "/dev/stderr";	    
	    exit(0);
	    
	}   
    }




    #####
    # 4. if we get here, the line is consistent with previoulsy treated lines:  

    if(trIDwithChrom in lastLine){
	if($7=="+"){
	    start=lastEnd[trIDwithChrom]+1;
	    end=$4-1;       
	}
	else{
	    start=$5+1;
	    end=lastStart[trIDwithChrom]-1;
	    
	}
	if(!noID){
	    print $1,$2,"intron",start,end,".",$7,".","transcript_id_with_chr="trIDwithChrom;
	}
	else{
	    print $1,$2,"intron",start,end,".",$7,".","transcript_id_with_chr=.";	    
	}
    }
    

    


    #####
    # 5. remebering this exon as the last exon forthis transcript
    
    lastStart[trIDwithChrom]=$4;
    lastEnd[trIDwithChrom]=$5;
    lastLine[trIDwithChrom]=$0;
}
