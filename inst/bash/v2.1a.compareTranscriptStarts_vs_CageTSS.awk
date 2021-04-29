# by Hagen;

function max(a,b){
    if(a>b){return(a);}
    return(b);
}
function min(a,b){
    if(a<b){return(a);}
    return(b);
}
function allOptionsThere(){  
    if(!exonMapping){
	print "ERROR: -v exonMapping=<string> must be provided; exiting" > "/dev/stderr";
	exit;
    }
    if(!readsVsGenes){
	print "ERROR: -v readsVsGenes=<string> must be provided; exiting" > "/dev/stderr";
	exit;
    }
    if(!anno){
	print "ERROR: -v anno=<string> must be provided; exiting" > "/dev/stderr";
	exit;
    }
    # from when on do we actually use
    if(!(cageCutoff)){
	print "ERROR: -v cageCutoff must be provided;  exiting" > "/dev/stderr";
	exit;	
    }
    if(!(unzipCommand)){
	print "ERROR: -v unzipCommand must be provided;  exiting" > "/dev/stderr";
	exit;	
    }

}

function comp(x){
    x=toupper(x) ;
    if(x=="A"){return("T");}
    if(x=="T"){return("A");} 
    if(x=="C"){return("G");} 
    if(x=="G"){return("C");} 
    if(x=="N"){return("N");}   
    
    print "ERROR:"x" is not a vlaid symbol" > "/dev/stderr";
    exit(0);
} 

function revcomp(x){
    y=""; 
    for(j=length(x);j>=1;j--){
	y=y""comp(substr(x,j,1));
    } 
    return(y);
}

function abs(x){
    if(x<0){x=(-1)*x;} 
    return(x);
}


BEGIN{

    # need to have the required options
    allOptionsThere();
    OFS="\t";
    if(verbose>=1)
	print "# 1. reading annotation" > "/dev/stderr";

    comm=unzipCommand" "anno;
    while(comm|getline){
	
	if($5>=cageCutoff){
	    xP1=$2+1;
	    for(i=$2+1;i<=$3;i++)
		TSS[$1"_"i"_"$6]=$1"_"xP1"_"$3"_"$6;	
	}
    } 

    if(verbose>=1)
	print "# 2. reading readsVsGenes" > "/dev/stderr";
  
    comm=unzipCommand" "readsVsGenes; 
    while(comm | getline){
	if($2=="none" || $3~/problem/ || $4~/problem/){skippedReads[$1]=1; continue;}
	read2Gene[$1]="\""$2"\";";
    }

    if(verbose>=1)
	print "# 3. reading exonMapping" > "/dev/stderr";

    comm=unzipCommand" "exonMapping;
    while(comm|getline){
	split($10,a,"\""); 
	r=a[2]; 
	if(!(r in read2Gene)){if(!(r in skippedReads)){print "ERROR"r > "/dev/stderr"; exit;}else{continue;}}    
	if(r in read2ReadEnd){
	    if($7=="+"){read2ReadEnd[r]=min(read2ReadEnd[r],$4);}
	    if($7=="-"){read2ReadEnd[r]=max(read2ReadEnd[r],$5);}   
	}
	else{
	    if($7=="+"){read2ReadEnd[r]=$4;}
	    if($7=="-"){read2ReadEnd[r]=$5;}
	    readStrand[r]=$7;
	    readChrom[r]=$1;
	} 
    }

    if(verbose>=1)
	print "# 4. finding the closest TSS" > "/dev/stderr";
    
    readCounter=0;
    print "#readID\tminDUpstream\tminDUpstreamTSS\tminDDownstream\tminDDownstreamTSS\tGene\treadChrom\treadStrand"; 
    for(r in read2ReadEnd){
	readCounter++;
	if(readCounter % 1000000 == 0 && verbose>=2){
	    print readCounter > "/dev/stderr";
	}
	if(!(r in read2Gene && r in readStrand)){print "ERROR2" > "/dev/stderr"; exit;}
	minDDownstream=1000000000;
	minDDownstreamTSS="none--";

	minDUpstream=1000000000;
	minDUpstreamTSS="none--";

	TSS["none--"]="none--";
	#print r"\t"read2ReadEnd[r];
	#for(s in TSS){
	#    print "TSS\t"s;
	#}
	
	for(i=read2ReadEnd[r]-50;i<=read2ReadEnd[r]+50;i++){
	    s=readChrom[r]"_"i"_"readStrand[r];
	    if(s in TSS){
		
		# closest upstream
		if(readStrand[r]=="+" && i<=read2ReadEnd[r] && abs(read2ReadEnd[r]-i)<minDUpstream){
		    minDUpstream=abs(read2ReadEnd[r]-i);
		    minDUpstreamTSS=s;			     
		}
		if(readStrand[r]=="-" && i>=read2ReadEnd[r] && abs(read2ReadEnd[r]-i)<minDUpstream){
		    minDUpstream=abs(read2ReadEnd[r]-i);
		    minDUpstreamTSS=s;			     
		}
		# closest downstream
		if(readStrand[r]=="+" && i>=read2ReadEnd[r] && abs(read2ReadEnd[r]-i)<minDDownstream){
		    minDDownstream=abs(read2ReadEnd[r]-i);
		    minDDownstreamTSS=s;	
		}		
		if(readStrand[r]=="-" && i<=read2ReadEnd[r] && abs(read2ReadEnd[r]-i)<minDDownstream){
		    minDDownstream=abs(read2ReadEnd[r]-i);
		    minDDownstreamTSS=s;	
		}

	    }	    
	}
	print r"\t"minDUpstream"\t"TSS[minDUpstreamTSS]"\t"minDDownstream"\t"TSS[minDDownstreamTSS]"\t"read2Gene[r]"\t"readChrom[r]"\t"readStrand[r];          
    }
}

 
