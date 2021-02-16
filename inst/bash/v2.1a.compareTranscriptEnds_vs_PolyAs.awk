# goal: remove polyA-tails and remove revcomp sequence to RNA if necessary
# by Hagen;

# exampel of ccs read data: 
#nuvol:~ htilgner$ gzcat /Users/htilgner/data/trios/input/pacBioIsoforms/v1/yorubian/Gm19238/round1/ccs/fa/m120413_003026_00126_c100277812550000001523007807041240_s2_p0.ccs.fasta.gz | head -2
#>m120413_003026_00126_c100277812550000001523007807041240_s2_p0/57/ccs
#TTTTTTTTGAAGGTTCTCAGGTCTTTATTTGCTCTCTCAACTTCCAGGAATTGACTTATTTAATTAATCC

# exampel of subread data
#nuvol:~ htilgner$ gzcat /Users/htilgner/data/trios/input/pacBioIsoforms/v1/yorubian/Gm19238/round1/subreads/fa/GM19238_poly-AcDNA_opt_smrtanalysis_common_jobs_016_016525_data_filtered_subreads.fasta.gz| head -2
#>m120527_035935_00126_c100333442550000001523018909161283_s1_p0/7/0_334
#TTTGCGACCTTTCGCCACGACGTCAGCGCGCTTTCCCGGGCAGACGCCTC


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
		polyA[$1"_"i"_"$6]=$1"_"xP1"_"$3"_"$6;
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
	    if($7=="+"){read2ReadEnd[r]=max(read2ReadEnd[r],$5);}
	    if($7=="-"){read2ReadEnd[r]=min(read2ReadEnd[r],$4);}   
	}
	else{
	    if($7=="+"){read2ReadEnd[r]=$5;}
	    if($7=="-"){read2ReadEnd[r]=$4;}
	    readStrand[r]=$7;
	    readChrom[r]=$1;
	} 
    }

    if(verbose>=1)
	print "# 4. finding the closest polyA" > "/dev/stderr";
    
    readCounter=0;
    print "#readID\tminDUpstream\tminDUpstreamPolyA\tminDDownstream\tminDDownstreamPolyA\tGene\treadChrom\treadStrand"; 
    for(r in read2ReadEnd){
	readCounter++;
	if(readCounter % 1000000 == 0 && verbose>=2){
	    print readCounter > "/dev/stderr";
	}
	if(!(r in read2Gene && r in readStrand)){print "ERROR2" > "/dev/stderr"; exit;}
	minDDownstream=1000000000;
	minDDownstreampolyA="none--";

	minDUpstream=1000000000;
	minDUpstreampolyA="none--";
	
	polyA["none--"]="none--";
        #print r"\t"read2ReadEnd[r];
	#for(s in polyA){
	#    print "polyA\t"s;
	#}
	
	for(i=read2ReadEnd[r]-50;i<=read2ReadEnd[r]+50;i++){
	    s=readChrom[r]"_"i"_"readStrand[r];
	    if(s in polyA){
		
		# closest upstream
		if(readStrand[r]=="+" && i<=read2ReadEnd[r] && abs(read2ReadEnd[r]-i)<minDUpstream){
		    minDUpstream=abs(read2ReadEnd[r]-i);
		    minDUpstreampolyA=s;			     
		}
		if(readStrand[r]=="-" && i>=read2ReadEnd[r] && abs(read2ReadEnd[r]-i)<minDUpstream){
		    minDUpstream=abs(read2ReadEnd[r]-i);
		    minDUpstreampolyA=s;			     
		}
		# closest downstream
		if(readStrand[r]=="+" && i>=read2ReadEnd[r] && abs(read2ReadEnd[r]-i)<minDDownstream){
		    minDDownstream=abs(read2ReadEnd[r]-i);
		    minDDownstreampolyA=s;	
		}		
		if(readStrand[r]=="-" && i<=read2ReadEnd[r] && abs(read2ReadEnd[r]-i)<minDDownstream){
		    minDDownstream=abs(read2ReadEnd[r]-i);
		    minDDownstreampolyA=s;	
		}

	    }	    
	}
	print r"\t"minDUpstream"\t"polyA[minDUpstreampolyA]"\t"minDDownstream"\t"polyA[minDDownstreampolyA]"\t"read2Gene[r]"\t"readChrom[r]"\t"readStrand[r];          
    }
}

 
