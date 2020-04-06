# by Hagen Tilgner for mapping of pacBio reads : Jan 8th, 2011
# expects to get 
# - an annotation file sorted by position piped in
# -v sortedAnno=
# -v intronGFF=
# -v feature=
# -v transIdColumn=
# -v task=



function d(x,y){
    split(x,X,"_"); 
    split(y,Y,"_");

    if(X[1]!=Y[1]){
	print "ERROR:"X[1]"!="Y[1] > "/dev/stderr";
	exit(0);
    }
    if(X[4]!=Y[4]){
	print "ERROR:"X[4]"!="Y[4] > "/dev/stderr";
	exit(0);
    }    

    if(X[3]<Y[2]){
	return(Y[2]-X[3]-1);
    }
    if(Y[3]<X[2]){
	return(X[2]-Y[3]-1);
    }
    return(-1);
      
} 

# treating the annotation
BEGIN{

    print "### executing pB.getGenes.awk" > "/dev/stderr";
    if(!sortedAnno){print "ERROR: no value for sortedAnno" > "/dev/stderr";exit(0);}
    if(!intronGFF){print "ERROR: no value for intronGFF" > "/dev/stderr";exit(0);}
    if(!feature){print "ERROR: no value for feature" > "/dev/stderr";exit(0);}
    if(!transIdColumn){print "ERROR: no value for transIdColumn" > "/dev/stderr";exit(0);}
    if(!geneIdColumn){print "ERROR: no value for geneIdColumn" > "/dev/stderr";exit(0);}
    

    print "## A. parsing file2="intronGFF > "/dev/stderr"
    comm="cat "intronGFF;
    while(comm | getline){	

	split($9,a,"transcript_id_with_chr=");
	split(a[2],b,"@");
	readID=b[1];	
	read2Gene[readID]="";

	l1=$4-1;
	key1=$1"_"l1"_"$7;
	
	if(key1 in endSite2Read){endSite2Read[key1]=endSite2Read[key1]";"readID}
	else{endSite2Read[key1]=readID;}
       
	l2=$5+1;
	key2=$1"_"l2"_"$7;
	
	if(key2 in startSite2Read){startSite2Read[key2]=startSite2Read[key2]";"readID}
	else{startSite2Read[key2]=readID;}

    }


    


    print "## B. parsing annotation: " > "/dev/stderr"
    comm="cat "sortedAnno;
    while(comm | getline){

	if(feature && $3!=feature){continue;}
        if(lastColumn4 && $4<lastColumn4){
            print "ERROR:lastColumn4="lastColumn4"; and $0="$0 > "/dv/stderr";
            exit(0);
        }
	if($7!="+" && $7!="-"){
	    print "ERROR: cannot deal with strand="$7 > "/dev/stderr";
	    exit(0);
	}
	if($transIdColumn in strand && strand[$transIdColumn]!=$7){
	    print "ERROR: strands do no match:"strand[$transIdColumn]"!="$7 > "/dev/stderr";
	    exit(0);
	}
	if($transIdColumn in chr && chr[$transIdColumn]!=$1){
	    print "ERROR: chroms do no match:"chr[$transIdColumn]"!="$1 > "/dev/stderr";
	    exit(0);
	}
	split($geneIdColumn,G,"\"");
	gene[$transIdColumn]=G[2];
	n[$transIdColumn]++;
	exon[$transIdColumn"\t"n[$transIdColumn]]=$1"_"$4"_"$5"_"$7;	
	strand[$transIdColumn]=$7;
	chr[$transIdColumn]=$1;
	lastColumn4=$4;
    }
    
    print "## C. going over all annotated transcripts: " > "/dev/stderr";
    for(tr in strand){
    	
	# going over all exons of this transcript
	for(i=1;i<=n[tr];i++){
	    
	    split(exon[tr"\t"i],a,"_");
	    keyStart=a[1]"_"a[2]"_"a[4];
	    keyEnd=a[1]"_"a[3]"_"a[4];
	    

	    if(keyEnd in endSite2Read && i<n[tr] && !(gene[tr]"\t"keyEnd"\t""end" in geneSpliceSitePair) ){				
		m=split(endSite2Read[keyEnd],b,";");
		for(j=1;j<=m;j++){
		    read2Gene[b[j]]=read2Gene[b[j]]";"gene[tr];
		}
		#sitesEnd[keyEnd"\t"gene[tr]]=1;
	    }
	    if(i<n[tr]){
		geneSpliceSitePair[gene[tr]"\t"keyEnd"\t""end"]=1;
	    }
	    if(keyStart in startSite2Read && i>1 && !(gene[tr]"\t"keyStart"\t""start" in geneSpliceSitePair) ){
		m=split(startSite2Read[keyStart],b,";");
		for(j=1;j<=m;j++){
		    read2Gene[b[j]]=read2Gene[b[j]]";"gene[tr];
		}
		#sitesStart[keyStart"\t"gene[tr]]=1;
	    }
	    if(i>1){
		geneSpliceSitePair[gene[tr]"\t"keyStart"\t""start"]=1;		
	    }
       			    
	}	    
    }
    print "## D. counting the number of splice sites per gene:" > "/dev/stderr";
    for(k in geneSpliceSitePair){
	x=split(k,a,"\t");
	if(x!=3){
	    print "ERROR: "k" is not a valid key in geneSpliceSitePair" > "/dev/stderr";
	    exit(0);
	}
	spliceSiteNumber[a[1]]++;

	if(a[2]"\t"a[3] in spliceSite2Gene){
	    problematicGene[a[1]]=""; #print a[1] > problematicGeneOutFile;
	    problematicGene[spliceSite2Gene[a[2]"\t"a[3]]]=""; #print spliceSite2Gene[a[2]"\t"a[3]] > problematicGeneOutFile;
	}
	else{
	    spliceSite2Gene[a[2]"\t"a[3]]=a[1];
	}
    }
    

    print "## E. checking whether a read has equal numbers of splice sites with multiple genes: " > "/dev/stderr"    
    for(r in read2Gene){
	#if(r!="m121212_085615_00126_c100418772550000001523036412191280_s1_p0/47403/ccs.path1"){
	#    continue;
	#}
	for(k in h){
	    #print "deleting "k" from h"; 
	    delete h[k];
	}
	m=split(read2Gene[r],b,";");
	#print "m="m;
	for(i=2;i<=m;i++){
	    #print "i="i"; b[i]="b[i]
	    h[b[i]]++;
	}	
	g="none";
	v=-1;
	problemRead="fineRead";
	problemGene="fineGene";
	for(k in h){
	    if(h[k]>v){
		problemRead="fineRead";
		if(k in problematicGene){
		    problemGene="problematicGene";
		}
		else{
		    problemGene="fineGene";
		}
		g=k;
		v=h[k];
		continue;
	    }
	    if(h[k]==v){
		problemRead="problematicRead";
		if(k in problematicGene){
		    problemGene="problematicGene";
		}
		g=g"@;@"k;
		v=h[k];
	    }
	}	
	print r"\t"g"\t"problemRead"\t"problemGene;
    }



}


