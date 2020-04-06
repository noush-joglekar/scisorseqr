# by Hagen Tilgner for mapping of pacBio reads : Jan 1st, 2011
# expects -v intronTypeFile=??? (teh file that gives the splice site consensus) -v bestMatchFile=??? (the mapping file) -v outputCorrectFile=??? -v outputUnClassifiableFile=??? 
function reverseStrand(x){
    n=split(x,a,"\t"); 
    if(a[7]!="+" && a[7]!="-"){print "STRAND-ERROR:"x > "/dev/stderr";} 
    if(a[7]=="+"){strand="-";}else{strand="+";}  
    s=a[1]"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[5]"\t"a[6]"\t"strand; 
    for(i=8;i<=NF;i++){s=s"\t"$i;} 
    return(s);     
} 

function getID(x){
    split($9,a,"transcript_id_with_chr="); 
    split(a[2],b,"@"); 
    return(b[1]);
} 

function getType(x){
    split(x,a,"="); 
    return(a[2]); 
} 


BEGIN{
    comm="zcat "intronTypeFile; 
    while(comm | getline){
	trID=getID($9); 
	if(!(trID in cStrand)){cStrand[trID]=0;} 
	if(!(trID in wStrand)){wStrand[trID]=0;} 
	if(!(trID in notClassifiable)){notClassifiable[trID]=0;}    
	t=getType($10);   
	if(t=="GTAG" || t=="GCAG" || t=="ATAC"){cStrand[trID]++; continue;}  
	if(t=="CTAC" || t=="CTGC" || t=="GTAT"){wStrand[trID]++; continue;}      
	notClassifiable[trID]++;
    }     

    while(getline<bestMatchFile>0){
	split($10,a,"\"");  
	readID=a[2];   
	if(readID in wStrand){
	    if(cStrand[readID]>0 && (wStrand[readID]==0 && notClassifiable[readID]==0)){print $0 > outputCorrectFile ;continue;}  
	    if(wStrand[readID]>0 && (cStrand[readID]==0 && notClassifiable[readID]==0)){print reverseStrand($0) > outputCorrectFile ;continue;}
	    if(notClassifiable[readID]>0){print $0 > outputUnClassifiableFile ;continue;}
	}   
    }
}