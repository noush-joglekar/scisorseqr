function comp(x){
    x=toupper(x) ;
    if(x=="A"){return("T");}
    if(x=="T"){return("A");} 
    if(x=="C"){return("G");} 
    if(x=="G"){return("C");} 
    if(x=="N"){return("N");}   
} 
function revcomp(x){
    y=""; 
    for(i=length(x);i>=1;i--){y=y""comp(substr(x,i,1));} 
    return(y);
}

{seq=$2;}


END{
    while(getline<file>0){
	if($1!=chr){continue;} 
	key=$1"_"$4"_"$5"_"$7; 
	if(key in c){ 
	    print $0"\tsite_consensus="c[key]"\tsite_consensus_ext="c2[key]; 
	    continue;
	} 
	compare=toupper(substr(seq,$4,2))""toupper(substr(seq,$5-1,2)); 
	compare2=toupper(substr(seq,$4-5,2+10))"NNNNN"toupper(substr(seq,$5-1-5,2+10));  
	if($7=="-"){
	    compare=revcomp(compare);
	} 
	print $0"\tsite_consensus="compare"\tsite_consensus_ext="compare2; 
	c[key]=compare; 
	c2[key]=compare2;   
    }
}