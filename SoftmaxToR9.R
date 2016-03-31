#–––––––––––––––––––03/02/2011 by Quentin GEISSMANN q.geissmann@sheffield.ac.uk–––––––––––––––––––#
# Modified by Joe Gallagher (joedgallagher@gmail.com)
#
# V9(02/02/2014) fixes: enabled cutting and reading of Softmax files that have been prematurely stopped
# V8(16/09/2013) fixes: enabled reading of 'endpoint' (single read) Softmax files
# V7(29/02/2012) fixes: enabled reading of files that don't have a "Time(xx,xx,xx)" in the end
# V6(29/05/2011) fixes: NA for more than 48h records.
#
# This script transforms a Softmax, V-5.1 ouput file to a dataframe
# The dataframe is then save to "RdataFrame-FileName.csv"
# Look to the end of this file for examples

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
library(chron)

fun<-function(a){
	time<-60*60*as.numeric(unlist(strsplit(a,"\\:")))[1]+60*as.numeric(unlist(strsplit(a,"\\:")))[2]+as.numeric(unlist(strsplit(a,"\\:")))[3]
	return(time)
	}
		
		
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
sofToR<-function(file,nr=8,nc=12,rnames=1:nr,cnames=1:nc){
#–––––––––––––––––––––––––––––––––––
#the function returns a dataframe
#–––––––––––––––––––––––––––––––––––
#file = file from Softmax	(string)
#nr = number of rows		default = 8	(scalar)
#nr = number of columns 	default = 12	(scalar)
#rnames = names of rows		default = 1:8	(vector of string,factors,integer... | length=nr)
#cnames = names of columns	default = 1:12	(vector of string,factors,integer... | length=nc)
#–––––––––––––––––––––––––––––––––––
	#according to the exportation format, skips 5 or 3 lines
	skp<-ifelse(scan(file,what="character",skip=2,nmax=1)=="~End",5,3)					
	#vector of all data
	v<-scan(file,what="character",skip=skp) 

	if(length(which(v == "Time(hh:mm:ss)"))== 1){
		Vmax<-which(v=="Time(hh:mm:ss)")-1
	} else {
		Vmax<-which(v=="~End")-1
	}
	v<-v[1:Vmax]

#processing for endpoint reads
	if(Vmax==97){
		data.Temp <- v[1]
		#data.time <- rep(0,each=nr*nc)
		v<-v[2:Vmax]
		#convert OD to numeric
		v <- as.numeric(v)
		data.Temp <- as.numeric(data.Temp)
		#output
		data.rnames<-rep(rnames,each=nc,length.out=length(v))
		data.cnames<-rep(cnames,length.out=length(v))
		d <- data.frame("OD"=v, "row"=data.rnames, "column"=data.cnames, "Temp"=data.Temp)	
		return(d)
	}

#processing for kinetic reads
	else{

#test for early termination of run and cut data if appropriate
	if(any(v=="0.00") == TRUE){
		warning("Softmax file appears to be incomplete. Processing data up to last full read only.")
		new_end <- match("0.00", v) - 2
		v <- v[1:new_end]	
	}		

	S<-2+nr*nc
	it<-c((0:(length(v)/S))*S)+1										#index of times in v
	it<-it[-length(it)]											#delete last (border effect)
	for(i in it) if(nchar(v[i])<=5)v[i]<-sprintf("00:%s",v[i])
	time<-vector()
	time<-sapply(v[it],fun)
	Temp<-as.numeric(v[it+1])										#Temperature as a factor
	v<-as.numeric(v[c(-it,-it-1)])										#Remove time and temperature from data. OD is numeric

	data.time<-rep(time,each=nr*nc,length.out=length(v))
	data.Temp<-rep(Temp,each=nr*nc,length.out=length(v))
	data.rnames<-rep(rnames,each=nc,length.out=length(v))
	data.cnames<-rep(cnames,length.out=length(v))
	
	d<-data.frame("OD"=v,"time"=data.time,"Temp"=data.Temp,"row"=data.rnames,"column"=data.cnames)
	#write.csv(d,sprintf("RdataFrame-%s.csv",file))
	return(d)
	}
}
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
#–––––––––––––––––––––––––––––––––––#–––––––––––––––––––––––––––––––––––#–––––––––––––––––––––––––––––––––––#–––––––––––––––––––––––––––––––––––#

###Example
#pexiganan<-256
#for(i in 2: 11)pexiganan[i]<-pexiganan[i-1]/2
#pexiganan<-c(pexiganan,0)
#d<-sofToR("exple.txt",rnames=c("Ec.AA1","Ec.AA2","Ec.AF2","Ec.AF2","Ec.AF3","Ec.AG1","Ec.AG2","Ec.AA2"),cnames=pexiganan)

