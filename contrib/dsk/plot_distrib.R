#!/usr/bin/Rscript

#usage :  ./plot_distrib.sh  file.histo
args <- commandArgs(TRUE)
if (length(args)==0) {
print("Usage:")
print("./plot_distrib.R  file.histo")
quit()
}
name <-args[1]
pdf(paste(name,"_large.pdf",sep=""),width=7,height=7)



a<-read.table(name)
x<-a[,1]
y<-a[,2]
plot(x,y,type='l',log="y",xlim=c(0,5000),main = "Coverage distribution",xlab="coverage",ylab="Nb kmers",)


pdf(paste(name,".pdf",sep=""),width=7,height=7)
plot(x,y,type='l',log="y",xlim=c(0,200),main = "Coverage distribution",xlab="coverage",ylab="Nb kmers",)
