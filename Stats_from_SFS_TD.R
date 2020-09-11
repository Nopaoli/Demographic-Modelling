
### Read the sfs (in vector format, txt file as space separated as output by angsd)
### In this example, your working directory is the folder where all sfs files for 3spine are
### We expect all sfs files to be called "POPNAME.sfs"

file.names <- dir("./", pattern ="*.sfs")
out.tab        <- data.frame(matrix(ncol=7,nrow=length(file.names)))
names(out.tab) <-c( "POP","N", "nsites","S","pi", "theta","TajimaD")

for(i in 1:length(file.names)){
  # Load the SFS in vector format
  dat <- scan(file.names[i])
  # get lenght of sfs (i.e. nseq+1)
  n=length(dat)
 
  
  
  ### Now calculate the n-1 harmmonic number (needed for theta calculation) 
  ### Remember that Watterson theta = K/n-1nth harmonic number, where K is the number fo segregating sites 
  ### and  the number of chromosomes in a pop
  numChrom=n-1
  harmonicNumber = 0
  for (j in 1:(numChrom - 1)) {
    harmonicNumber = harmonicNumber + 1.0/j
  }
  # now get the terms to calculate Tajima's D
  
  a1 = sum(1/seq(1,numChrom))
  a2 = sum(1/seq(1,numChrom)^2)
  b1 = (numChrom+1)/(3*(numChrom-1))
  b2 = 2*(numChrom^2 + numChrom + 3)/(9*numChrom * (numChrom-1))
  c1 = b1 - 1/a1
  c2 = b2 - (numChrom+2)/(a1*numChrom) + a2/a1^2
  e1 = c1/a1
  e2 = c2/(a1^2+a2)
  
  # get allele frequency weights for calculating pi
  p = seq(0,numChrom)/numChrom
  # now calculate all stats
  out.tab$pi[i] <- numChrom/(numChrom-1)*2*sum(dat*p*(1-p))/sum(dat)
  out.tab$theta[i]<-sum(dat[2:(length(dat)-1)])/(harmonicNumber) /sum(dat)
  out.tab$S[i]<-sum(dat[2:(length(dat)-1)])
  out.tab$nsites[i]<-sum(dat)
  out.tab$POP[i]<-file.names[i]
  out.tab$POP[i]<-gsub(".sfs", "", out.tab$POP[i])
  out.tab$N[i]<-length(dat)-1
  # now calculate Tajima's D
  #C = sqrt((c1/a1)*out.tab$S[i] + c2/(a1^2 + a2) * out.tab$S[i]*(out.tab$S[i]-1))
  P<-numChrom/(numChrom-1)*2*sum(dat*p*(1-p))
  W<-(out.tab$S[i]/(harmonicNumber))
  
  out.tab$TajimaD[i]<-(P - W)/sqrt(e1*out.tab$S[i]+e2*out.tab$S[i]*(out.tab$S[i]-1))
  
}

write.table(out.tab,"Diversity_stats.txt", quote = F, row.names = F, sep = "\t")
