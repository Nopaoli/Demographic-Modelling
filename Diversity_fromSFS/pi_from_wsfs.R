################################################################################################################
##
## Paolo Momigliano August 2020
##
## Script to estimate pi from windowed 1D SFS
##  The required input file consists of 
## - window coordinates in the first rows in format Scaf start end (tab separated)
## - the 1D-SFS in the following rows (one row per window)
## for example: 
## CHR1 0 10000 
##
## Command line input required to run the script is:
## - input file name 
##
## Output is a table with five columns:
## - window coordinates in format scaf  start   end
## - pi 
## - number of sites in the window (sum of the SFS)
##
## To run the script:
##
################################################################################################################


# get command line arguments.
comarg <- commandArgs()
input  <- comarg[6]



# read count data per class (do not have to create matrix, as we simply multiply the vectors)
dat <- read.table(input, h=F)

n=length(dat[1,4:ncol(dat)])

# Calculate n-1th harmonic number to use for Whatterson's theta estimation
numChrom=n-1
harmonicNumber = 0
for (j in 1:(numChrom - 1)) {
  harmonicNumber = harmonicNumber + 1.0/j
}

# creates weights for each class. This assumes the spectrum is unfolded. 
# The "weights" object stores allele frequencies for each spectrum category

p = seq(0,numChrom)/numChrom

# now get the terms to calculate Tajima's D

a1 = sum(1/seq(1,numChrom))
a2 = sum(1/seq(1,numChrom)^2)
b1 = (numChrom+1)/(3*(numChrom-1))
b2 = 2*(numChrom^2 + numChrom + 3)/(9*numChrom * (numChrom-1))
c1 = b1 - 1/a1
c2 = b2 - (numChrom+2)/(a1*numChrom) + a2/a1^2
e1 = c1/a1
e2 = c2/(a1^2+a2)

## iterate through sfs per window to compute window-wise dxy

out.tab        <- data.frame(matrix(ncol=8,nrow=nrow(dat)))
names(out.tab) <-c( "scaffold","win_start","win_end","S","pi","ThetaW","nsites", "TajimaD")

##
for (i in 1:nrow(dat)){
  
  out.tab$scaffold[i] <- as.character(dat[i,1])
  out.tab$win_start[i] <- dat[i,2]
  out.tab$win_end[i] <- dat[i,3]
  d <- dat[i,4:ncol(dat)]
  out.tab$nsites[i] <- sum(d)
  
  # Number of variable sites is the sum of the SFS excluding corners
  out.tab$S[i]<-sum(d[2:(ncol(d)-1)])
  
  # Calculate Theta as K/an (where an is the harmonic numebr of n-1)

  #out.tab$ThetaW[i]<-((out.tab$S[i]/(harmonicNumber))/out.tab$nsites[i])
  W<-(out.tab$S[i]/(harmonicNumber))
  out.tab$ThetaW[i]<- W/out.tab$nsites[i]
  # calculate pi
  #out.tab$pi[i] <- numChrom/(numChrom-1)*2*sum(d*p*(1-p))/
  P<-numChrom/(numChrom-1)*2*sum(d*p*(1-p))
  out.tab$pi[i] <- P/out.tab$nsites[i]
  # now calculate Tajima's D
  #C = sqrt((c1/a1)*out.tab$S[i] + c2/(a1^2 + a2) * out.tab$S[i]*(out.tab$S[i]-1))
  
  out.tab$TajimaD[i]<-(P - W)/sqrt(e1*out.tab$S[i]+e2*out.tab$S[i]*(out.tab$S[i]-1))
  
}


## write output

out <- paste(gsub(".sfs","",basename(input)),"_pi.txt", sep="" )
write.table(out.tab, out, col.names=T, row.names=F, quote=F, sep="\t")

