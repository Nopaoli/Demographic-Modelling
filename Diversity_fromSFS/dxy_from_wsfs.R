################################################################################################################
##
## Paolo Momigliano, August 2020 (Based on a script originally provided by Reto Burri, though the expected input is different)
##
## This script calculate dxy based on the 2D-SFS calculated by ANGSD in windows across the genome
## It assumes an input structure like this: 
## ecah row starts with the scaffold ID, start and end of window, followed by the 2D-SFS
##
## Command line input required to run the script is:
## - input file name in format group1-group2.windowSize.*
## - number of individuals in group1
## - number of indivifuals in group2
##
## Output is a table with three columns:
## - window coordinates in format scaf:start-end
## - dxy
## - number of sites in the window
##
################################################################################################################


# get command line arguments.
comarg <- commandArgs()
input  <- comarg[6]
n      <- as.numeric(comarg[7])
m      <- as.numeric(comarg[8])
#outdir <- comarg[9]

# set number of states in each group
n <- n * 2 + 1
m <- m * 2 + 1


# read count data per class (do not have to create matrix, as we simply multiply the vectors)
dat <- read.table(input, h=F)


## create weights, i.e. (p1*q2)+(p2*q1) values for each count class

# allele frequency vectors for each group
P1 <- seq(0, 1, by=1/(n-1))
P2<- seq (0, 1, by=1/(m-1))
# create all pairwise combinations
	# For comparisons in which not the same number of individuals are included the sequence in which these values are generated is crucial to match the SFS/counts contained in d.
	# The angsd manual states: For 2dsfs the results is a single line, assuming we have n categories in population1 and m categories in population2, then the first m values will be the SFS for the first category in population1, etc.
weights <- vector(length=n*m)
c <- 1
for (i in 1:n){
  p1=P1[i]
  for (j in 1:m){
    p2=P2[j]
    weights[c] <- p1*(1-p2)+p2*(1-p1)
    c <- c+1
  }
}


## iterate through sfs per window to compute window-wise dxy

out.tab        <- data.frame(matrix(ncol=5,nrow=nrow(dat)))
names(out.tab) <-c( "scaffold","win_start","win_end","dxy", "nsites")
  

##
for (i in 1:nrow(dat)){
  
  out.tab$scaffold[i] <- as.character(dat[i,1])
  out.tab$win_start[i] <- dat[i,2]
  out.tab$win_end[i] <- dat[i,3]
  d <- dat[i,4:ncol(dat)]
  out.tab$nsites[i] <- sum(d)
  out.tab$dxy[i] <- sum(d * weights)/out.tab$nsites[i]
  
}
## write output

## write output

out <- paste(gsub(".sfs","",basename(input)),"_dxy.txt", sep="" )
write.table(out.tab, out, col.names=T, row.names=F, quote=F, sep="\t")

