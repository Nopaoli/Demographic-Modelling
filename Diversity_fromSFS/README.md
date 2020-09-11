# Scripts to estimate several one and two population statistics based on 1D and 2D sfs estimated in windows across the genome

## NS_BS_250kb_win_dxy.sh

This scripts estimates the 1D and 2D sfs in windows (in this specific caes 250kb) across the genome. it has the advantange that , while slow, the script avoid creating any large file, and uses only a few MB of memory to analyse data from RAD or whole genome data. it will output, from every population, a tab delimited file with the 1D SFS for each window, as well as 2D sfs for the population pair. The script generates unfolded SFS, but since polaration of the SFS is irrelevant to the statics calculate, the use of an outgroup is not necessary and the SFS can be polarized using the refernce sequence. 

## Diversity_from_wsfs.R
this script takes the output of the previous script and calcuate, for each window, pi, Whatterson's theta, the number of segregating sites and Tajima's D

## dxy_from_wsfs.R

This script calculate dxy based on the windowed 2D-SFS. It's based on a script provided by Reto Burri


