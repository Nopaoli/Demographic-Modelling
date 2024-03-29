

# Scripts to estimate several one and two population statistics based on 1D and 2D sfs estimated in windows across the genome

## NS_BS_250kb_win_dxy.sh

This scripts estimates the 1D and 2D sfs in windows (in this specific case 250kb) across the genome. It has the advantage that , while slow, the script avoids creating any large file, and uses only a few MB of memory to analyse  RAD or whole genome data. It will output, from every population, a tab delimited file with the 1D SFS for each window, as well as 2D sfs for the population pair. The script generates unfolded SFSs, but since polarization of the SFS is irrelevant to the statistics calculated, the use of an outgroup is not necessary and the SFS can be polarized using the refernce sequence. Some of these stats (pi, Tajima's D) can be estimated directely in ANGSD as well (this script provides nearly identical results). *d*<sub>xy</sub> can't, and other approaches to calculate *d*<sub>xy</sub> require having a VCF with all variant and invariant sites, which are huge and generally a pian in the neck to generate and work with. Furthermore, hard-calling variants require high coverage depth, and there are issues when using data with different depth, from different sources, with different amounts of missing data.  

You will need to change most of the variables defined at the beginning to your specicifc system and population names, but once you do that it should work on anything.

## Diversity_from_wsfs.R
This script takes the output of the previous script and calcuate, for each window, pi, Whatterson's theta, the number of segregating sites and Tajima's *D*. No need to estimate sample sizes, since this infomation is already included in the 1D-SFS

## dxy_from_wsfs.R

This script calculate *d*<sub>xy</sub> based on the windowed 2D-SFS. It's based on a script orginally provided to me by Reto Burri, but has been modified to work on the input generated from the first script. You need to provide info on the number of individual in each population, in the same order used to generate the 2d-SFS

## Stats_from_SFS_TD.R

This script does something different, Rather than expecting the SFS for each window in separate lines, it expect one SFS file for each population, and outputs basic summary stats for all populations in a tab delimited text file. It assumes every SFS file for each population is present in the same folder and named "POP.sfs". Each file has a single line, exactely as outputted by ANGSD. 
