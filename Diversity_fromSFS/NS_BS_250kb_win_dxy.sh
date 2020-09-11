#!/bin/bash -l
# author: momiglia
#SBATCH --account=project_2000465
#SBATCH -o SM_winSFS.ou
#SBATCH -e SM_winSFS.err
#SBATCH -p small
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH -t 3-00:00:00
#SBATCH --mem-per-cpu=2000



module load bioconda
source activate angsd

########################################################################################################
#### This will read a list of 2 kb windows as ANGSD regions. This is made, for example, like this: 	####
####	bedtools makewindows -w 250000 -g GCA_003186165.1_ASM318616v1_genomic.fna.fai > wind_250kb		####
####	awk '{print  $1,":",$2,"-",$3 }' OFS="" wind_250kb > wind_250kb.rf								####
####	Then we generate temporary SAF files, and add the sfs for each windows (both 1D for each	####
####	pop and 2d. From these, we can estimates both dxy and pi with no bias using the other 		####
####	scripts in the folder.																		####
####	All filters and commands are given as variables to make the scripts easier to re-use.		####
####	The script takes time but requires virtually no memory. All summary stats can then be 		####
####	calculated from the SFSs, so you need only to run this one and you can get pi, Dxy, D etc.	####
########################################################################################################




#### Bamlist prefix , assume your Bamlists are called "Bamlis_POP"
BAMLIST="/scratch/project_2000465/Turbot_2b/All_Trimmed/ANGSD/DXY_Paolo_scripts/bamlists/Bamlist_"
####	path to Genome
GENOME="/scratch/project_2000465/Turbot_2b/Assembly/GCA_003186165.1_ASM318616v1_genomic.fna"
####	path to Folder or SAF
SAF="/scratch/project_2000465/Turbot_2b/All_Trimmed/ANGSD/DXY_Paolo_scripts/SAF"
####	Path to folder of SFS
SFS="/scratch/project_2000465/Turbot_2b/All_Trimmed/ANGSD/DXY_Paolo_scripts/SFS"
####	Path to file containing windows
WINDOWS="/scratch/project_2000465/Turbot_2b/Assembly/win250.rf"
####	Filters for calculate SAF
FILTERS=" -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -C 50"
####	Angsd commands
TODO="-dosaf 1 -GL 2"
####	Path to output Dxy

### Number of cores used if given by the "-P" flag, in this case 12
### First remove the 250kb_windows SFSs if there are, since the script will append results

rm -f ${SFS}/NS-BS_250kwin.sfs
rm -f ${SFS}/BS_250kwin.sfs
rm -f ${SFS}/NS_250kwin.sfs



while read a;do
 #### use a the widow ID number generator to label windowed temporary SAF files to be deleted at end of loop
 echo $a | sed "s,:, ,g" | sed "s,-, ,g" > $SFS/WIN_$a
 ####	Generate SAF for the window for both populations
 angsd -b ${BAMLIST}NS -out ${SAF}/NS_win_$a -doSaf 1 -GL 2 -r $a $FILTERS -anc $GENOME -ref $GENOME -P 12
 angsd -b ${BAMLIST}BS -out ${SAF}/BS_win_$a -doSaf 1 -GL 2 -r $a $FILTERS -anc $GENOME -ref $GENOME -P 12
 #### Now gerenate the 1d and 2d SFS and apped results (will have the same order as your window bed file)
 realSFS ${SAF}/NS_win_$a.saf.idx ${SAF}/BS_win_$a.saf.idx -P 12 > ${SFS}/NS-BS_${a}_250kwin.sfs
 paste -d ' ' $SFS/WIN_$a ${SFS}/NS-BS_${a}_250kwin.sfs >>  ${SFS}/NS-BS_250kwin.sfs
 realSFS ${SAF}/NS_win_$a.saf.idx -P 12 > ${SFS}/NS_${a}_250kwin.sfs
 paste -d ' ' $SFS/WIN_$a ${SFS}/NS_${a}_250kwin.sfs >>  ${SFS}/NS_250kwin.sfs
 realSFS ${SAF}/BS_win_$a.saf.idx -P 12 > ${SFS}/BS_${a}_250kwin.sfs
 paste  -d ' ' $SFS/WIN_$a ${SFS}/BS_${a}_250kwin.sfs >>  ${SFS}/BS_250kwin.sfs
 ### now remove temporary SAF files
 rm -f ${SAF}/*$a*
 rm -f ${SFS}/*$a* 
done < $WINDOWS

#### 	Now some (hopeully very few)  windows will have no results, because of missing data (ANGSD being able to calculate the SFS)
####	These will have the colums with the window info, but no SFS. We want to generate a window SFSs file with no missing data
####	The first entry of the SFS is on the fourth column (the first three are CHR, window start and window end
####	We use awk to get only lines with no missing data on the fourth column
awk -F" " '{if ($4) print $0;}' ${SFS}/BS_250kwin.sfs > ${SFS}/BS_250kwin_noMissing.sfs
awk -F" " '{if ($4) print $0;}' ${SFS}/NS_250kwin.sfs > ${SFS}/NS_250kwin_noMissing.sfs
awk -F" " '{if ($4) print $0;}' ${SFS}/NS-BS_250kwin.sfs > ${SFS}/NS-BS_250kwin_noMissing.sfs

####	Also make a list of the missing data
awk -F" " '{if ($4=="") print $0;}' ${SFS}/BS_250kwin.sfs > ${SFS}/BS_250kwin_missing.sfs
awk -F" " '{if ($4=="") print $0;}' ${SFS}/NS_250kwin.sfs > ${SFS}/NS_250kwin_missing.sfs
awk -F" " '{if ($4=="") print $0;}' ${SFS}/NS-BS_250kwin.sfs > ${SFS}/NS-BS_250kwin_missing.sfs
