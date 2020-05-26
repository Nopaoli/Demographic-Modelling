#!/bin/bash -l
# author: momiglia
#SBATCH --constraint="snb|hsw"
#SBATCH -o MS_fit.ou
#SBATCH -e MS_fit.err
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -t 72:00:00
#SBATCH --mem-per-cpu=9000



module load bioconda
source activate moments
export PYTHONPATH=$PYTHONPATH:$HOME/appl_taito/moments

cd /wrk/momiglia/DONOTREMOVE/MS_fit/2b_SIM/64_SYM/
while read a b c d e f g h k l; do $b $c $d $e $f $g $h $k $l > "$a".sim; done < IM_MS_TABLE.txt
