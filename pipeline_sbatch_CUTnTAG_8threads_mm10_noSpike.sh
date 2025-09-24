#!/bin/bash
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=4G


# usage
# run_CUTnTAG_mm10_noSpike.sh <full path to file>.NO.EXTENSION
# example:
# ./run_CUTnTAG_mm10_noSpike.sh /media/teamgreenberg/Stagiaire/2021-Line-CUT-TAG/201230_X591_FCHF2YLCCX2_L4_CF1-FLAG
# extension in these cases MUST BE "_R1.fastq.gz"

# load functions coded in separate files
source /home/julien/sra2bw_functions.sh

# ---------------- CONFIGURE VARIABLES ----------------- #
# You can range threads from 10 to 20 depending on the availability | And thread mem up to 4G

FLAG=1540
MIN_MAPQ=10
THREADS=8
THREAD_MEM=4G
BLACKLIST="/ribes/01_Lab_Common_Resources/02_Lab_HPscreens/mm10/mm10_blackList_ENCFF547MET.bed"



# ------------------ CALL FUNCTIONS -------------------- #

# takes as input the full path to the fastq file
setupVariables $1
trimReads 5
run_fastQC
alignBowtie2_cutNtag /ribes/01_Lab_Common_Resources/02_Lab_HPscreens/mm10/mm10
groomSam
trueStats
bamToBigWigDeeptoolsCPMsmoothKeepDup 1 0
rename_cleanup


