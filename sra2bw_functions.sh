#!/bin/bash
source /home/julien/epipax.config

# usage
# fastqToBigWig.sh <full path to file>.NO.EXTENSION
# fastqToBigWig.sh /brcwork/FASTQ/ICM_ChIP/BC_ICM_H3K4me3
# fastqToBigWig.sh /brcwork/FASTQ/ICM_ChIP/BC_ICM_RNA_rep1
# extension in these cases MUST BE "_R1.fastq.gz"

# UPDATE 28 Feb 2018: remove --local option for bowtie2 alignments
# UPDATE 31 May 2018: added illuminaclip + adapters.fa for trim
# UPDATE 06 Jul 2018: reinstate --local option for bowtie2 alignments
# UPDATE 11 Sep 2020: use lessons learned from Tiff's script to improve my functions
# UPDATE 25 Feb 2021: add function for CUTnTAG alignments as per Henikoff lab recommendations
# UPDATE 25 Feb 2021: add pre-alignment PCR duplicate read removal. Only using for CnT for now
# UPDATE 17 Mar 2021: add multithreading to TrueStats
# UPDATE 07 Jun 2021: do not remove duplicate alignments for miseq enzymatic-converted DNA sequencing (ALB&AMS)
# UPDATE 06 Sep 2021: add function removePCRduplicatesPreAlignment and function calculateErrorFilter
# UPDATE 06 Sep 2021: add trim version 5, where clumpify is run by default (it's trimV4 but with new functions)
# UPDATE 06 Dec 2021: add trim version and reference genome name used to FILE name
# UPDATE 15 Dec 2021: add condition to handle read lengths <50
# UPDATE 31 Jan 2022: change CUTnTAG alignment to automatically detect short (<25bp) reads for end-to-end alnmnt
# UPDATE 11 May 2022: remove previous update on 31 Jan. There was a mistake on my mistake

# TODO incorporate XS field in STAR alignment (include SAM intron motif argument): why?
# TODO add reference genome name to file names

DEPENDENCIES=($FASTERQDUMP "$TRIMMOMATIC" "$STAR" $BISMARK $BOWTIE2 $SAMTOOLS "$PICARD" awk $BAM2FASTQ $BEDGRAPHTOBW $BAMCOVERAGE $FASTQC $FASTQSCREEN $STRINGTIE $CLUMPIFY)


# -------------------------- DEFINE FUNCTIONS ------------------------------- #

# ------- PIPELINE AND VARIABLE SETUP, AND FILE HANDLING FUNCTIONS ---------- #

# this is where the 'global' variables used in almost all functions are set
function setupVariables {
	echo "SCRATCH: $SCRATCH"
	PATH_TO_FILE=$1
	FILE=$(basename "$PATH_TO_FILE")
	mkdir -p $SCRATCH/$FILE
	cd $SCRATCH/$FILE
	LOG_FILE="$FILE"_"$(date '+%y-%m-%d')"_log.txt

	# check whether data are paired-end or single-end
	if [[ ! -z $(ls "$PATH_TO_FILE"*_R2.fastq.gz) ]]; then #true if variable exists
		PAIRED_END=true
		printProgress "[setupVariables] Data are paired-end."
		ln -s "$PATH_TO_FILE"_R1.fastq.gz ./
		ln -s "$PATH_TO_FILE"_R2.fastq.gz ./
	else
		PAIRED_END=false
		printProgress "[masterDownload] Data are single-end."
		ln -s "$PATH_TO_FILE".fastq.gz ./

	fi

	checkDependencies
	echo "PAIRED: $PAIRED_END"
	printProgress "[setUpLogFile] Script started at [$(date)]"
}


function printProgress {
	echo $1 at [`date`] | tee -a $LOG_FILE #using tee will also show in stdout
}

# Checking dependencies of the functions
function checkDependencies {
	printProgress "[checkDependencies] Checking Dependencies [$(date)]"
	EXIT=0
	for COMMAND in "${DEPENDENCIES[@]}"; do
		printProgress "[setup checkDependencies] $COMMAND..."
		command -v $COMMAND > /dev/null 2>&1 || {
			echo -e >&2 "\t\t$COMMAND not found!"
			EXIT=1
		}
	done

	if [[ $EXIT = 1 || $DEPEND = 1 ]] ; then
		exit 1
	fi
}

# exit pipeline if input file does not exist
function checkFileExists {
    if [ ! -f "$1" ]; then
        echo "Error: File $1 does not exist" | tee -a $LOG_FILE
        exit 1
    fi
}

# checks whether at least 1 read alignment exists for input bam file
function checkBamExists () {
	if [[ ! -f $1 || $($SAMTOOLS view $1 | wc -l | cut -d ' ' -f4) == 0 ]]; then
		echo "ERROR: Bam file $1 does not exist or contains no aligned reads" | tee -a $LOG_FILE
		exit 1
	fi
}

# call reads as PCR duplicates if they have less than ErrorFilter
# assumes 1% sequencing error rate
function calculateErrorFilter () {
	if $PAIRED_END; then
		READ_LENGTH=$(zcat "$FILE"_R1.fastq.gz | head -n 2 | tail -n 1 | wc -m)
	else
		READ_LENGTH=$(zcat "$FILE".fastq.gz | head -n 2 | tail -n 1 | wc -m)
	fi

	if (( $READ_LENGTH<50 )); then SUBS_COUNT=0; else SUBS_COUNT=$(echo "scale=2;$READ_LENGTH * 0.01 + 0.5" | bc | cut -f1 -d '.'); fi
	printProgress "Reads are length: $READ_LENGTH"
	printProgress "Reads with same sequence $SUBS_COUNT or fewer substitutions will be marked as PCR duplicates"
}

# ----------------------- PRE-ALIGNMENT FUNCTIONS --------------------------- #

# identify and remove duplicates before alignment and before read trimming
function removePCRduplicatesPreAlignment () {
	if $PAIRED_END; then
		checkFileExists "$FILE"_R1.fastq.gz
		checkFileExists "$FILE"_R2.fastq.gz
	else # single-end
		checkFileExists "$FILE".fastq.gz
	fi
	# run fastqc on raw FASTQs
	for x in *fastq.gz; do ln -s $x ${x//.fastq.gz/_raw.fastq.gz}; done
	printProgress "[fastQC] started on raw FASTQs"
	$FASTQC --outdir ./ --format fastq --threads $THREADS "$FILE"*_raw.fastq.gz
	for x in *_raw.fastq.gz; do unlink $x; done
	printProgress "[fastQC] on raw FASTQs finished succesfully. Will remove PCR duplicates now"
	calculateErrorFilter
	if $PAIRED_END; then
		($CLUMPIFY in1="$FILE"_R1.fastq.gz in2="$FILE"_R2.fastq.gz out1="$FILE"_dedupe_R1.fastq.gz out2="$FILE"_dedupe_R2.fastq.gz reorder dedupe=t k=19 passes=6 subs=$SUBS_COUNT -Xmx"$THREAD_MEM") &>> "$FILE"_clumpLog.txt
		tail -n 6 "$FILE"_clumpLog.txt >> "$FILE"_trimLog.txt
		# unlink removes links and not files. good for not fucking up and removing your raw FASTQs
		unlink "$FILE"_R1.fastq.gz
		unlink "$FILE"_R2.fastq.gz
			if [ -f "$FILE"_R1.fastq.gz ]; then
				printProgress "Error: File "$FILE"_R1.fastq.gz was not backed up. Quitting cumplify stage and pipeline"
				exit 1
			else
				mv "$FILE"_dedupe_R1.fastq.gz  "$FILE"_R1.fastq.gz
				mv "$FILE"_dedupe_R2.fastq.gz  "$FILE"_R2.fastq.gz
			fi
	else # Single-end
		checkFileExists "$FILE".fastq.gz
		($CLUMPIFY in="$FILE".fastq.gz out="$FILE"_dedupe.fastq.gz reorder dedupe=t k=19 passes=6 subs=$SUBS_COUNT -Xmx"$THREAD_MEM") &>> "$FILE"_clumpLog.txt
		tail -n 6 "$FILE"_clumpLog.txt >> "$FILE"_trimLog.txt
		checkFileExists "$FILE"_dedupe.fastq.gz
		# unlink removes links and not files. good for not fucking up and removing your raw FASTQs
		unlink "$FILE".fastq.gz
			if [ -f "$FILE".fastq.gz ]; then
				printProgress "Error: File "$FILE".fastq.gz was not backed up. Quitting cumplify stage and pipeline"
				exit 1
			else
				mv "$FILE"_dedupe.fastq.gz "$FILE".fastq.gz
			fi
	fi
#	rm *clump*temp*.fastq.gz
	printProgress "Removal of PCR duplicates complete!"
}




# choose trimming version when calling the trimReads_run function
function trimReads {

	INPUT_TRIMVER=$1
	# TRIM VERSIONS:
	# 0 : do not trim
	# 1 : default trimming parameters
	# 2 : for shorter reads
	# 3 : for PBAT with random octomers

	# skip trimming
	if [[ $INPUT_TRIMVER == "0" ]]; then
		TRIMVER="trimV0"
		mv "$FILE"_R1.fastq.gz "$FILE"_"$TRIMVER"_R1.fastq.gz
		mv "$FILE"_R2.fastq.gz "$FILE"_"$TRIMVER"_R2.fastq.gz

	# for standard HTS alignment
	elif [[ $INPUT_TRIMVER == "1" ]]; then
		TRIMVER="trimV1"
		TRIM_PARAMS="ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36"
		trimReads_run

	# for SHORT read alignment
	elif [[ $INPUT_TRIMVER == "2" ]];then
		TRIMVER="trimV2"
		TRIM_PARAMS="ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL:2:30:10 SLIDINGWINDOW:4:20 MINLEN:20"
		trimReads_run

	# for PBAT random octomer
	elif [[ $INPUT_TRIMVER == "3" ]];then
		TRIMVER="trimV3"
		TRIM_PARAMS="HEADCROP:8 ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36"
		trimReads_run

	# for CUT and TAG
	elif [[ $INPUT_TRIMVER == "4" ]];then
		TRIMVER="trimV4"
		TRIM_PARAMS="ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL:2:30:10 SLIDINGWINDOW:4:20 MINLEN:20"


		# for CUTnTAG, we identify and remove duplicates before alignment and before read trimming
		if $PAIRED_END; then
			($CLUMPIFY in1="$FILE"_R1.fastq.gz in2="$FILE"_R2.fastq.gz out1="$FILE"_dedupe_R1.fastq.gz out2="$FILE"_dedupe_R2.fastq.gz reorder dedupe=t k=19 passes=6) 2>> "$FILE"_trimLog.txt
			# unlink removes links and not files. good for not fucking up and removing your raw FASTQs
			unlink "$FILE"_R1.fastq.gz
			unlink "$FILE"_R2.fastq.gz
			if [ -f "$FILE"_R1.fastq.gz ]; then
				echo "Error: File "$FILE"_R1.fastq.gz was not backed up. Quitting CUTnTAG pipeline" | tee -a $LOG_FILE
				exit 1
			else
				mv "$FILE"_dedupe_R1.fastq.gz  "$FILE"_R1.fastq.gz
				mv "$FILE"_dedupe_R2.fastq.gz  "$FILE"_R2.fastq.gz
			fi

		else
			($CLUMPIFY in="$FILE".fastq.gz out="$FILE"_dedupe.fastq.gz reorder dedupe=t k=19 passes=6) 2>> "$FILE"_trimLog.txt
			# unlink removes links and not files. good for not fucking up and removing your raw FASTQs
			unlink "$FILE".fastq.gz
			if [ -f "$FILE".fastq.gz ]; then
				echo "Error: File "$FILE".fastq.gz was not backed up. Quitting CUTnTAG pipeline" | tee -a $LOG_FILE
				exit 1
			else
				mv "$FILE"_dedupe.fastq.gz "$FILE".fastq.gz
			fi
		fi
		rm *clump*temp*.fastq.gz
		trimReads_run

	# For standard HTS alignment but with pre-alignment PCR removal
	elif [[ $INPUT_TRIMVER == "5" ]];then
		TRIMVER="trimV5"
		TRIM_PARAMS="ILLUMINACLIP:$ILLUMINA_ADAPATORS_ALL:2:30:10 SLIDINGWINDOW:4:20 MINLEN:24"
		removePCRduplicatesPreAlignment
		trimReads_run



	else
		printprogress "[trimReads] Please provide a parameter version for read trimming"
		exit 1
	fi
}


function trimReads_run {
	printProgress "[trimReads] Started using version $TRIMVER parameters"
	echo "$TRIM_PARAMS"
	if $PAIRED_END; then
		checkFileExists "$FILE"_R1.fastq.gz
		checkFileExists "$FILE"_R2.fastq.gz
		printProgress "[trimReads] Trimming "$FILE"_R1.fastq.gz and "$FILE"_R2.fastq.gz..."
		( $TRIMMOMATIC PE -threads $THREADS \
		-summary "$FILE"_trimSummary.txt \
		"$FILE"_R1.fastq.gz "$FILE"_R2.fastq.gz \
		"$FILE"_"$TRIMVER"_R1.fastq.gz "$FILE"_unpaired_trim_R1.fastq.gz "$FILE"_"$TRIMVER"_R2.fastq.gz "$FILE"_unpaired_trim_R2.fastq.gz \
		$TRIM_PARAMS ) 2>> "$FILE"_trimLog.txt

		# Keep unpaired filtered reads for WGBS data
		if [[ $FILE == *"Bisulfite-Seq"* ]] || [[ $FILE == *"RRBS"* ]] || [[ $FILE == *"BSSeq"* ]] || [[ $FILE == *"PBAT"* ]] || [[ $FILE == *"DNAme"* ]] || [[ $FILE == *"WGBS"* ]]; then
			return 0
		else
			rm "$FILE"_unpaired_trim_R1.fastq.gz "$FILE"_unpaired_trim_R2.fastq.gz
		fi

	else #Single-End
		checkFileExists "$FILE".fastq.gz
		printProgress "[trimReads] Trimming "$FILE".fastq.gz..."
		( $TRIMMOMATIC SE -threads $THREADS \
		-summary "$FILE"_trimSummary.txt \
		"$FILE".fastq.gz \
		"$FILE"_"$TRIMVER".fastq.gz \
		$TRIM_PARAMS ) 2>> "$FILE"_trimLog.txt
	fi

	# Regardless of PE or SE, do this:
	cat "$FILE"_trimLog.txt "$FILE"_trimSummary.txt >> $LOG_FILE
	printProgress "[trimReads] Finished successfully"
}


function run_fastQC {
	# list all FASTQ files in order to run on trimmed fastq as well
	# i.e. see how well the trim performed
	local FASTQ_LIST=$(ls "$FILE"*fastq.gz)
	printProgress "[fastQC] started on $FASTQ_LIST"
	$FASTQC --outdir ./ --format fastq --threads $THREADS $FASTQ_LIST
	printProgress "[fastQC] finished succesfully"
}


# check HTS library against multiple model organisms for contamination
# run with command "--bisulfite" for WGBS data
function run_fastqScreen {
	if $PAIRED_END; then # only run on R1. no PE option available at this time
		printProgress "[fastqScreen] started on "$FILE"_"$TRIMVER"_R1.fastq.gz"
		$FASTQSCREEN --threads $THREADS "$FILE"_"$TRIMVER"_R1.fastq.gz
	else
                printProgress "[fastqScreen] started on "$FILE"_"$TRIMVER".fastq.gz"
		$FASTQSCREEN --threads $THREADS "$FILE"_"$TRIMVER".fastq.gz
	fi
	printProgress "[fastqScreen] finished successfully"
}


function run_fastqScreen_convertedDNA {
	if $PAIRED_END; then # only run on R1. no PE option available at this time
		printProgress "[fastqScreen_convertedDNA] started on "$FILE"_"$TRIMVER"_R1.fastq.gz"
		$FASTQSCREEN --bisulfite --conf /home/teamgreenberg/bin/miniconda2/share/fastq-screen-0.14.0-0/fastq_screen_bis.conf --threads $THREADS "$FILE"_"$TRIMVER"_R1.fastq.gz
	else
		printProgress "[fastqScreen_convertedDNA] started on "$FILE"_"$TRIMVER".fastq.gz"
		$FASTQSCREEN --bisulfite --conf /home/teamgreenberg/bin/miniconda2/share/fastq-screen-0.14.0-0/fastq_screen_bis.conf --threads $THREADS "$FILE"_"$TRIMVER".fastq.gz
	fi
	printProgress "[fastqScreen] finished successfully"
}



# ------------------------- ALIGNMENT FUNCTIONS ----------------------------- #

function alignBowtie2 {
	# input trimmed fastq, output SAM
	PATH_TO_GENOME=$1
	GENOME_NAME=$(basename $1)
	printProgress "[bowtie2 alignment] started with reference "$PATH_TO_GENOME""

	if $PAIRED_END; then
		checkFileExists "$FILE"_"$TRIMVER"_R1.fastq.gz
		checkFileExists "$FILE"_"$TRIMVER"_R2.fastq.gz

		($BOWTIE2 -x "$PATH_TO_GENOME" \
		-p $THREADS -1 "$FILE"_"$TRIMVER"_R1.fastq.gz -2 "$FILE"_"$TRIMVER"_R2.fastq.gz \
		--local -S "$FILE"_raw.sam) 2>> "$FILE"_"$GENOME_NAME"_bt2_alignLog.txt

	else # single-end
		checkFileExists "$FILE"_"$TRIMVER".fastq.gz
		($BOWTIE2 -x "$PATH_TO_GENOME" \
		-p $THREADS -U "$FILE"_"$TRIMVER".fastq.gz \
		--local -S "$FILE"_raw.sam) 2>> "$FILE"_"$GENOME_NAME"_bt2_alignLog.txt
	fi
	printProgress "[alignment] Bowtie2 alignment finished successfully"
}

function alignBowtie2_cutNtag {
	if $PAIRED_END ; then
		printProgress "[bowtie2 alignment] Running CUT and TAG pipeline on paired-end samples"
	else
		printProgress "[bowtie2 alignment] Trying to run CUT and TAG pipeline on single-end sample, which doesn't work! Exiting"
		exit 1
	fi

	# from calculateErrorFilter function
	printProgress "[bowtie2 alignment] starting using long read (>24bp) settings"
	BOWTIE2_CUTNTAG_PARAMETERS="--local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700"

	# input trimmed fastq, output SAM
	PATH_TO_GENOME=$1
	GENOME_NAME=$(basename $1)
	printProgress "[bowtie2 alignment] for cut n tag libraries with reference "$PATH_TO_GENOME""

	checkFileExists "$FILE"_"$TRIMVER"_R1.fastq.gz
	checkFileExists "$FILE"_"$TRIMVER"_R2.fastq.gz

	($BOWTIE2 -x "$PATH_TO_GENOME" \
	-p $THREADS -1 "$FILE"_"$TRIMVER"_R1.fastq.gz -2 "$FILE"_"$TRIMVER"_R2.fastq.gz \
	$BOWTIE2_CUTNTAG_PARAMETERS \
	-S "$FILE"_raw.sam) 2>> "$FILE"_"$GENOME_NAME"_bt2_alignLog.txt

	printProgress "[alignment] Bowtie2 alignment for CUTnTAG libraries finished successfully"
}

function alignSTAR {
	local PATH_TO_GENOME=$1
	GENOME_NAME=$(basename $1)
	printProgress "[STAR alignment] started with reference "$PATH_TO_GENOME""

	if $PAIRED_END; then
		checkFileExists "$FILE"_"$TRIMVER"_R1.fastq.gz
		checkFileExists "$FILE"_"$TRIMVER"_R2.fastq.gz

		$STAR --runMode alignReads \
		--genomeDir "$PATH_TO_GENOME" \
		--runThreadN $THREADS \
		--limitBAMsortRAM $THREAD_MEM \
		--readFilesCommand zcat \
		--readFilesIn "$FILE"_"$TRIMVER"_R1.fastq.gz "$FILE"_"$TRIMVER"_R2.fastq.gz  \
		--outFilterType BySJout \
		--outSAMtype SAM

	else # single-end
		checkFileExists "$FILE"_"$TRIMVER".fastq.gz

		$STAR --runMode alignReads \
		--genomeDir "$PATH_TO_GENOME" \
		--runThreadN $THREADS \
		--limitBAMsortRAM $THREAD_MEM \
		--readFilesCommand zcat \
		--readFilesIn "$FILE"_"$TRIMVER".fastq.gz \
		--outFilterType BySJout \
		--outSAMtype SAM

	fi
	mv Log.final.out "$FILE"_"$GENOME_NAME"_STAR_alignLog.txt
	mv Aligned.out.sam "$FILE"_raw.sam
	printProgress "[STAR alignment] finished successfully"
}



# ----------------------- POST ALIGNMENT FUNCTIONS -------------------------- #

# convert raw SAM to a PCR duplicate-marked, coordinate-sorted BAM
# NOT compatible with DNAme data
function groomSam {
	checkBamExists "$FILE"_raw.sam
	printProgress "[sam grooming] started"

	$PICARD CleanSam -I "$FILE"_raw.sam -O "$FILE"_cleaned.sam --TMP_DIR ./
	rm "$FILE"_raw.sam
	$PICARD SamFormatConverter -I "$FILE"_cleaned.sam -O "$FILE"_unsorted.bam --TMP_DIR ./ 
	rm "$FILE"_cleaned.sam
	$PICARD SortSam -SO coordinate -I "$FILE"_unsorted.bam -O "$FILE"_sorted.bam --TMP_DIR ./ 
	rm "$FILE"_unsorted.bam
	$PICARD MarkDuplicates -I "$FILE"_sorted.bam -O "$FILE".bam -M "$FILE"_markDupeMetrics.txt -REMOVE_DUPLICATES false --TMP_DIR ./ 
#	cat "$FILE"_markDupeMetrics.txt >> $LOG_FILE
	rm "$FILE"_sorted.bam
	$SAMTOOLS index -@ $THREADS $FILE.bam

	printProgress "[sam grooming] finished successfully"
}


# to be run on duplicate-marked BAM files
# NOT compatible with DNAme data
function trueStats {

        TRUESTATS_INPUT=${1:-$FILE}
	checkBamExists "$TRUESTATS_INPUT".bam
	printProgress "[trueStats] Collecting alignment statistics for $TRUESTATS_INPUT.bam"

	echo "" >> "$TRUESTATS_INPUT"_flagstats.txt
	echo ""$TRUESTATS_INPUT" trueStats" >> "$TRUESTATS_INPUT"_flagstats.txt
	echo "filtering "$TRUESTATS_INPUT" for mapped reads" >> "$TRUESTATS_INPUT"_flagstats.txt
	$SAMTOOLS view -@ $THREADS -b -F 4 "$TRUESTATS_INPUT".bam > "$TRUESTATS_INPUT"_mapped.bam
	$SAMTOOLS flagstat -@ $THREADS "$TRUESTATS_INPUT"_mapped.bam >> "$TRUESTATS_INPUT"_flagstats.txt
	echo "" >> "$TRUESTATS_INPUT"_flagstats.txt
	echo "filtering "$TRUESTATS_INPUT" for MAPQ > "$MIN_MAPQ"" >> "$TRUESTATS_INPUT"_flagstats.txt
	$SAMTOOLS view -@ $THREADS -bq "$MIN_MAPQ" "$TRUESTATS_INPUT"_mapped.bam > "$TRUESTATS_INPUT"_mapped_MAPQ.bam
	$SAMTOOLS flagstat -@ $THREADS "$TRUESTATS_INPUT"_mapped_MAPQ.bam >> "$TRUESTATS_INPUT"_flagstats.txt
	rm "$TRUESTATS_INPUT"_mapped_MAPQ.bam "$TRUESTATS_INPUT"_mapped.bam
	cat "$TRUESTATS_INPUT"_flagstats.txt >> $LOG_FILE
}


function run_preseq {
	checkBamExists "$FILE".bam
	printProgress "[preseq] started on $FILE.bam"
	$PRESEQ lc_extrap -D -B "$FILE".bam -o "$FILE"_preseq.txt
	checkFileExists "$FILE"_preseq.txt
	printProgress "[preseq] finished successfully"
#	cat "$FILE"_preseq.txt >> $LOG_FILE
}


function bamToBigWig {
	local CHROM_SIZES=$1
	checkFileExists $CHROM_SIZES
	checkBamExists "$FILE".bam

	#OPTIONAL filtering step. Set FLAG and MIN_MAPQ=0 to skip
	printProgress "[bamToWig] started with parameters F="$FLAG" and q="$MIN_MAPQ" to generate bigWigs"
	$SAMTOOLS view -bh -F "$FLAG" -q "$MIN_MAPQ" "$FILE".bam > "$FILE"_F"$FLAG"_q"$MIN_MAPQ".bam
	$BEDTOOLS genomecov -ibam "$FILE"_F"$FLAG"_q"$MIN_MAPQ".bam -bg -split > "$FILE"_F"$FLAG"_q"$MIN_MAPQ"_unsorted.bedGraph
	sort -k1,1 -k2,2n "$FILE"_F"$FLAG"_q"$MIN_MAPQ"_unsorted.bedGraph > "$FILE"_F"$FLAG"_q"$MIN_MAPQ".bedGraph
	$BEDGRAPH_TO_BIGWIG "$FILE"_F"$FLAG"_q"$MIN_MAPQ".bedGraph $CHROM_SIZES "$FILE"_F"$FLAG"_q"$MIN_MAPQ".bw
	printProgress "[bamToWig] finished successfully"
}


function bamToBigWigRPM {
	local CHROM_SIZES=$1
	checkFileExists $CHROM_SIZES
	checkBamExists "$FILE".bam

	#OPTIONAL filtering step. Set FLAG and MIN_MAPQ=0 to skip
	printProgress "[bamToWigRPM] started with parameters F="$FLAG" and q="$MIN_MAPQ" to generate RPM bigWigs"
	$SAMTOOLS view -bh -F "$FLAG" -q "$MIN_MAPQ" "$FILE".bam > "$FILE"_F"$FLAG"_q"$MIN_MAPQ".bam
	local READ_COUNT=$(samtools view -c "$FILE"_F"$FLAG"_q"$MIN_MAPQ".bam)
	local SCALING_FACTOR=$(echo "scale=25;1000000/$READ_COUNT" | bc)
	printProgress "detected "$READ_COUNT" filtered reads"
	printProgress "scaling by "$SCALING_FACTOR""
	$BEDTOOLS genomecov -ibam "$FILE"_F"$FLAG"_q"$MIN_MAPQ".bam -bg -split -scale $SCALING_FACTOR > "$FILE"_F"$FLAG"_q"$MIN_MAPQ"_unsorted_RPM.bedGraph
	sort -k1,1 -k2,2n "$FILE"_F"$FLAG"_q"$MIN_MAPQ"_unsorted_PRM.bedGraph > "$FILE"_F"$FLAG"_q"$MIN_MAPQ"_RPM.bedGraph
	$BEDGRAPH_TO_BIGWIG "$FILE"_F"$FLAG"_q"$MIN_MAPQ"_RPM.bedGraph $CHROM_SIZES "$FILE"_F"$FLAG"_q"$MIN_MAPQ"_RPM.bw
	printProgress "[bamToWigRPM] finished successfully"
}


function bamToBigWigDeeptoolsCPMsmooth {
	local BIN_SIZE=$1
	local SMOOTH_LEN=$2
	printProgress "[bamToWigCPMsmooth] Started with bin size $BIN_SIZE and smoothing over $SMOOTH_LEN bp"

	($BAMCOVERAGE \
	--binSize $BIN_SIZE \
	--smoothLength $SMOOTH_LEN \
	--minMappingQuality $MIN_MAPQ \
	-p $THREADS \
	--normalizeUsing CPM \
	--outFileFormat bigwig \
	--ignoreDuplicates \
	--blackListFileName $BLACKLIST \
	--ignoreForNormalization chrX chrM chrY \
	-b "$FILE".bam \
	--outFileName "$FILE"_rmDup_q"$MIN_MAPQ"_b"$BIN_SIZE"_s"$SMOOTH_LEN"_CPM.bw) 2>> $LOG_FILE
	printProgress "[bamToWigCPMsmooth] finished successfully"
}

function bamToBigWigDeeptoolsCPMsmoothKeepDup {
	local BIN_SIZE=$1
	local SMOOTH_LEN=$2

	printProgress "[bamToWigCPMsmooth] Started with bin size: $BIN_SIZE and smoothing over $SMOOTH_LEN bp and KEEPING DUPLICATES"

	($BAMCOVERAGE \
	--binSize $BIN_SIZE \
	--smoothLength $SMOOTH_LEN \
	--minMappingQuality $MIN_MAPQ \
	-p $THREADS \
	--normalizeUsing CPM \
	--outFileFormat bigwig \
	--blackListFileName $BLACKLIST \
	--ignoreForNormalization chrX chrM chrY \
	-b "$FILE".bam \
	--outFileName "$FILE"_kpDup_q"$MIN_MAPQ"_b"$BIN_SIZE"_s"$SMOOTH_LEN"_CPM.bw) 2>> $LOG_FILE
	printProgress "[bamToWigCPMsmooth] finished successfully"
}






function bamToBigWigDeeptoolsNoNormSmooth {
	local BIN_SIZE=$1
	local SMOOTH_LEN=$2
	printProgress "[bamToWigNoNormSmooth] Started with bin size $BIN_SIZE and smoothing over $SMOOTH_LEN bp"

	($BAMCOVERAGE \
	--binSize $BIN_SIZE \
	--smoothLength $SMOOTH_LEN \
	--minMappingQuality $MIN_MAPQ \
	-p $THREADS \
	--normalizeUsing None \
	--outFileFormat bigwig \
	--ignoreDuplicates \
	--blackListFileName $BLACKLIST \
	--ignoreForNormalization chrX chrM chrY \
	-b "$FILE".bam \
	--outFileName "$FILE"_rmDup_q"$MIN_MAPQ"_b"$BIN_SIZE"_s"$SMOOTH_LEN"_noNorm.bw) 2>> $LOG_FILE
	printProgress "[bamToWigNoNormsmooth] finished successfully"
}





function bamToBigWigDeeptoolsCPMsmooth_stranded {
	local BIN_SIZE=$1
	local SMOOTH_LEN=$2
	printProgress "[bamToWigCPMsmoothStranded] Started with bin size $BIN_SIZE and smoothing over $SMOOTH_LEN bp and stranded"

	($BAMCOVERAGE \
	--filterRNAstrand forward \
	--binSize $BIN_SIZE \
	--smoothLength $SMOOTH_LEN \
	--minMappingQuality $MIN_MAPQ \
	-p $THREADS \
	--normalizeUsing CPM \
	--outFileFormat bigwig \
	--ignoreDuplicates \
	--blackListFileName $BLACKLIST \
	--ignoreForNormalization chrX chrM chrY \
	-b "$FILE".bam \
	--outFileName "$FILE"_rmDup_q"$MIN_MAPQ"_b"$BIN_SIZE"_s"$SMOOTH_LEN"_CPM_forward.bw) 2>> $LOG_FILE

	($BAMCOVERAGE \
	--filterRNAstrand reverse \
	--scaleFactor -1 \
	--binSize $BIN_SIZE \
	--smoothLength $SMOOTH_LEN \
	--minMappingQuality $MIN_MAPQ \
	-p $THREADS \
	--normalizeUsing CPM \
	--outFileFormat bigwig \
	--ignoreDuplicates \
	--blackListFileName $BLACKLIST \
	--ignoreForNormalization chrX chrM chrY \
	-b "$FILE".bam \
	--outFileName "$FILE"_rmDup_q"$MIN_MAPQ"_b"$BIN_SIZE"_s"$SMOOTH_LEN"_CPM_reverse.bw) 2>> $LOG_FILE
	printProgress "[bamToWigCPMsmoothStranded] finished successfully"
}





function rRNAmetrics_mm10 {
	checkBamExists "$FILE".bam
	printProgress "[collecting rRNA metrics] started"
	$PICARD CollectRnaSeqMetrics --REF_FLAT /home/teamgreenberg/sra2bw/reference_genomes/mm10/mm10_RefSeq_refFlat.bed \
	--RIBOSOMAL_INTERVALS /home/teamgreenberg/sra2bw/reference_genomes/mm10/mm10_rRNA.bed \
	--I "$FILE".bam --O "$FILE"_RNAmetrics.txt --STRAND_SPECIFICITY FIRST_READ_TRANSCRIPTION_STRAND
	cat "$FILE"_RNAmetrics.txt >> "$FILE"_alignLog.txt
}


function run_stringTie {
	local GTF_FILE=$1
	checkBamExists "$FILE".bam
	printProgress "[stringTie] Performing de novo assembly of transcripts"
	$STRINGTIE "$FILE".bam \
	-G $GTF_FILE \
	-o "$FILE".gtf \
	-l $FILE \
	-p $THREADS \
	-A "$FILE"_abundance.gtf
	checkFileExists "$FILE".gtf
	printProgress "[stringTie] De novo transcriptome assembly completed successfully"
}


# ----------------------  DNA METHYLATION FUNCTIONS ------------------------ #

function alignBismark {
	# input trimmed fastq, output bam
	local PATH_TO_GENOME=$1
	GENOME_NAME=$(basename $1)

	BISMARK_ARGUMENTS="$PATH_TO_GENOME -p $BISMARK_THREAD --bowtie2 --bam " 	# assumes that samtools and bowtie2 are in the path

	if [[ $FILE == *"PBAT"* ]]; then
		BISMARK_ARGUMENTS_LOCAL=$BISMARK_ARGUMENTS"--pbat "
	else
		BISMARK_ARGUMENTS_LOCAL=$BISMARK_ARGUMENTS
	fi

	printProgress "[Bismark alignment] started with reference $PATH_TO_GENOME and arguments $BISMARK_ARGUMENTS_LOCAL"

	if $PAIRED_END; then
		checkFileExists "$FILE"_"$TRIMVER"_R1.fastq.gz
		checkFileExists "$FILE"_"$TRIMVER"_R2.fastq.gz

		printProgress "[Bismark alignment] Aligning "$FILE"_"$TRIMVER"_R1.fastq.gz and "$FILE"_"$TRIMVER"_R2.fastq.gz to genome..."
		# -un argument keeps mates that do not align, important for many WGBS libraries
		$BISMARK \
		$BISMARK_ARGUMENTS_LOCAL \
		-un \
		-1 "$FILE"_"$TRIMVER"_R1.fastq.gz \
		-2 "$FILE"_"$TRIMVER"_R2.fastq.gz

		mv "$FILE"_"$TRIMVER"_R1_bismark_bt2_pe.bam "$FILE"_PE_raw.bam
		checkBamExists "$FILE"_PE_raw.bam
		mv "$FILE"_"$TRIMVER"_R1_bismark_bt2_PE_report.txt "$FILE"_"$TRIMVER"_bismark_bt2_PE_report.txt

		#Unmapped sequences will be written to test_trimV1_R1.fastq.gz_unmapped_reads_1.fq.gz and test_trimV1_R2.fastq.gz_unmapped_reads_2.fq.gz


		printProgress "[Bismark alignment] Aligning "$FILE" unaligned read 1 and trimmed read 1 to genome..."
		# if read trimming did not occur, the cat function will print a warning
		# do not fear, it just didn't use an inexistent file
		cat "$FILE"_unpaired_trim_R1.fastq.gz "$FILE"_"$TRIMVER"_R1.fastq.gz_unmapped_reads_1.fq.gz > "$FILE"_UNPAIRED_R1.fastq.gz
		$BISMARK \
		$BISMARK_ARGUMENTS_LOCAL \
		"$FILE"_UNPAIRED_R1.fastq.gz

		rm  "$FILE"_unpaired_trim_R1.fastq.gz "$FILE"_"$TRIMVER"_R1.fastq.gz_unmapped_reads_1.fq.gz   "$FILE"_UNPAIRED_R1.fastq.gz
		# outputs: $FILE"_UNPAIRED_R1_bismark_bt2.bam"

		printProgress "[Bismark alignment] Aligning "$FILE" unaligned read 2 and trimmed read 2 to genome..."
		cat "$FILE"_unpaired_trim_R2.fastq.gz "$FILE"_"$TRIMVER"_R2.fastq.gz_unmapped_reads_2.fq.gz > "$FILE"_UNPAIRED_R2.fastq.gz
		# do not run with --pbat option if it exists. Run default in order to align to all 4 reference genome strands
		$BISMARK \
		$BISMARK_ARGUMENTS \
		"$FILE"_UNPAIRED_R2.fastq.gz

		rm  "$FILE"_unpaired_trim_R2.fastq.gz "$FILE"_"$TRIMVER"_R2.fastq.gz_unmapped_reads_2.fq.gz   "$FILE"_UNPAIRED_R2.fastq.gz
		# outputs: $FILE"_UNPAIRED_R2_bismark_bt2.bam"

		printProgress "[masterAlign Bismark] Combining unpaired SE read Bams..."
		# combine SE bams
		$SAMTOOLS merge \
		--threads $THREADS \
		"$FILE"_SE_raw.bam \
		"$FILE"_UNPAIRED_R1_bismark_bt2.bam "$FILE"_UNPAIRED_R2_bismark_bt2.bam

		rm "$FILE"_UNPAIRED_R1_bismark_bt2.bam "$FILE"_UNPAIRED_R2_bismark_bt2.bam
		checkBamExists "$FILE"_SE_raw.bam

		printProgress "[masterAlign Bismark] Done alignments. Now to clean them up and extract information..."


		# remove duplicate reads (if not RRBS) and extract DNA methylation information
		cleanAndExtractBismark "$FILE"_PE_raw.bam
		cleanAndExtractBismark "$FILE"_SE_raw.bam

	else #Single-End
		checkFileExists "$FILE"_"$TRIMVER".fastq.gz

		printProgress "[Bismark alignment] started with reference $PATH_TO_GENOME and arguments $BISMARK_ARGUMENTS_LOCAL"
		$BISMARK \
		$BISMARK_ARGUMENTS_LOCAL \
		"$FILE"_"$TRIMVER".fastq.gz

		mv "$FILE"_"$TRIMVER"_bismark_bt2.bam "$FILE"_SE_raw.bam
		checkBamExists "$FILE"_SE_raw.bam

		cleanAndExtractBismark "$FILE"_SE_raw.bam
	fi

	BISMARK_ARGUMENTS_LOCAL="" # reset variable to avoid compounding it (e.g. --pbat --pbat --pbat in the case of 3 replicate samples)
	printProgress "[masterAlign Bismark] Alignment of $FILE completed!"
	echo '' | tee -a >> $LOG_FILE && echo '' | tee -a >> $LOG_FILE

}


function cleanAndExtractBismark {

	# This is called by AlignBismark. do not call in entire pipeline script!
	# EXTRACT STEP TAKES A LONG TIME. DO NOT CODE MORE STUFF INTO THIS FUNCTION
	# input "$FILE"_SE_raw.bam "$FILE"_PE_raw.bam
	# NOTE TO JRA (HIMSELF): dismark deduplicate is equivalent to picard markduplicates plus samtools filtering (with -F 1024)

	local CLEANBIS_INPUT=${1:-CLEANBIS_INPUT}
	local DEDUPED=${CLEANBIS_INPUT//.bam/.deduplicated.bam}
	local CLEANBIS_OUTPUT=${DEDUPED//.deduplicated.bam/_rmDup.bam}

	printProgress "[cleanAndExtractBismark] started with file $CLEANBIS_INPUT"

	if [[ $CLEANBIS_INPUT == *"_PE_"* ]]; then
		if [[ $CLEANBIS_INPUT == *RRBS* || $CLEANBIS_INPUT == *MiSEQ* ]]; then
			printProgress "[cleanAndExtractBismark] Data was generated by RRBS. Duplicate read removal is not recommended for these data. Skipping. Ugly but keeping Deduplicated in name for downstream compatibility"
			cp $CLEANBIS_INPUT $CLEANBIS_OUTPUT
		else
			printProgress "[cleanAndExtractBismark] Deduplicating reads from PE file $CLEANBIS_INPUT"
			$BISMARK_DEDUPLICATE -p \
			--bam $CLEANBIS_INPUT

			mv $DEDUPED $CLEANBIS_OUTPUT
			cat ${DEDUPED//.deduplicated.bam/.deduplication_report.txt} >> $LOG_FILE
		fi
		printProgress "[cleanAndExtractBismark] Extracting DNA methylation data from PE BAM $CLEANBIS_OUTPUT"
		checkBamExists $CLEANBIS_OUTPUT
		$BISMARK_METH_EXTRACT -p \
		--multicore $THREADS \
		--buffer_size 80% \
		--comprehensive --merge_non_CpG --bedGraph --counts --cytosine_report --gzip \
		--genome_folder $PATH_TO_GENOME \
		$CLEANBIS_OUTPUT

		#test_WGBS_SE_raw_rmDup.CpG_report.tx
		checkFileExists ${CLEANBIS_OUTPUT/.bam/.CpG_report.txt.gz}
		printProgress "[cleanAndExtractBismark] Done extracting DNA methylation data from PE BAM $CLEANBIS_OUTPUT"



	else #_*SE*_
		if [[ $CLEANBIS_INPUT == *RRBS* || $CLEANBIS_INPUT == *MiSEQ* ]]; then
			printProgress "[cleanAndExtractBismark] Data was generated by RRBS. Duplicate read removal is not recommended for these data. Skipping. Ugly but keeping Deduplicated in name for downstream compatibility"
			cp $CLEANBIS_INPUT $CLEANBIS_OUTPUT
		else
			printProgress "[cleanAndExtractBismark] Deduplicating reads from SE file $CLEANBIS_INPUT"
			$BISMARK_DEDUPLICATE -s \
			--bam $CLEANBIS_INPUT

			mv $DEDUPED $CLEANBIS_OUTPUT
			cat ${DEDUPED//.deduplicated.bam/.deduplication_report.txt} >> $LOG_FILE
		fi
		printProgress "[cleanAndExtractBismark] Extracting DNA methylation data from SE BAM $CLEANBIS_OUTPUT"
		checkBamExists $CLEANBIS_OUTPUT
		$BISMARK_METH_EXTRACT -s \
		--multicore $THREADS \
		--buffer_size 80% \
		--comprehensive --merge_non_CpG --bedGraph --counts --cytosine_report --gzip \
		--genome_folder $PATH_TO_GENOME \
		$CLEANBIS_OUTPUT

		checkFileExists ${CLEANBIS_OUTPUT/.bam/.CpG_report.txt.gz}
		printProgress "[cleanAndExtractBismark] Done extracting DNA methylation data from SE BAM $CLEANBIS_OUTPUT"
	fi
}


function collectBismarkStats {

	# collapse SE PE files to avoid doubling stats files
	printProgress "[collapseSEPEbismark] Combining $FILE SE and PE files"
	if $PAIRED_END; then

		cat "$FILE"_PE_raw_rmDup.CpG_report.txt.gz "$FILE"_SE_raw_rmDup.CpG_report.txt.gz > "$FILE"_SEPE_raw_rmDup.CpG_report.txt.gz
		$SAMTOOLS merge \
		--threads $THREADS \
		"$FILE"_SEPE_raw.bam \
		"$FILE"_SE_raw.bam "$FILE"_PE_raw.bam

		checkBamExists "$FILE"_SEPE_raw.bam

#		cat Non_CpG_context_"$FILE"_PE_raw_rmDup.txt Non_CpG_context_"$FILE"_SE_raw_rmDup.txt > Non_CpG_context_"$FILE"_raw_rmDup_SEPE.txt
#		cat CpG_context_"$FILE"_PE_raw_rmDup.txt CpG_context_"$FILE"_SE_raw_rmDup.txt > CpG_context_"$FILE"_raw_rmDup_SEPE.txt
		# filtering these files should greatly reduce the total size of the temporary directory
		# unless I am mistaking, this is the only time we need the Non_CpG and CpG_context files
		zcat Non_CpG_context_"$FILE"_PE_raw_rmDup.txt.gz Non_CpG_context_"$FILE"_SE_raw_rmDup.txt.gz | grep -e "J02459.1" > Non_CpG_context_"$FILE"_raw_rmDup_SEPE.txt
		zcat CpG_context_"$FILE"_PE_raw_rmDup.txt.gz CpG_context_"$FILE"_SE_raw_rmDup.txt.gz | grep -e "J02459.1" > CpG_context_"$FILE"_raw_rmDup_SEPE.txt
		checkFileExists Non_CpG_context_"$FILE"_raw_rmDup_SEPE.txt
		checkFileExists CpG_context_"$FILE"_raw_rmDup_SEPE.txt

		rm "$FILE"_PE_raw_rmDup.CpG_report.txt.gz "$FILE"_SE_raw_rmDup.CpG_report.txt.gz \
		"$FILE"_SE_raw.bam "$FILE"_PE_raw.bam "$FILE"_PE_raw_rmDup.bam "$FILE"_SE_raw_rmDup.bam \
		Non_CpG_context_"$FILE"_PE_raw_rmDup.txt.gz Non_CpG_context_"$FILE"_SE_raw_rmDup.txt.gz # add dash here
#		CpG_context_"$FILE"_PE_raw_rmDup.txt.gz CpG_context_"$FILE"_SE_raw_rmDup.txt.gz
	fi

	printProgress "[collectBismarkStats] Collecting bisulphite conversion rate and coverage for files $(ls *raw.bam)"

	printProgress "[collectBismarkStats] Counting the number of lamba methylated (Z,X,H) and unmethylated (z,x,h) Cs"
	echo "Conversion rate (%) = converted Cs / (converted+unconverted C)s" >> "$FILE"_countingCs.txt

	# we still have to filter on J0 for the SE libraries
	cat CpG_context_"$FILE"*.txt Non_CpG_context_"$FILE"*.txt \
	| grep -e "J02459.1" \
	| cut -f5 | sort | uniq -c \
	>> "$FILE"_countingCs.txt
	rm Non_CpG_context_"$FILE"*.txt
	rm CpG_context_"$FILE"*.txt

	if [[ $(wc -l "$FILE"_countingCs.txt | awk '{print $1}') == 1 ]]; then
		printProgress "[collectBismarkStats] No lambda-phase spike in DNA found. Uh oh!"
	else
		awk '{ if ($2=="Z" || $2=="H" || $2=="X") summet+=$1; else sumunmet+=$1} END {print "Conversion rate: "sumunmet/(sumunmet+summet)"\tUnconverted (methylated): "summet"\tConverted (unmethylated): "sumunmet}' \
		"$FILE"_countingCs.txt \
		>> $LOG_FILE

	fi
	cat "$FILE"_countingCs.txt >> $LOG_FILE


	if $PAIRED_END; then
		$PICARD SortSam -I "$FILE"_SEPE_raw.bam -O "$FILE"_SEPE_sort.bam -SO coordinate --VALIDATION_STRINGENCY SILENT
		$PICARD MarkDuplicates -I "$FILE"_SEPE_sort.bam -O "$FILE"_SEPE_mrkDup.bam -M trash.txt --VALIDATION_STRINGENCY SILENT
		trueStats "$FILE"_SEPE_mrkDup
		$SAMTOOLS view -bh -F 1540 "$FILE"_SEPE_mrkDup.bam > "$FILE"_SEPE_rmDup.bam

		$SAMTOOLS coverage \
		"$FILE"_SEPE_rmDup.bam \
		| grep -v "random" \
		| grep -v "Un" \
		>> $LOG_FILE

		$SAMTOOLS coverage \
		"$FILE"_SEPE_rmDup.bam \
		| grep -v "random" \
		| grep -v "Un" \
		| grep -e "chr[1-9]" \
		| awk '{sum+=$7} END { print "Autosomal avg coverage = ",sum/NR}' \
		>> $LOG_FILE

		rm "$FILE"_SEPE_sort.bam "$FILE"_SEPE_mrkDup.bam "$FILE"_SEPE_rmDup.bam trash.txt

	else # single-end
		$PICARD SortSam -I "$FILE"_SE_raw.bam -O "$FILE"_SE_sort.bam -SO coordinate --VALIDATION_STRINGENCY SILENT
		$PICARD MarkDuplicates -I "$FILE"_SE_sort.bam -O "$FILE"_SE_mrkDup.bam -M trash.txt --VALIDATION_STRINGENCY SILENT
		trueStats "$FILE"_SE_mrkDup
		$SAMTOOLS view -bh -F 1540 "$FILE"_SE_mrkDup.bam > "$FILE"_SE_rmDup.bam
		$SAMTOOLS coverage \
		"$FILE"_SE_rmDup.bam \
		| grep -v "random" \
		| grep -v "Un" \
		>> $LOG_FILE

		$SAMTOOLS coverage \
		"$FILE"_SE_rmDup.bam \
		| grep -v "random" \
		| grep -v "Un" \
		| grep -e "chr[1-9]" \
		| awk '{sum+=$7} END { print "Autosomal avg coverage = ",sum/NR}' \
		>> $LOG_FILE

		rm "$FILE"_SE_sort.bam "$FILE"_SE_mrkDup.bam "$FILE"_SE_rmDup.bam trash.txt
	fi
	printProgress "[collectBismarkStats] done!"
}


# Is this really what this does?
### merge a cytosine report into a CpG site report
function mergeTwoStrandMethylation {

	#input  "$FILE"_SEPE_rmDup.CpG_report.txt
	if $PAIRED_END; then
		INPUT_CPG_REPORT=$FILE"_SEPE_raw_rmDup.CpG_report.txt.gz"
	else
		INPUT_CPG_REPORT=$FILE"_SE_raw_rmDup.CpG_report.txt.gz"
	fi

#	INPUT_CPG_REPORT="$FILE"_*_raw_rmDup.CpG_report.txt
#	echo $INPUT_CPG_REPORT

	printProgress "[mergeTwoStrandMethylation] Started mergeTwoStrandMethylation on $INPUT_CPG_REPORT"
	checkFileExists $INPUT_CPG_REPORT

	zcat $INPUT_CPG_REPORT | sort -k1,1 -k2,2n > "$INPUT_CPG_REPORT"_sort
	awk '
	BEGIN {
			FS = "\t"
			OFS = "\t"
			CHR = ""
			FIRST_POS = 0
			METHYL = 0
			UNMETHYL = 0
			FIRST_TRI = ""
			START = 0
	}
	{
		if ($2 == FIRST_POS || $2 == FIRST_POS + 1) {
			METHYL += $4
			UNMETHYL += $5
		} else {
			if (START != 0 ) {
				if (METHYL + UNMETHYL > 0) {
					printf CHR "\t" FIRST_POS "\t" 
					printf "%6f\t", METHYL / (METHYL + UNMETHYL) * 100.0
					print METHYL, UNMETHYL, $6, FIRST_TRI
				} else {
					print CHR, FIRST_POS, ".", METHYL, UNMETHYL, $6, FIRST_TRI
				}
			}
			START = 1
			CHR = $1
			FIRST_POS = $2
			METHYL = $4
			UNMETHYL = $5
			FIRST_TRI = $7
		}
	}
	END {
		if (METHYL + UNMETHYL > 0) {
			printf CHR "\t" FIRST_POS "\t" 
			printf "%6f\t", METHYL / (METHYL + UNMETHYL) * 100.0
			print METHYL, UNMETHYL, $6, FIRST_TRI
		} else {
			print CHR, FIRST_POS, "NA", METHYL, UNMETHYL, $6, FIRST_TRI
		}
	}' "$INPUT_CPG_REPORT"_sort > ${INPUT_CPG_REPORT//.txt.gz/_mergeTwoStrands.txt}
	rm "$INPUT_CPG_REPORT"_sort

	checkFileExists ${INPUT_CPG_REPORT//.txt.gz/_mergeTwoStrands.txt}
        printProgress "[mergeTwoStrandMethylation] mergeTwoStrandMethylation on $INPUT_CPG_REPORT complete!"
}


function convertMethylationToBigWig {


	#input  "$FILE"_SEPE_rmDup.CpG_report_mergeTwoStrands.txt
	if $PAIRED_END; then
		INPUT_TWOSTRANDMERGED_CPG_REPORT=$FILE"_SEPE_raw_rmDup.CpG_report_mergeTwoStrands.txt"
	else
		INPUT_TWOSTRANDMERGED_CPG_REPORT=$FILE"_SE_raw_rmDup.CpG_report_mergeTwoStrands.txt"
	fi

	METH_TO_BIGWIG_MIN_DEPTH=$1
	CHROM_SIZES=$2
#	This doesn't work for some reason...
#	INPUT_TWOSTRANDMERGED_CPG_REPORT=${3:-INPUT_TWOSTRANDMERGED_CPG_REPORT}
	checkFileExists "$INPUT_TWOSTRANDMERGED_CPG_REPORT"

	printProgress "[CpGReportToBigWig] Started on file $INPUT_TWOSTRANDMERGED_CPG_REPORT with min coverage x$METH_TO_BIGWIG_MIN_DEPTH"
	awk -v MIN_DEPTH=$METH_TO_BIGWIG_MIN_DEPTH '
		BEGIN {
			print "#chr\tstart\tend\tstrand\tvalue"
			FS = "\t"
			OFS = "\t"
			if (MIN_DEPTH < 1)
				MIN_DEPTH = 1
		}
	$4 + $5 >= MIN_DEPTH {
			print $1, $2-1, $2, $3
	}' $INPUT_TWOSTRANDMERGED_CPG_REPORT \
	| sort -k1,1 -k2,2n \
	>  $INPUT_TWOSTRANDMERGED_CPG_REPORT"_x"$METH_TO_BIGWIG_MIN_DEPTH

	bedGraphToBigWig $INPUT_TWOSTRANDMERGED_CPG_REPORT"_x"$METH_TO_BIGWIG_MIN_DEPTH $CHROM_SIZES ${INPUT_TWOSTRANDMERGED_CPG_REPORT/.txt/_x$METH_TO_BIGWIG_MIN_DEPTH.bw}

	# count number of covered CpGs
	printProgress "[CpGReportToBigWig] Counting the number of covered CpGs by $METH_TO_BIGWIG_MIN_DEPTH"
	CPG_COUNT=$( grep -v "J02459.1" $INPUT_TWOSTRANDMERGED_CPG_REPORT"_x"$METH_TO_BIGWIG_MIN_DEPTH |  wc -l  |  awk '{print $1}' )
	echo "$CPG_COUNT CpGs covered by $METH_TO_BIGWIG_MIN_DEPTH reads in $INPUT_TWOSTRANDMERGED_CPG_REPORT" >> $LOG_FILE
	rm $INPUT_TWOSTRANDMERGED_CPG_REPORT"_x"$METH_TO_BIGWIG_MIN_DEPTH
	printProgress "[CpGReportToBigWig] $INPUT_TWOSTRANDMERGED_CPG_REPORT x$METH_TO_BIGWIG_MIN_DEPTH CpG_report -> bigWig done!"
	echo "" | tee -a $LOG_FILE

}



function cleanupBismark  () {

	printProgress "[file cleanup bismark] Moving relevant output files. Good job!"
	mkdir -p $OUTPUT_STATS_FOLDER/$FILE
	mkdir -p $OUTPUT_BIGWIG_FOLDER
	mkdir -p $OUTPUT_CPG_REPORT_FOLDER
	mkdir -p $OUTPUT_BAM_FOLDER

        for x in $FILE*; do
                mv $x ${x//$FILE/$FILE"_"$TRIMVER"_"$GENOME_NAME}
        done


	mv *CpG_report_mergeTwoStrands.txt $OUTPUT_CPG_REPORT_FOLDER/
	rm *CpG_report.txt.gz
	rm CpG_context*

	cp *.txt $OUTPUT_STATS_FOLDER/$FILE/
	cp *.html $OUTPUT_STATS_FOLDER/$FILE/
	mv *.bw $OUTPUT_BIGWIG_FOLDER/

	rm -r $SCRATCH/$FILE
}


function cleanup {

	printProgress "[file cleanup] started. Good job!"
        mkdir -p $OUTPUT_STATS_FOLDER/$FILE
        mkdir -p $OUTPUT_BIGWIG_FOLDER
        mkdir -p $OUTPUT_BAM_FOLDER
	mkdir -p $OUTPUT_DE_NOVO_ASSEMBLY_FOLDER

	cp *.txt $OUTPUT_STATS_FOLDER/$FILE/
	cp *.html $OUTPUT_STATS_FOLDER/$FILE/

	mv *.bw $OUTPUT_BIGWIG_FOLDER/
	mv "$FILE".bam* $OUTPUT_BAM_FOLDER/
	mv *.gtf $OUTPUT_DE_NOVO_ASSEMBLY_FOLDER/

	rm -r $SCRATCH/$FILE

#	OLD:
#	="/home/teamgreenberg/sra2bw/bams"
#	="/home/teamgreenberg/sra2bw/bigWigs"
#	="/home/teamgreenberg/sra2bw/stats/$FILE"
}

function rename_cleanup {

	printProgress "[file cleanup] started. Good job!"
	mkdir -p $OUTPUT_STATS_FOLDER/$FILE
	mkdir -p $OUTPUT_BIGWIG_FOLDER
	mkdir -p $OUTPUT_BAM_FOLDER
	mkdir -p $OUTPUT_DE_NOVO_ASSEMBLY_FOLDER

        cp *.html $OUTPUT_STATS_FOLDER/$FILE/
        cp *.txt $OUTPUT_STATS_FOLDER/$FILE/

	for x in $FILE*; do
		mv $x ${x//$FILE/$FILE"_"$TRIMVER"_"$GENOME_NAME}
	done

	mv *.bw $OUTPUT_BIGWIG_FOLDER/
	mv "$FILE"_"$TRIMVER"_"$GENOME_NAME".bam* $OUTPUT_BAM_FOLDER/
	mv *.gtf $OUTPUT_DE_NOVO_ASSEMBLY_FOLDER/
	rm -r $SCRATCH/$FILE

#       OLD:
#       ="/home/teamgreenberg/sra2bw/bams"
#       ="/home/teamgreenberg/sra2bw/bigWigs"
#       ="/home/teamgreenberg/sra2bw/stats/$FILE"
}


