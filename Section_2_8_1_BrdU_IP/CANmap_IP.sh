#!/bin/bash -l
scriptVersion="1.1"

#~ shell script to map Illuminia paired-end or single reads to reference, sort resulting bam file, mark duplicates 
#~ (for paired-end), index the bam file, extract number of reads mapped in bins (5', first in pair for paired end)

#~ Authors: Conrad Nieduszynski, Dzmitry Batrakou

#~ example command for a single read file:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~   CANmap.sh -g sacCer3 -U singleReadFile.fastq.gz -s T9475_G2 (-t -l -w 1000 -c 16)  ~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~ example command for paired-end read files:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~  CANmap.sh -g sacCer3 -1 firstMateFile.fastq.gz -2 secondMateFile.fastq.gz -s T9475_G2 (-t -l -w 1000 -c 16)  ~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ mandatory parameters:                                                            ~~~~~~~
#~ -g (--genome)        reference genome, indexed with bowtie                       ~~~~~~~
#~ ********************   EITHER   (paired end reads)   *************************** ~~~~~~~
#~ -1 (--first)         file(s) with first mate pair sequences in fastq format   ** ~~~~~~~
#~ -2 (--second)        file(s) with second mate pair sequences in fastq format  ** ~~~~~~~
#~ *********************   OR   (single end reads)    ***************************** ~~~~~~~
#~ -U (--reads)         file(s) with single end sequence reads in fastq format   ** ~~~~~~~
#~ ******************************************************************************** ~~~~~~~
#~ -s (--sample)        base name for the resulting files                           ~~~~~~~
#~ ________________________________________________________________________________ ~~~~~~~
#~ optional parameters:                                                             ~~~~~~~
#~ -w (--window)        the window size in basepairs, defaults to 1000 bp           ~~~~~~~
#~ -c (--cores)         number of cores to use, defaults to 16                      ~~~~~~~
#~ -t (--test)          testing mode: keep temp files, don't move raw read files &  ~~~~~~~
#~                                      also won't execute the resulting bash file  ~~~~~~~
#~ -l (--local)         run locally - you need to have all programs installed       ~~~~~~~
#~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~ surround multiple read files with quotes and separate by commas (no spaces!):
#~			"firstFile.fastq.gz,secondFile.fastq.gz"

#~ if running locally, you need to have all the necessary tools installed and available for the script. There are several options:
#~   1) all tools are installed into your search path (you still need to worry about $BOWTIE2_INDEXES);
#~   2) installation folders are added manually to your search path (again, mind $BOWTIE2_INDEXES);
#~   3) variables $BOWTIE, $BOWTIE2_INDEXES, $SAMTOOLS, $BEDTOOLS, $PICARDTOOLS are specified in ~/.profile
#~   4) the variables are specified in this script (see below).





##  Specify and uncomment as needed
# BOWTIE="/path/to/bowtie2"
# BOWTIE2_INDEXES="/path/to/bowtie2/index"
# SAMTOOLS="/path/to/samtools"
# BEDTOOLS="/path/to/bedtools"
# PICARDTOOLS="/path/to/picard.jar" ## only required for paired-end reads

# set default options
WORKING_DIRECTORY=`pwd`
if [ ! -w $WORKING_DIRECTORY ] ; then WORKING_DIRECTORY="/home/"${USER} ; echo "Not a writable directory. The files will be written to /home/${USER} instead." ; fi
theNUMBERofCORESdefault=16
theWINDOWsize=1000

# collect command line options
while [ "$1" != "" ]; do
	case $1 in
		-g | --genome )		shift				# name of reference genome (must be present in bowtie index directory
							theGENOME=$1
									;;
		-1 | --first )		shift				# file(s) with first mate pair sequences in fastq format
							firstMATEfiles=$1
									;;
		-2 | --second)		shift				# file(s) with second mate pair sequences in fastq format
							secondMATEfiles=$1	
									;;
		-U | --reads)		shift				# file(s) with single end sequence reads in fastq format
							readFILES=$1
									;;
		-s | --sample)		shift				#name for resulting files
							sampleName=$1
									;;
		-w | --window )		shift				#the window size in basepairs, defaults to 1000 bp
							theWINDOWsize=$1
									;;
		-c | --cores )		shift				# maximum number of cores available, in theory more should be faster
							theNUMBERofCORES=$1
									;;
		-t | --test )		testMode=TRUE			# testing mode: keep temp files, don't move raw read files, don't execute the bash script
									;;
		-l | --local )		runLocally=TRUE			# run without queue controls and module load
									;;
        esac
        shift
done


if [ "$runLocally" != TRUE ] ; then
	## running on the server
	# set environment variable for location of bowtie indexes on the server
	export BOWTIE2_INDEXES="/data/nieduszynski/RESOURCES/REFERENCE_GENOMES"
	# load required modules and set variables
	module load SLURM/5.08.6
	module load BOWTIE/2.2.5 && BOWTIE="bowtie2"
	module load SAMTOOLS/1.3.1 && SAMTOOLS="samtools"
	module load PICARDTOOLS/6.15 && PICARDTOOLS="/software/BIOINFORMATICS/PICARDTOOLS/picard-06.15/dist/picard.jar"
    module load BEDTOOLS/2.26.0 && BEDTOOLS="bedtools"
    module load BAMTOOLS && BAMTOOLS="bamtools"
	echo "Required modules are loaded"
else 
	## running locally, so need to check that all programs and paths are accessible
	if ! hash bowtie2 ; then
		if [ -z "$BOWTIE" ] ; then
			echo "Could not find bowtie. You must define BOWTIE variable (e.g. in your ~/.profile file)"
			exit 3
		fi
	else
		BOWTIE="bowtie2"
	fi
	if [ -z "$BOWTIE2_INDEXES" ] ; then
		echo "Path to bowtie genome indices is not set. You must define BOWTIE2_INDEXES variable (e.g. in your ~/.profile file)"
		exit 3
	fi
	if ! hash samtools ; then
		if [ -z "$SAMTOOLS" ] ; then
			echo "Could not find samtools. You must define SAMTOOLS variable (e.g. in your ~/.profile file)"
			exit 3
		fi
	else
		SAMTOOLS="samtools"
	fi
	if ! hash bedtools ; then
		if [ -z "$BEDTOOLS" ] ; then
			echo "Could not find bedtools. You must define BEDTOOLS variable (e.g. in your ~/.profile file)"
			exit 3
		fi
	else
		BEDTOOLS="bedtools"
	fi
	theNUMBERofCORESdefault=4
fi

if [ -z "$secondMATEfiles" ]; then
	if [ -z "$readFILES" ]; then
		echo "Missing read file(s)"
		exit 4
	else
		echo "Single end sequence file(s) detected: $readFILES"
	fi
else
	if [ -z "$PICARDTOOLS" ] ; then
		echo "Could not find picardtools. You must define PICARDTOOLS variable (e.g. in your ~/.profile file)"
		exit 3
	fi
	matePAIR=yes
	echo "Mate pair detected: $firstMATEfiles and $secondMATEfiles"
fi

# create subdirectories
mkdir $WORKING_DIRECTORY/${sampleName}
mkdir $WORKING_DIRECTORY/${sampleName}/processed
mkdir $WORKING_DIRECTORY/${sampleName}/raw
mkdir $WORKING_DIRECTORY/${sampleName}/tmp
# set up variables
sampleDir=$WORKING_DIRECTORY/${sampleName}
processed=$WORKING_DIRECTORY/${sampleName}/processed
raw=$WORKING_DIRECTORY/${sampleName}/raw
tmp=$WORKING_DIRECTORY/${sampleName}/tmp
theBASHfile=${sampleName}_CANmap.bash
echo "The number of cores to use is ${theNUMBERofCORES:=$theNUMBERofCORESdefault}"
if [ "$testMode" = TRUE ] ; then
	echo "Running in testing mode"
fi
if [ "$runLocally" = TRUE ] ; then
	echo "Running locally"
fi

# define the comands
if [ -n "$matePAIR" ] ; then
    command1="$BOWTIE -N 0 -p $theNUMBERofCORES -q --phred33 -x $theGENOME -U \"$readFILES\" -S ${sampleName}.sam"
#	command2="mv {$firstMATEfiles,$secondMATEfiles} ${raw} && chmod 0444 ${raw}/*"
	command3="$SAMTOOLS sort -l 9 -m 750M -O bam -o ${sampleName}.bam -T {sampleName}.sort.tmp -@ $theNUMBERofCORES ${sampleName}.sam"
#	command4="rm ${tmp}/${sampleName}.sam"
	command5="java -XX:ParallelGCThreads=$((theNUMBERofCORES-1)) -jar $PICARDTOOLS MarkDuplicates I=${sampleName}.bam O=${sampleName}.bam M=${sampleName}_picard_metrics.txt TMP_DIR=${tmp}"
	else
	command1="$BOWTIE -N 0 -p $theNUMBERofCORES -q --phred33 -x $theGENOME -U \"$readFILES\" -S ${sampleName}.sam"
#	command2="mv {$readFILES} ${raw} && chmod 0444 ${raw}/*"
#    command3="grep XM:i:0 ${tmp}/${sampleName}.sam > ${tmp}/${sampleName}.noMismatch.sam"

	command4="$SAMTOOLS sort -l 9 -m 750M -O bam -o ${sampleName}.bam -T ${sampleName}.sort.tmp -@ $theNUMBERofCORES ${sampleName}.sam"
#	command4="rm ${sampleName}.sam"
	command5="bamtools filter -tag XM:0 -in ${sampleName}.bam -out ${sampleName}_noMismatch.bam"
fi


command6="$SAMTOOLS index ${sampleName}_noMismatch.bam"

# check if the genome divided into windows already exists, and if not - make one:
#command7="test ! -e ${BOWTIE2_INDEXES}/genomeWindows/${theGENOME}.${theWINDOWsize}bp.bed && $SAMTOOLS idxstats ${sampleName}.bam | awk 'BEGIN {OFS=\"\\t\"} {if (\$2>0) print (\$1,\$2)}' > ${theGENOME}.txt  && $BEDTOOLS makewindows -g ${theGENOME}.txt -w $theWINDOWsize > ${BOWTIE2_INDEXES}/genomeWindows/${theGENOME}.${theWINDOWsize}bp.bed"

# check if the genome divided into windows already exists, and if not - make one:


# determine the number of 5’ ends of reads from the sample that map to each position in the genome (you should check the -F and -f flags in the samtools comment and the grep -v XS:i: which we use to remove reads mapping to multiple genomic locations (how these are recorded depends on the mapping software you used)
if [ -n "$matePAIR" ]; then
	command9="$SAMTOOLS view -h -@ $theNUMBERofCORES -q 30 -F 3840 -f 64 -L ${BOWTIE2_INDEXES}/genomeWindows/${theGENOME}.${theWINDOWsize}bp.bed ${sampleName}_noMismatch.bam | grep -v XS:i: | $SAMTOOLS view -@ $theNUMBERofCORES -b - | $BEDTOOLS genomecov -5 -d -ibam stdin | awk 'BEGIN {OFS=\"\\t\"} {if (\$3>0) print \$1,\$2,\$2,\"name\",\$3}' > ${sampleName}.bed"
	else
	command9="$SAMTOOLS view -h -@ $theNUMBERofCORES -q 30 -F 3840 -L ${BOWTIE2_INDEXES}/genomeWindows/${theGENOME}.${theWINDOWsize}bp.bed ${sampleName}_noMismatch.bam | grep -v XS:i:| $SAMTOOLS view -@ $theNUMBERofCORES -b - | $BEDTOOLS genomecov -5 -d -ibam stdin | awk 'BEGIN {OFS=\"\\t\"} {if (\$3>0) print \$1,\$2,\$2,\"name\",\$3}' > ${tmp}/${sampleName}.tmp.bed"
fi
# this sums reads for the control sample within the specified windows (e.g. 1000 bp) and uses awk to convert to bed format
    command10="$BEDTOOLS map -a ${BOWTIE2_INDEXES}/genomeWindows/${theGENOME}.${theWINDOWsize}bp.bed -b ${tmp}/${sampleName}.tmp.bed -null 0 -o sum | awk 'BEGIN {OFS=\"\\t\"} {if (\$4>0) print \$1,\$2,\$3,\"name\",\$4}' > ${sampleName}.bed"
# clean up
#command10="rm -r ${tmp}"

# write the bash file
	echo "#!/bin/bash" > $theBASHfile
	echo "#CANmap version $scriptVersion"
if [ $"$runLocally" != TRUE ] ; then
echo "#SBATCH --job-name=CANmap
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$theNUMBERofCORES
#SBATCH --mem-per-cpu=1
#SBATCH --time=24:00:00
#SBATCH --output=${sampleDir}/${sampleName}_outfile.%j
#SBATCH --error=${sampleDir}/${sampleName}_errfile.%j" >> $theBASHfile
fi

echo $command1
echo $command1 >> $theBASHfile
if [ "$testMode" != TRUE ] ; then
	echo $command2
	echo $command2 >> $theBASHfile
fi
echo $command3
echo $command3 >> $theBASHfile
if [ "$testMode" != TRUE ] ; then
	echo $command4
	echo $command4 >> $theBASHfile
fi
echo $command5
echo $command5 >> $theBASHfile
echo $command6
echo $command6 >> $theBASHfile
echo $command7
echo $command7 >> $theBASHfile
echo $command8
echo $command8 >> $theBASHfile
echo $command9
echo $command9 >> $theBASHfile

if [ "$testMode" != TRUE ] ; then
	echo $command10
	echo $command10 >> $theBASHfile
fi
# execute the bash file
if [ "$testMode" != TRUE ] ; then
	if [ "$runLocally" != TRUE ] ; then
		sbatch $theBASHfile
	else
		/bin/bash $theBASHfile
	fi
fi
