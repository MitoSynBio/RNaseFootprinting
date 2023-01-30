#!/bin/bash

USAGEMSG="usage: $(basename $0) [-p protein name] [-i input data directory] [-w working directory]

Calculate C-score for WT and KO samples.

-p protein name			The name of the protein of interest.
-d data directory		The input data directory that contains the 5' end normalised wig files 
				(required input name format: PROTEIN_COND_TREAT.STRAND.5p.wig where:
				PROTEIN = name of protein of interest
				COND = WT/KO
				TREAT = UNTX/A/TI/IF
				STRAND = fwd/rev
-w working directory		The working directory where all outputs will go

The output will be stored in a directories named cscore & cutoff in the working directory.
"

# This command prints usage and the name of the script
[ $# -lt 5 ] && echo "$USAGEMSG" >&2 && exit 1

while getopts "p:d:w:" flag
do
	case ${flag} in
	(p) PROTEIN=${OPTARG};;
	(d) DATA=${OPTARG};;
	(w) WD=${OPTARG};;
	(*) echo "$0: error - unrecognized option $1" >&2; exit 1;;
	esac
done

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

SCRIPTS=$SCRIPT_DIR/scripts
CSCORE=$WD/cscore
CUTOFF=$WD/cutoff
mkdir $CSCORE

############################
# Calculate C-scores
############################
for COND in WT KO 
do
	for STRAND in fwd rev
	do
		for TREAT in UNTX A TI IF
		do
			INPUT=${DATA}/${PROTEIN}_${COND}_${TREAT}.${STRAND}.5p.wig
			OUTPUT=${CSCORE}/${PROTEIN}_${COND}_${TREAT}.${STRAND}.5p.wig.tmp

			# Remove wig header
			awk 'NR>2 {print $0}' $INPUT > $OUTPUT
		done

		UNTX=${CSCORE}/${PROTEIN}_${COND}_UNTX.${STRAND}.5p.wig.tmp
		A=${CSCORE}/${PROTEIN}_${COND}_A.${STRAND}.5p.wig.tmp
		TI=${CSCORE}/${PROTEIN}_${COND}_TI.${STRAND}.5p.wig.tmp
		IF=${CSCORE}/${PROTEIN}_${COND}_IF.${STRAND}.5p.wig.tmp

		OUTPUT=${CSCORE}/${PROTEIN}_${COND}.${STRAND}.Cscore.bed

		# Calculate C-score
		paste $UNTX $A $TI $IF | cut -f1,2,4,6,8 | $SCRIPTS/Cscorelog.pl - | awk 'BEGIN {FS="\t"; OFS="\t"} {if (($1-2) == -1){print "chrM",16298, 16299,$6} else{print "chrM",$1-2,$1-1,$6} }' | sort -k2,2n > $OUTPUT

		# Convert C-score BED to wig format
		WIG=${CSCORE}/${PROTEIN}_${COND}.${STRAND}.Cscore.wig
		cut -f3 $OUTPUT > $WIG
		sed -i '1ivariableStep\tchrom=chrM' $WIG
		sed -i '1itrack type=wiggle_0 name="'$WIG'" description="'$WIG'" visibility=full color=0,0,0 autoScale=on' $WIG

		rm $UNTX $A $TI $IF
	done
done

############################
# Calculate C-score cut-offs
############################
mkdir $CUTOFF

for COND in WT KO
do
	FWD=${CSCORE}/${PROTEIN}_${COND}.fwd.Cscore.bed
	REV=${CSCORE}/${PROTEIN}_${COND}.rev.Cscore.bed
	OUTPUT=${CUTOFF}/${PROTEIN}_${COND}.Cscore.bed.stat
	CUTOFFOUT=${CUTOFF}/${PROTEIN}_${COND}.Cscore.bed.cutoff

	cat $FWD $REV | cut -f4 > $OUTPUT
	Rscript ${SCRIPTS}/Cscorecutoff.R $OUTPUT $CUTOFFOUT
done

