#!/bin/bash

USAGEMSG="usage: $(basename $0) [-p protein name] [-w working directory] [-m min length] [-M max length] [-f flank length] [-i iterations]

Create empirical null distribution and calculate FDR.

-p protein name			The name of the protein of interest. 
-w working directory		The working directory where all outputs will go
-m minlength			The minimum footprint length
-M maxlength			The maximum footprint length
-f flank			The flanking length
-i iterations			The number of iterations for cscore shuffling (e.g. 1000)

The output will be stored in a directory named null_dist in the working directory. BEDtools <2.27.0 must be in PATH (2.25.0 & 2.26.0 tested).
"

# This command prints usage and the name of the script
[ $# -lt 5 ] && echo "$USAGEMSG" >&2 && exit 1

while getopts "p:w:m:M:f:i:" flag
do
	case ${flag} in
	(p) PROTEIN=${OPTARG};;
	(w) WD=${OPTARG};;
	(m) MIN=${OPTARG};;
	(M) MAX=${OPTARG};;
	(f) FLANK=${OPTARG};;
	(i) ITERATIONS=${OPTARG};;
	(*) echo "$0: error - unrecognized option $1" >&2; exit 1;;
	esac
done

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Test bedtools
#BEDTOOLS=/media/Data/tools/bedtools/bedtools-2.26.0/bin/bedtools
BEDTOOLS=`which bedtools`

## Test for BEDtools version
currBT="$($BEDTOOLS --version | sed 's/bedtools v//g')"
maxBT="2.26.0"

if [ "$(printf '%s/n' "$maxBT" "$currBT" | sort -V | tail -n1)" = "$maxBT" ]; then
	echo "BEDtools is <= 2.26.0."
else
	echo "BEDtools > 2.26.0. Please use BEDtools <= 2.26.0."
	exit 1
fi

SCRIPTS=$SCRIPT_DIR/scripts
CSCORE=$WD/cscore
CUTOFF=$WD/cutoff
FOOTPRINTS=$WD/footprints
COMPARE=$WD/comparison
NULL=$WD/null_dist
mkdir $NULL

####################################################################################
### Construction of empirical null model for footprinting analysis ###
####################################################################################

## Cscore shuffling both KO & WT samples
for COND in WT KO
do
	OUTDIR=$NULL/${COND}_null/cscore

	mkdir -p $OUTDIR

	for STRAND in fwd rev
	do
		INPUT=${CSCORE}/${PROTEIN}_${COND}.${STRAND}.Cscore.bed

		for ((i=1;i<=ITERATIONS;i++))
		do
			OUTPUT=${OUTDIR}/${PROTEIN}_${COND}.${STRAND}.Cscore.shuf_$i.bed
			Rscript $SCRIPTS/Cscore_shuffle.R $INPUT $OUTPUT
		done
	done
done

## Strand splitting - Experimental Fscores
for COND in WT KO
do
	INPUT=$COMPARE/${PROTEIN}_${COND}.footprint_Fscores.bed
	OUT_FWD=$COMPARE/${PROTEIN}_${COND}.footprint_Fscores.fwd.bed
	OUT_REV=$COMPARE/${PROTEIN}_${COND}.footprint_Fscores.rev.bed

	awk 'BEGIN{FS="\t"; OFS="\t"} {if ($6=="+" && $2>$3) {print $1,$2,$3,(16299-$2)+$3} else if ($6=="+" && $3>$2) {print $1,$2,$3,$3-$2} else {}}' $INPUT > $OUT_FWD
	awk 'BEGIN{FS="\t"; OFS="\t"} {if ($6=="-" && $2>$3) {print $1,$2,$3,(16299-$2)+$3} else if ($6=="-" && $3>$2) {print $1,$2,$3,$3-$2} else {}}' $INPUT > $OUT_REV
done


##############################################
## Null for footprints called on WT samples ##
##############################################

DIST=${NULL}/WT_null/dists
mkdir -p $DIST
mkdir ${NULL}/WT_null/comparison

### FWD ###
awk -v DIST="$DIST" -v PROTEIN="$PROTEIN" '{filename = sprintf(DIST"/"PROTEIN"_WT.footprints.fwd.null_%d.bed", NR); print >filename; close(filename)}' $COMPARE/${PROTEIN}_WT.footprint_Fscores.fwd.bed

FWD_LEN=`ls $DIST/${PROTEIN}_WT.footprints.fwd.null_*.bed | wc -l`

for ((i=1;i<=FWD_LEN;i++))
do
	for ((j=1;j<=ITERATIONS;j++))
	do
		INPUT=${DIST}/${PROTEIN}_WT.footprints.fwd.null_$i.bed
		WT_C_FWD=${NULL}/WT_null/cscore/${PROTEIN}_WT.fwd.Cscore.shuf_$j.bed
		KO_C_FWD=${NULL}/KO_null/cscore/${PROTEIN}_KO.fwd.Cscore.shuf_$j.bed
		WT_FWD_TMP=${DIST}/${PROTEIN}_WT.footprints.fwd.null_$i.WT_cscores.tmp
		KO_FWD_TMP=${DIST}/${PROTEIN}_WT.footprints.fwd.null_$i.KO_cscores.tmp

		# For each experimental footprint, extract Cscore from KO & WT from shuffled Cscore
		${SCRIPTS}/footprint_compare.pl $INPUT $WT_C_FWD $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'+'", $5,$7,$8}' >> $WT_FWD_TMP
		${SCRIPTS}/footprint_compare.pl $INPUT $KO_C_FWD $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'+'", $5,$7,$8}' >> $KO_FWD_TMP
	done
	OUTPUT=${NULL}/WT_null/comparison/${PROTEIN}_WT.footprints.fwd.null_$i.log2FC.bed

	paste $WT_FWD_TMP $KO_FWD_TMP | awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$3,$4,$5,$6,$14,log($14/$5)/log(2)}' - > $OUTPUT
	rm $WT_FWD_TMP $KO_FWD_TMP
done


### REV ###
awk -v DIST="$DIST" -v PROTEIN="$PROTEIN" '{filename = sprintf(DIST"/"PROTEIN"_WT.footprints.rev.null_%d.bed", NR); print >filename; close(filename)}' $COMPARE/${PROTEIN}_WT.footprint_Fscores.rev.bed

REV_LEN=`ls $DIST/${PROTEIN}_WT.footprints.rev.null_*.bed | wc -l`

for ((i=1;i<=REV_LEN;i++))
do
	for ((j=1;j<=ITERATIONS;j++))
	do
		INPUT=${DIST}/${PROTEIN}_WT.footprints.rev.null_$i.bed
		WT_C_REV=${NULL}/WT_null/cscore/${PROTEIN}_WT.rev.Cscore.shuf_$j.bed
		KO_C_REV=${NULL}/KO_null/cscore/${PROTEIN}_KO.rev.Cscore.shuf_$j.bed
		WT_REV_TMP=${DIST}/${PROTEIN}_WT.footprints.rev.null_$i.WT_cscores.tmp
		KO_REV_TMP=${DIST}/${PROTEIN}_WT.footprints.rev.null_$i.KO_cscores.tmp

		# For each experimental footprint, extract Cscore from KO & WT from shuffled Cscore
		${SCRIPTS}/footprint_compare.pl $INPUT $WT_C_REV $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'-'", $5,$7,$8}' >> $WT_REV_TMP
		${SCRIPTS}/footprint_compare.pl $INPUT $KO_C_REV $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'-'", $5,$7,$8}' >> $KO_REV_TMP
	done
	OUTPUT=${NULL}/WT_null/comparison/${PROTEIN}_WT.footprints.rev.null_$i.log2FC.bed

	paste $WT_REV_TMP $KO_REV_TMP | awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$3,$4,$5,$6,$14,log($14/$5)/log(2)}' - > $OUTPUT
	rm $WT_FWD_TMP $KO_FWD_TMP
done


##############################################
## Null for footprints called on KO samples ##
##############################################

DIST=${NULL}/KO_null/dists
mkdir -p $DIST
mkdir ${NULL}/KO_null/comparison

### FWD ###
awk -v DIST="$DIST" -v PROTEIN="$PROTEIN" '{filename = sprintf(DIST"/"PROTEIN"_KO.footprints.fwd.null_%d.bed", NR); print >filename; close(filename)}' $COMPARE/${PROTEIN}_KO.footprint_Fscores.fwd.bed

FWD_LEN=`ls $DIST/${PROTEIN}_KO.footprints.fwd.null_*.bed | wc -l`

for ((i=1;i<=FWD_LEN;i++))
do
	for ((j=1;j<=ITERATIONS;j++))
	do
		INPUT=${DIST}/${PROTEIN}_KO.footprints.fwd.null_$i.bed
		WT_C_FWD=${NULL}/WT_null/cscore/${PROTEIN}_WT.fwd.Cscore.shuf_$j.bed
		KO_C_FWD=${NULL}/KO_null/cscore/${PROTEIN}_KO.fwd.Cscore.shuf_$j.bed
		WT_FWD_TMP=${DIST}/${PROTEIN}_KO.footprints.fwd.null_$i.WT_cscores.tmp
		KO_FWD_TMP=${DIST}/${PROTEIN}_KO.footprints.fwd.null_$i.KO_cscores.tmp

		# For each experimental footprint, extract Cscore from KO & WT from shuffled Cscore
		${SCRIPTS}/footprint_compare.pl $INPUT $WT_C_FWD $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'+'", $5,$7,$8}' >> $WT_FWD_TMP
		${SCRIPTS}/footprint_compare.pl $INPUT $KO_C_FWD $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'+'", $5,$7,$8}' >> $KO_FWD_TMP
	done
	OUTPUT=${NULL}/KO_null/comparison/${PROTEIN}_KO.footprints.fwd.null_$i.log2FC.bed

	paste $WT_FWD_TMP $KO_FWD_TMP | awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$3,$4,$5,$6,$14,log($14/$5)/log(2)}' - > $OUTPUT
	rm $WT_FWD_TMP $KO_FWD_TMP
done


### REV ###
awk -v DIST="$DIST" -v PROTEIN="$PROTEIN" '{filename = sprintf(DIST"/"PROTEIN"_KO.footprints.rev.null_%d.bed", NR); print >filename; close(filename)}' $COMPARE/${PROTEIN}_KO.footprint_Fscores.rev.bed

REV_LEN=`ls $DIST/${PROTEIN}_KO.footprints.rev.null_*.bed | wc -l`

for ((i=1;i<=REV_LEN;i++))
do
	for ((j=1;j<=ITERATIONS;j++))
	do
		INPUT=${DIST}/${PROTEIN}_KO.footprints.rev.null_$i.bed
		WT_C_REV=${NULL}/WT_null/cscore/${PROTEIN}_WT.rev.Cscore.shuf_$j.bed
		KO_C_REV=${NULL}/KO_null/cscore/${PROTEIN}_KO.rev.Cscore.shuf_$j.bed
		WT_REV_TMP=${DIST}/${PROTEIN}_KO.footprints.rev.null_$i.WT_cscores.tmp
		KO_REV_TMP=${DIST}/${PROTEIN}_KO.footprints.rev.null_$i.KO_cscores.tmp

		# For each experimental footprint, extract Cscore from KO & WT from shuffled Cscore
		${SCRIPTS}/footprint_compare.pl $INPUT $WT_C_REV $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'-'", $5,$7,$8}' >> $WT_REV_TMP
		${SCRIPTS}/footprint_compare.pl $INPUT $KO_C_REV $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'-'", $5,$7,$8}' >> $KO_REV_TMP
	done
	OUTPUT=${NULL}/KO_null/comparison/${PROTEIN}_KO.footprints.rev.null_$i.log2FC.bed

	paste $WT_REV_TMP $KO_REV_TMP | awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$3,$4,$5,$6,$14,log($14/$5)/log(2)}' - > $OUTPUT
	rm $WT_FWD_TMP $KO_FWD_TMP
done


### Collate ###
## WT Called ##
# FWD #
WT_FWD_LIST=`ls ${NULL}/WT_null/comparison/${PROTEIN}_WT.footprints.fwd.null_*.log2FC.bed`
WT_OUT_FWD=${NULL}/WT_null/comparison/${PROTEIN}_WT.footprints.fwd.null_scores.log2FC.all.bed
cat $WT_FWD_LIST >> $WT_OUT_FWD
sort -k2,2n -o $WT_OUT_FWD $WT_OUT_FWD

# REV #
WT_REV_LIST=`ls ${NULL}/WT_null/comparison/${PROTEIN}_WT.footprints.rev.null_*.log2FC.bed`
WT_OUT_REV=${NULL}/WT_null/comparison/${PROTEIN}_WT.footprints.rev.null_scores.log2FC.all.bed
cat $WT_REV_LIST >> $WT_OUT_REV
sort -k2,2n -o $WT_OUT_REV $WT_OUT_REV

OUT=${NULL}/WT_null/comparison/${PROTEIN}.WT_called.null_scores.all.bed
cat $WT_OUT_FWD $WT_OUT_REV | sort -k2,2n > $OUT
sed -i '1iChr\tStart\tEnd\tName\tFscore\tStrand\tKO_Fscore\tlog2FC_Fscore' $OUT

## KO Called ##
# FWD #
KO_FWD_LIST=`ls ${NULL}/KO_null/comparison/${PROTEIN}_KO.footprints.fwd.null_*.log2FC.bed`
KO_OUT_FWD=${NULL}/KO_null/comparison/${PROTEIN}_KO.footprints.fwd.null_scores.log2FC.all.bed
cat $KO_FWD_LIST >> $KO_OUT_FWD
sort -k2,2n -o $KO_OUT_FWD $KO_OUT_FWD

# REV #
KO_REV_LIST=`ls ${NULL}/KO_null/comparison/${PROTEIN}_KO.footprints.rev.null_*.log2FC.bed`
KO_OUT_REV=${NULL}/KO_null/comparison/${PROTEIN}_KO.footprints.rev.null_scores.log2FC.all.bed
cat $KO_REV_LIST >> $KO_OUT_REV
sort -k2,2n -o $KO_OUT_REV $KO_OUT_REV

OUT=${NULL}/KO_null/comparison/${PROTEIN}.KO_called.null_scores.all.bed
cat $KO_OUT_FWD $KO_OUT_REV | sort -k2,2n > $OUT
sed -i '1iChr\tStart\tEnd\tName\tFscore\tStrand\tKO_Fscore\tlog2FC_Fscore' $OUT


############################
# Calculate FDR
############################
FDR=$WD/FDR
mkdir $FDR

$SCRIPTS/FDR.R $COMPARE/${PROTEIN}.WT_called.footprint_scores.bed ${NULL}/WT_null/comparison/${PROTEIN}.WT_called.null_scores.all.bed $FDR/${PROTEIN}.WT_called.scored_footprints.txt
$SCRIPTS/FDR.R $COMPARE/${PROTEIN}.KO_called.footprint_scores.bed ${NULL}/KO_null/comparison/${PROTEIN}.KO_called.null_scores.all.bed $FDR/${PROTEIN}.KO_called.scored_footprints.txt



