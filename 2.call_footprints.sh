#!/bin/bash

USAGEMSG="usage: $(basename $0) [-p protein name] [-w working directory] [-m min length] [-M max length] [-f flank length]

Call footprints for WT and KO samples.

-p protein name			The name of the protein of interest.
-w working directory		The working directory where all outputs will go
-m minlength			The minimum footprint length
-M maxlength			The maximum footprint length
-f flank			The flanking length

The output will be stored in a directories named footprints and comparison in the working directory. BEDtools <2.27.0 must be in PATH (2.25.0 & 2.26.0 tested).
"

# This command prints usage and the name of the script
[ $# -lt 5 ] && echo "$USAGEMSG" >&2 && exit 1

while getopts "p:w:m:M:f:" flag
do
	case ${flag} in
	(p) PROTEIN=${OPTARG};;
	(w) WD=${OPTARG};;
	(m) MIN=${OPTARG};;
	(M) MAX=${OPTARG};;
	(f) FLANK=${OPTARG};;
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
mkdir $FOOTPRINTS

############################
# Call footprints
############################
for COND in WT KO
do
	CUTOFFIN=${CUTOFF}/${PROTEIN}_${COND}.Cscore.bed.cutoff
	FWD=${CSCORE}/${PROTEIN}_${COND}.fwd.Cscore.bed
	REV=${CSCORE}/${PROTEIN}_${COND}.rev.Cscore.bed

	TMP=${FOOTPRINTS}/${PROTEIN}_${COND}.footprints.bed.tmp
	OUTPUT=${FOOTPRINTS}/${PROTEIN}_${COND}.footprints.bed

	# Forward strand
	cat $CUTOFFIN | uniq | cat - | while read line; do ${SCRIPTS}/Rnasecutoff.pl $FWD $line | sort -k2,2n | ${SCRIPTS}/Rnasefootprintlog.pl - $FWD $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'+'", $5,$7, $8}' >> $TMP; done

	# Reverse strand
	cat $CUTOFFIN | uniq | cat - | while read line; do ${SCRIPTS}/Rnasecutoff.pl $REV $line | sort -k2,2n | ${SCRIPTS}/Rnasefootprintlog.pl - $REV $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'-'", $5,$7, $8}' >> $TMP; done

	sort -k2,2n $TMP | uniq > $OUTPUT
	rm $TMP


	TMP=${FOOTPRINTS}/${PROTEIN}_${COND}.footprints.merged
	MERGED=${FOOTPRINTS}/${PROTEIN}_${COND}.footprints.merged.bed
	TMP2=${FOOTPRINTS}/${PROTEIN}_${COND}.footprints.merged.tmp

	# Reformat & merge footprints
	# Version of bedtools MUST be earlier than 2.27.0 (2.25.0 & 2.26.0 tested) as changes to the output were made that breaks this portion of the script
	awk 'BEGIN{FS="\t";OFS="\t"}{if($3>$2){print $0}}' $OUTPUT | sort -k2,2n - > $TMP
	$BEDTOOLS merge -s -c 1 -o count -i $TMP |  awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$3,"I",$5,$4}' - | sed 's/[ \t]*$//' | $BEDTOOLS intersect -s -wa -wb -a - -b $TMP > $MERGED

	CALLED=${FOOTPRINTS}/${PROTEIN}_${COND}.called_footprints.bed
	# Pick footprints with smallest F-score
	awk 'BEGIN{FS="\t"; OFS="\t"} NR == FNR { ID=$1 $2 $3 $4; if(array[ID] == ""){ array[ID]=$11} else {if (array[ID] > $11){array[ID]=$11}} } NR > FNR {ID=$1 $2 $3 $4; if (array[ID] == $11 ){print $7,$8,$9,$10,$11,$12,$13,$14,$15}}' $MERGED $MERGED > $TMP2

	# Fix for circular chrM
	awk 'BEGIN{FS="\t";OFS="\t"}{if($3<$2 && $6 == "+"){print $0}}' $OUTPUT | sort -k5,5n | head -1 >> $TMP2
	awk 'BEGIN{FS="\t";OFS="\t"}{if($3<$2 && $6 == "-"){print $0}}' $OUTPUT | sort -k5,5n | head -1 >> $TMP2
	sort -k2,2n -o $CALLED $TMP2

	rm $TMP $TMP2 $MERGED
done

############################
# Compare footprints
############################
mkdir $COMPARE

for COND in WT KO
do
	INPUT=${FOOTPRINTS}/${PROTEIN}_${COND}.called_footprints.bed
	FWD=${FOOTPRINTS}/${PROTEIN}_${COND}.called_footprints.fwd.bed
	REV=${FOOTPRINTS}/${PROTEIN}_${COND}.called_footprints.rev.bed

	# Split footprints by strand
	awk 'BEGIN{FS="\t"; OFS="\t"} {if ($6=="+" && $2>$3) {print $1,$2,$3,(16299-$2)+$3} else if ($6=="+" && $3>$2) {print $1,$2,$3,$3-$2} else {}}' $INPUT > $FWD
	awk 'BEGIN{FS="\t"; OFS="\t"} {if ($6=="-" && $2>$3) {print $1,$2,$3,(16299-$2)+$3} else if ($6=="-" && $3>$2) {print $1,$2,$3,$3-$2} else {}}' $INPUT > $REV
done

# WT called footprints
WT_FWD=${FOOTPRINTS}/${PROTEIN}_WT.called_footprints.fwd.bed
WT_REV=${FOOTPRINTS}/${PROTEIN}_WT.called_footprints.rev.bed
WT_FWD_TMP=${COMPARE}/${PROTEIN}_WT.footprint_Fscores.fwd.tmp
WT_REV_TMP=${COMPARE}/${PROTEIN}_WT.footprint_Fscores.rev.tmp
WT_MERGE=${COMPARE}/${PROTEIN}_WT.footprint_Fscores.bed

KO_C_FWD=${CSCORE}/${PROTEIN}_KO.fwd.Cscore.bed
KO_C_REV=${CSCORE}/${PROTEIN}_KO.rev.Cscore.bed

# Retrieve F-scores for equivalent regions in KO data set
${SCRIPTS}/footprint_compare.pl $WT_FWD $KO_C_FWD $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'+'", $5,$7,$8}' > $WT_FWD_TMP
${SCRIPTS}/footprint_compare.pl $WT_REV $KO_C_REV $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'-'", $5,$7,$8}' > $WT_REV_TMP

cat $WT_FWD_TMP $WT_REV_TMP > $WT_MERGE
sort -k2,2n -o $WT_MERGE $WT_MERGE
rm $WT_FWD_TMP $WT_REV_TMP


# KO called footprints
KO_FWD=${FOOTPRINTS}/${PROTEIN}_KO.called_footprints.fwd.bed
KO_REV=${FOOTPRINTS}/${PROTEIN}_KO.called_footprints.rev.bed
KO_FWD_TMP=${COMPARE}/${PROTEIN}_KO.footprint_Fscores.fwd.tmp
KO_REV_TMP=${COMPARE}/${PROTEIN}_KO.footprint_Fscores.rev.tmp
KO_MERGE=${COMPARE}/${PROTEIN}_KO.footprint_Fscores.bed

WT_C_FWD=${CSCORE}/${PROTEIN}_WT.fwd.Cscore.bed
WT_C_REV=${CSCORE}/${PROTEIN}_WT.rev.Cscore.bed

# Retrieve F-scores for equivalent regions in WT data set
${SCRIPTS}/footprint_compare.pl $KO_FWD $WT_C_FWD $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'+'", $5,$7,$8}' > $KO_FWD_TMP
${SCRIPTS}/footprint_compare.pl $KO_REV $WT_C_REV $MIN $MAX $FLANK | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,"footprints_" $4 "_" $5 "_" $6, $6, "'-'", $5,$7,$8}' > $KO_REV_TMP

cat $KO_FWD_TMP $KO_REV_TMP > $KO_MERGE
sort -k2,2n -o $KO_MERGE $KO_MERGE
rm $KO_FWD_TMP $KO_REV_TMP


# Create footprint data arrays
# WT called
paste ${FOOTPRINTS}/${PROTEIN}_WT.called_footprints.bed $WT_MERGE | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$14,$16,log($14/$5)/log(2)}' > ${COMPARE}/${PROTEIN}.WT_called.footprint_scores.bed
sed -i '1iChr\tStart\tEnd\tName\tFscore\tStrand\tCentreCscore\tKOFscore\tKOCentreCscore\tlog2FC_Fscore' ${COMPARE}/${PROTEIN}.WT_called.footprint_scores.bed

# KO called
paste ${FOOTPRINTS}/${PROTEIN}_KO.called_footprints.bed $KO_MERGE | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$14,$16,log($5/$14)/log(2)}' > ${COMPARE}/${PROTEIN}.KO_called.footprint_scores.bed
sed -i '1iChr\tStart\tEnd\tName\tFscore\tStrand\tCentreCscore\tWTFscore\tWTCentreCscore\tlog2FC_Fscore' ${COMPARE}/${PROTEIN}.KO_called.footprint_scores.bed

