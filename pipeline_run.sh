#!/bin/bash

USAGEMSG="usage: $(basename $0) [-p protein name] [-d input data directory] [-w working directory] [-m min kength] [-M max length] [-f flank length] [-i iterations]

RNase footprinting pipeline. This will run the complete pipeline (scripts 1 - 3).

-p protein name			The name of the protein of interest.
-d data directory		The input data directory that contains the 5' end normalised wig files 
				(required input name format: PROTEIN_COND_TREAT.STRAND.5p.wig where:
				PROTEIN = name of protein of interest
				COND = WT/KO
				TREAT = UNTX/A/TI/IF
				STRAND = fwd/rev
-w working directory		The working directory where all pipeline outputs will go
-m minlength			The minimum footprint length
-M maxlength			The maximum footprint length
-f flank			The flanking length
-i iterations			The number of iterations for cscore shuffling (e.g. 1000)

The pipeline output will be stored in a directory in the specified working directory.
"

# This command prints usage and the name of the script
[ $# -lt 5 ] && echo "$USAGEMSG" >&2 && exit 1

while getopts "p:d:w:m:M:f:i:" flag
do
	case ${flag} in
	(p) PROTEIN=${OPTARG};;
	(d) DATA=${OPTARG};;
	(w) WD=${OPTARG};;
	(m) MIN=${OPTARG};;
	(M) MAX=${OPTARG};;
	(f) FLANK=${OPTARG};;
	(i) ITERATIONS=${OPTARG};;
	(*) echo "$0: error - unrecognized option $1" >&2; exit 1;;
	esac
done

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

mkdir -p $WD

echo 'Calculating C-scores & cutoffs'
$SCRIPT_DIR/1.calculate_cscore.sh -p $PROTEIN -d $DATA -w $WD || echo 'error' && exit
echo 'Calling footprints & calculating F-scores'
$SCRIPT_DIR/2.call_footprints.sh -p $PROTEIN -w $WD -m $MIN -M $MAX -f $FLANK || echo 'error' && exit
echo 'Creating empirical null distribution & calculating FDR'
$SCRIPT_DIR/3.null_dist_FDR.sh -p $PROTEIN -w $WD -m $MIN -M $MAX -f $FLANK -i $ITERATIONS || echo 'error' && exit
echo 'RNase footprinting complete'


