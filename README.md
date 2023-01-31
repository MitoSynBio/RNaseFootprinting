# RNase Footprinting Pipeline

## Run whole RNase footprinting pipeline
Usage: `./pipeline_run.sh [-p protein name] [-d input data directory] [-w working directory] [-m min kength] [-M max length] [-f flank length] [-i iterations]`

RNase footprinting pipeline. This will run the complete pipeline (scripts 1 - ). 

`-p protein name` - The name of the protein of interest.<br>
`-d data directory` - The input data directory that contains the 5' end normalised wig files<br>
`-w working directory` - The working directory where all outputs will<br>
`-m minlength` - The minimum footprint length<br>
`-M maxlength` - The maximum footprint length<br>
`-f flank` - The flanking length<br>
`-i iterations` - The number of iterations for cscore shuffling (e.g. 1000)<br>

Required input name format: PROTEIN_COND_TREAT.STRAND.5p.wig, where:<br>
PROTEIN = name of protein of interest<br>
COND = WT/KO<br>
TREAT = UNTX/A/TI/IF<br>
STRAND = fwd/rev<br>

The pipeline output will be stored in a directory in the specified working directory.
<br>
<br>
<br>
#### OR, run the pipeline in sections:
## 1. Calculate C-scores
Usage: `./1.calculate_cscore.sh [-p protein name] [-d input data directory] [-w working directory]`

Calculate C-score for WT and KO samples.

`-p protein name` - The name of the protein of interest.<br>
`-d data directory` - The input data directory that contains the 5' end normalised wig files<br>
`-w working directory` - The working directory where all outputs will go<br>

Required input name format: PROTEIN_COND_TREAT.STRAND.5p.wig, where:<br>
PROTEIN = name of protein of interest<br>
COND = WT/KO<br>
TREAT = UNTX/A/TI/IF<br>
STRAND = fwd/rev<br>

The output will be stored in a directory named **cscore** in the working directory.
<br>
<br>

## 2. Call footprints

Usage: `./2.call_footprints.sh [-p protein name] [-w working directory] [-m min length] [-M max length] [-f flank length]`

Call footprints for WT and KO samples.

`-p protein name` - The name of the protein of interest.<br>
`-w working directory` - The working directory where all outputs will go<br>
`-m minlength` - The minimum footprint length<br>
`-M maxlength` - The maximum footprint length<br>
`-f flank` - The flanking length<br>

The output will be stored in a directory named **footprints** in the working directory. BEDtools <2.27.0 must be in $PATH (2.25.0 & 2.26.0 tested).<br>
<br>
<br>

## 3. Calculate FDR

Create empirical null distribution and calculate FDR.

Usage: .`/3.null_dist_FDR.sh [-p protein name] [-w working directory] [-m min length] [-M max length] [-f flank length] [-i iterations]`

`-p protein name` - The name of the protein of interest.<br>
`-w working directory` - The working directory where all outputs will go<br>
`-m minlength` - The minimum footprint length<br>
`-M maxlength` - The maximum footprint length<br>
`-f flank` - The flanking length<br>
`-i iterations` - The number of iterations for cscore shuffling (e.g. 1000)<br>

The output will be stored in a directory named **null_dist** in the working directory.<br>

