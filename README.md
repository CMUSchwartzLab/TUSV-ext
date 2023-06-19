# TUSV-ext: Reconstructing tumor phylogenetic tree with SVs, CNVs and SNVs

This tool aims to deconvolve multi-regional bulk sequencing samples based on variant calls including SVs, SNVs and CNVs and infer a comprehensive tumor phylogenetic tree.

![pipeline_diagram](https://user-images.githubusercontent.com/44757116/178637052-cdf9be5a-687c-4d13-b8db-e886b506646c.png)

This work is published in Xuecong Fu, Haoyun Lei, Yifeng Tao, Russell Schwartz, Reconstructing tumor clonal lineage trees incorporating single-nucleotide variants, copy number alterations and structural variations, Bioinformatics, Volume 38, Issue Supplement_1, July 2022, Pages i125–i133, https://doi.org/10.1093/bioinformatics/btac253


## Requirements

The program is written in Python2. Python3 version will come out soon.

The packages you will need to install are listed below.

* `numpy`
* `graphviz`
* `ete2`
* `biopython`
* `gurobipy`

### Gurobi License

To obtain a Gurobi license, you can sign up as an academic user here [https://www.gurobi.com/downloads/end-user-license-agreement-academic/](https://www.gurobi.com/downloads/end-user-license-agreement-academic/) and follow the instructions for downloading a license. 

## Running

`python tusv-ext.py`

The script `tusv-ext.py` takes as input a single directory containing one or multiple `.vcf` files. Go here [https://samtools.github.io/hts-specs/VCFv4.2.pdf](https://samtools.github.io/hts-specs/VCFv4.2.pdf) for specifications on the `.vcf` format. Each `.vcf` file should contain SV breakpoints, CNVs and SNVs with their processed copy numbers from one sample of a patient. 

Inputs:
* `-i` the input directory containing vcf files of different samples from one patient
* `-o` the output directory with deconvoluted results
* `-n` number of leaves to infer in phylogenetic tree
* `-c` maximum copy number allowed for any breakpoint or segment on any node
* `-t` maximum number of coordinate-descent iterations (program can finish sooner if convergence is reached)
* `-r` number of random initializations of the coordinate-descent algorithm
* `-col` binary flag whether to collapse the redundant nodes
* `-leaf` binary flag whether to assume only leaf nodes are in the mixed samples or not
* `-sv_ub` approximate maximum number of subsampled breakpoints of structural variants, -1 if you don't want to do the subsampling and include all breakpoints
* `-const` maximum number of total subsampled breakpoints and SNVs

The following inputs are optional:

* `-b` (recommended) binary flag to automatically set hyper-parameters lambda 1 and lambda 2
* `-l` lambda 1 hyper-parameter. controls phylogenetic cost
* `-a` lambda 2 hyper-parameter. controls breakpoint to segment consistancy
* `-m` maximum time (in seconds) for a single cordinate-descent iteration
* `-s` number of segments (in addition to those containing breakpoints) that are randomly kept for unmixing. default keeps all segments
* `-p` (not recommended) number of processors to use. uses all available processors by default

Outputs:
* `C.tsv` the C matrix which is variants copy number profiles of each clone
* `U.tsv` the U matrix which is the frequencies of each clone in each sample
* `T.dot` the inferred phylogenetic tree
