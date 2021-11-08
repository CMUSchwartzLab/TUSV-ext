# TUSV-ext: Reconstructing tumor phylogenetic tree with SVs, CNVs and SNVs

### Requirements

The program is written in Python2. Python3 version will come out soon.

The packages you will need to install are listed below.

* `numpy`
* `graphviz`
* `ete2`
* `biopython`
* `gurobipy`

### Gurobi License

To obtain a Gurobi license, you can sign up as an academic user here [http://www.gurobi.com/registration/academic-license-reg](http://www.gurobi.com/registration/academic-license-reg) and follow the instructions for downloading a license. 

## Running

* tusv.py
* multi_tusv.py

The script `tusv.py` takes as input a single directory containing one or multiple `.vcf` files. Go here [https://samtools.github.io/hts-specs/VCFv4.2.pdf](https://samtools.github.io/hts-specs/VCFv4.2.pdf) for specifications on the `.vcf` format. Sample input can be found in `tusv/data/example/patient1/`. Each `.vcf` file corresponds to a patient sample and all records inside each `.vcf` file contain the structural variants and copy number segments for that sample. The script `multi_tusv.py` is a wrapper script of `tusv.py` and takes as input a directory where each subdirectory contains one or multiple `.vcf` files. Both scripts take required inputs

* `-n` number of leaves to infer in phylogenetic tree
* `-c` maximum copy number allowed for any breakpoint or segment on any node
* `-t` maximum number of coordinate-descent iterations (program can finish sooner if convergence is reached)
* `-r` number of random initializations of the coordinate-descent algorithm

and optional inputs

* `-b` (recommended) binary flag to automatically set hyper-parameters lambda 1 and lambda 2
* `-l` lambda 1 hyper-parameter. controls phylogenetic cost
* `-a` lambda 2 hyper-parameter. controls breakpoint to segment consistancy
* `-m` maximum time (in seconds) for a single cordinate-descent iteration
* `-s` number of segments (in addition to those containing breakpoints) that are randomly kept for unmixing. default keeps all segments
* `-p` (not recommended) number of processors to use. uses all available processors by default
* `-d` (not recommended) file containing metadata information for output .vcf files

