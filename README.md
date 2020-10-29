# tusv - Tree Unmixing with Structural Variants
Model written in Python for deconvolution and phylogeny inference of bulk tumor genomic data using structural variants and copy number segments. This model requires you to obtain a Gurobi license (free for academics).

The pre-print of our paper can be found on bioRXiv here [https://www.biorxiv.org/content/early/2018/01/30/257014](https://www.biorxiv.org/content/early/2018/01/30/257014).

The BibTeX citation for this paper is below.

```
@article {Eaton257014,
  author = {Eaton, Jesse and Wang, Jingyi and Schwartz, Russell},
  title = {Deconvolution and phylogeny inference of structural variations in tumor genomic samples},
  year = {2018},
  doi = {10.1101/257014},
  publisher = {Cold Spring Harbor Laboratory},
  URL = {https://www.biorxiv.org/content/early/2018/01/30/257014},
  eprint = {https://www.biorxiv.org/content/early/2018/01/30/257014.full.pdf},
  journal = {bioRxiv}
}
```

## Setup

### Download

Navigate to a directory where you would like to keep the tusv project. Run the command

```
git clone https://github.com/jaebird123/tusv.git
```

If you do not have the `git` command, download and install it from here [https://git-scm.com/book/en/v2/Getting-Started-Installing-Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

### Install Packages

The packages you will need to install are listed below.

* `numpy`
* `graphviz`
* `ete2`
* `biopython`
* `gurobipy`

To install these, you will need the `pip` command which comes pre-installed with Python. If you do not have the `pip` command, download and install it from here [https://pip.pypa.io/en/stable/installing/](https://pip.pypa.io/en/stable/installing/). Then run the command `pip install <package_name>` for each of the above packages.

### Gurobi License

To obtain a Gurobi license, sign up as an academic user here [http://www.gurobi.com/registration/academic-license-reg](http://www.gurobi.com/registration/academic-license-reg) and follow the instructions for downloading a license. If you are outside of academia, I sincerely apologize as the license is very expensive.

To check that you are able to use Gurobi, after obtaining a license, navigate to the `model/` directory with
```
cd tusv/model/
```
and run the command
```
python test_solver.py
```

If you have successfully installed `gurobipy` and obtained a Gurobi license, the script will output the standard Gurobi solver text and finish executing after approximately 30 seconds.

## Running

Two python scripts can be run.

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

## Data

tusv requires mixed copy numbers of segments and breakpoints for a bulk tumor. We recommend using Weaver which can be found here [https://github.com/ma-compbio/Weaver](https://github.com/ma-compbio/Weaver) before running tusv.

## Authors, License, and Acknowledgments

Authors: Jesse Eaton and Jingyi Wang

License: MIT [https://github.com/jaebird123/tusv/blob/master/LICENSE](https://github.com/jaebird123/tusv/blob/master/LICENSE)

Acknowledgements: We thank Ashok Rajaraman and Jian Ma for helpful discussions and assistance with Weaver. Portions of this work have been funded by the U.S. National Institutes of Health via award R21CA216452 and Pennsylvania Dept. of Health Grant GBMF4554 4100070287. The Pennsylvania Department of Health specifically disclaims responsibility for any analyses, interpretations or conclusions. This work used the Extreme Science and Engineering Discovery Environment (XSEDE), which is supported by National Science Foundation grant number OCI-1053575. Specifically, it used the Bridges system, which is supported by NSF award number ACI-1445606, at the Pittsburgh Supercomputing Center (PSC). The results published here are in whole or part based upon data generated by The Cancer Genome Atlas managed by the NCI and NHGRI. Information about TCGA can be found at [http://cancergenome.nih.gov](http://cancergenome.nih.gov).
