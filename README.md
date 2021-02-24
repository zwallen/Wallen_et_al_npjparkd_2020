This github repository houses scripts used in the processing of 16S rRNA DNA sequences and statistical analysis detailed in 
**Wallen et al. Characterizing dysbiosis of gut microbiome in PD: Evidence for overabundance of opportunistic pathogens. *npj Parkinson's Disease* 6, Article number: 11 (2020).doi: 10.1038/s41531-020-0112-6.**

The following README gives an overview of the overall structure of the repository, and important notes on how to run the scripts.

## Directory tree for repository
```
Wallen_et_al_npjparkd_2020
|
|-- Dataset1_Scripts -- houses scripts used to process sequences and run statistical analyses for Dataset 1
|   |
|   |-- 0.Bioinformatics -- scripts used to process sequences from raw sequences to ASV table, ASV taxonomy assignments
|                            and phylogenetic trees
|
|-- Dataset2_Scripts -- houses scripts used to process sequences and run statistical analyses for Dataset 2
|   |
|   |-- 0.Bioinformatics -- scripts used to process sequences from raw sequences to ASV table, ASV taxonomy assignments
|                           and phylogenetic trees
|
|-- Joint_Analyses_Scripts -- houses scripts used to generate plots and run analyses that involved both datasets
|
|-- Metadata -- sample metadata for Dataset 1 and 2, used in creation of phyloseq objects
|  
|-- PhyloseqObjects -- phyloseq objects used in analyses that were reported in manuscript
|   |
|   |-- Dataset1 -- phyloseq object for Dataset 1
|   |
|   |-- Dataset2 -- phyloseq object for Dataset 2
|
|-- Script_Output -- directory to house output from bioinformatic and analysis scripts
|   |
|   |-- Dataset1_Output -- output from scripts in Dataset1_Scripts directory
|   |
|   |-- Dataset2_Output -- output from scripts in Dataset2_Scripts directory
|   |
|   |-- Joint_Analyses_Output -- output from scripts in Joint_Analyses_Scripts directory
|
|-- Sequences -- houses a script to download sequences from SRA, and then the raw 16S rRNA DNA sequences once downloaded
|   |
|   |-- Dataset1 -- created by script that downloads SRA sequences. Houses raw 16S rRNA DNA sequences for dataset 1
|   |
|   |-- Dataset2 -- created by script that downloads SRA sequences. Houses raw 16S rRNA DNA sequences for dataset 2
|       |
|       |-- 5176-HP-Pool_01 -- houses sequences for dataset 2 samples included in sequencing pool #1
|       |
|       |-- 5176-HP-Pool_02 -- houses sequences for dataset 2 samples included in sequencing pool #2
|       |
|       |-- 5176-HP-Pool_03 -- houses sequences for dataset 2 samples included in sequencing pool #3
|       |
|       |-- 5176-HP-Pool_04 -- houses sequences for dataset 2 samples included in sequencing pool #4
|       |
|       |-- 5176-HP-Pool_05 -- houses sequences for dataset 2 samples included in sequencing pool #5
|       |
|       |-- 5176-HP-Pool_06 -- houses sequences for dataset 2 samples included in sequencing pool #6
|
|-- Support_Files -- files called upon by certain bioinformatic and analyses scripts
```

## Important notes about this repository

#### The repository is structured to work out of the box
Once downloaded, all bioinformatic and analyses scripts should be able to be run as is, without any modification to the repository structure or scripts themselves. The only prerequisites is that raw sequences need to be downloaded from SRA before running the bioinformatic scripts. A script titled `Download_SRA_Sequences.sh` is located in the `Sequences/` directory to do this.

#### Phyloseq objects used in the manuscript are included in this repository
As stated in the directory tree, phyloseq objects used in the manuscript for datasets 1 and 2 are located in the `PhyloseqObjects/` directory. Running scripts in `Dataset1_Scripts/`, `Dataset2_Scripts/`, and `Joint_Analyses_Scripts/` directories using these phyloseq objects should give same results as reported in the manuscript (with the exception of analyses that use permutations or iterations, i.e. PERMANOVA and SparCC, but results should be very similar). 

*WARNING: As the scripts are currently set up, the default phyloseq objects that are included in this repository will be overwritten if re-running bioinformatics on raw sequences without renaming the phyloseq objects first.*

#### Results for analyses utilizing manuscript phyloseq objects are included in this repository
Results that are outputted by analyses scripts and that were reported in the manuscript are located in the `Script_Output/` directory. The only results that are not included in this directory are the visualizations of correlation networks (as these were not generated using a script, but using the GUI of Gephi), but unaltered versions of these figures can all be found in the supplementary material of the manuscript.

*WARNING: As the scripts are currently set up, the results that are included in the `Script_Output/` repository will be overwritten if re-running the analyses scripts without moving or renaming the result files first.*

#### There are three types of scripts stored in this repository
Scripts with the extension `.sh` are shell scripts that were used to implement bioinformatic command line programs. Scripts with the extension `.R` are R scripts written in R programming language used to perform bioinformatics, statistical analyses, and generate plots. Scripts that have a `.job` extension are shell scripts that were used to submit `.sh` or `.R` scripts to a SLURM scheduling system on a high performance computing cluster. Each `.sh` and `.R` script should have an accompanying `.job` script.

*Note: `.sh` and `.R` can be run without using the `.job` script*

#### Scripts in `Dataset1_Scripts/` and `Dataset2_Scripts/` directories are named to reflect the approximate order they were implemented in the study and what type of action or analysis they perform
Each script name begins with a number (to signify the order in which they were ran) followed by a title that gives an idea of what action/analysis each script performs. This naming format also applies to the scripts in the `0.Bioinformatics/` subdirectories. Each script contains a comment header that should give a little more detail about what action/analysis a script performs if that cannot be gleaned from the script file name.

*Note: the scripts in `Joint_Analyses_Scripts/` do not follow this naming format, they can be ran in any order*

#### Some directories contain `ErrorOut/` and `Output/` subdirectories
These directories are where the stderr (`ErrorOut/`) and stdout (`Output/`) are piped to when running scripts using a SLURM scheduler. README files are located in these directories, but only as place holders so they show up in the repository.
