# PacBio de novo HiFi assembly and QA workflow

A reproducible Snakemake workflow for the de novo assembly and quality assessment of PacBio HiFi reads. 

## Getting going

This workflow is set up for execution on a Linux HPC with SLURM scheduling. It also requires an installation of Snakemake, Conda and Singularity. To run, follow these steps:

1. Clone the repo
2. Install tools as detailed below. Some additional configuration may be required to interface with a different HPC.
3. Create a folder called `reference` in the base directory and save the reference genome of interest in there under the name `genome.fasta`. Repeat for mitochondria and chloroplast references, naming them `mitochondria.fasta` and `chloroplast.fasta`. 
4. Create a folder called `0_raw` that contains the raw HiFi reads in .subreads.bam format. 
5. Enter the prefix of this file into the `PREFIXES` variable at the very top of the Snakefile.
6. Activate Snakemake
7. Enter the desired file path (found in the `rule all:` section) with wildcards substituted for actual values 
8. Multiple values can be specified for each wildcard

eg. To generate HiFi reads using CCS processing with 40 chunks, the following command would suffice

`snakemake --profile profiles/slurm --use-conda --use-singularity 1_subreads/chunks/SRR9969479.m64011_190430_000122.ccs.40.bam`

The workflow joins all files generated from the prefixes in the `PREFIXES` parameter into a single file with the descriptor specified in the `SAMPLES` parameter. To run quast on a 30x randomly subset hifiasm genome assembly, the following would suffice:

`snakemake --profile profiles/slurm --use-conda --use-singularity quast/assemblytype_genome_assemblytool_hifiasm_readselect_random_prefix_SRR99694_kmer_31_cov_0.9_depth_30/report.tsv`

Multiple different files can also be produced at once by specifying multiple filepaths.

To see which files will be produced when a call is made, specify the --dryrun option.

The workflow is set up to assemble nuclear genomes, chloroplast and mitochondria so `genome` must be specified as the {ASS_TYPE} wildcard

## Installation of tools.
The majority of tools used are available via conda environments or as docker images. The following tools were not and require manual installation. The exact commands used are documented. The tools were all installed in a /tools folder on the same level as the pacbio_hifi/ directory. The paths to the executables should be hard coded into variables at the beginning of the Snakefile. 

### Hifiasm
```
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make
```

### Mummer v4.0
```
wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
tar -xf mummer-4.0.0beta2.tar.gz
rm mummer-4.0.0beta2.tar.gz
cd mummer-4.0.0beta2/
./configure
make
```

### LTR_FINDER_Parallel
```
git clone https://github.com/oushujun/LTR_FINDER_parallel.git
```

### Hicanu
```
git clone https://github.com/marbl/canu.git
cd canu
git checkout hicanu_rc
cd src
make
```

## Still ToDo

- select_longest doesn't fail when a coverage greater than is available is requested
- streamline merging of chunked ccs from separate files into a single rule
- when kmer coverage and kmer length are 0, don't extract chloroplast and mitochondrial reads from nuclear genome dataset
