# SurveyFlow
A flexible pipeline for surveying genomic data using fastp, Jellyfish, GenomeScope, and Smudgeplot.

## Project Description
`SurveyFlow` is a Python-based pipeline designed to streamline genomic data analysis by integrating the following tools:
- **fastp**: For quality control of raw sequencing data.
- **Jellyfish**: For k-mer counting.
- **GenomeScope**: For estimating genome size, heterozygosity, and repeat content.
- **Smudgeplot**: For visualizing genome structure and ploidy.

This pipeline is ideal for genomics researchers who need a comprehensive and automated workflow for analyzing sequencing data.

## Installation
To run `SurveyFlow`, you need to set up a Conda environment with the required dependencies. We recommend using Mamba for faster installation.

1. **Install Mamba** (if not already installed):
   ```bash
   conda install -c conda-forge mamba
2. **Create the environment**
   ```bash
   mamba create -n surveyflow python=3.9 fastp jellyfish "r-base>=4.1,<4.2" r-ggplot2 r-argparse genomescope smudgeplot -c bioconda -c conda-forge -y
3. **Activate the environment**
   ```bash
   conda activate surveyflow
4. **Clone the repository**
   ```bash
   git clone https://github.com/lucien-ho/SurveyFlow.git

## Usage
Run the pipeline using the following command:
```bash
python surveyFlow.py --r1 <R1.fastq.gz> --r2 <R2.fastq.gz> --threads <threads> --kmer <kmer_size> --prefix <output_prefix> --ploidy <ploidy> --size <genome_size>
