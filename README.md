# fullITS-multiplatform-eval
Benchmarking Nanopore, PacBio, and Illumina 2×500 for eukaryotic full-length ITS amplicon sequencing

## Minovar

The script produces polished consensus variant sequences from nanopore reads of amplicons.

The script takes as input
1) fastq reads, one or more (e.g. demultiplexed) files
2) a fasta file containing primer sequences
3) a text file specifying expected amplicon min and max lengths.

The primer fasta file hast to be in the following form suitable for cutadapt:

```
>ITS9mun_ITS4ngsUni
GTACACACCGCCCGTCG...GCATATHANTAAGSGSAGGCG
>SymF3v2_SymR3
TTTTCWACWAAYCAYAAARAYATYGG...TGATTYTTTGGNCAYCCWGAAGTTTA
```

Forward primer has to follow the reverse primer on the same line separated by '...' and reverse primer has to be reverse complemented (3'-5' direction).

The primer file can contain as many primers for different amplicons as needed, except that the amplicon regions should not overlap with each other. The overlapping amplicons should be processed separately.

The text file for amplicon lengths has to contain at least three tab separated columns.

```
ITS9mun_ITS4ngsUni	450	3001
SymF3v2_SymR3	600	900
```

The first column contains primer name that matches the name in the fasta file. Second column is min and the third column max expected amplicon length. This is mainly to reduce noise in the output. Considering the indel errors in nanopore reads, it is worth setting the min length lower and max higher by 5-10%.


Each fastq input file is processed separately to get initial consensus variant sequences for each amplicon. Cutadapt (>3.4) classifies reads to amplicons according to the primer sequences. Then there are two clustering steps with vsearch. First, the classified reads are sorted by average read quality in decreasing order and clustered by vsearch at 80% similarity (`vsearch --cluster_smallmem --usersort --id 0.8 --iddef`). Despite this low similarity, vsearch tends to separate much more similar variants (>95% similarity) into separate clusters (perhaps because it is not optimized for error profile of nanopore reads). The initial clustering similarity threshold is set low to avoid oversplitting and creating too many small clusters (many of which may belong to identical variants) that are not good enough to create initial consensus sequences. Minimum of 5 and maximum of 100 random reads (too high coverage can be detrimental to consensus accuracy) of a cluster are used to create a consensus sequence with abpoa. The consensus sequences are clustered in the second vsearch step at 93% similarity (`vsearch --cluster_fast --id 0.93 --iddef`) and the initial clusters are accordingly merged. This is followed by new consensus sequence creation with abpoa. Identical consensus sequences are then detected with seqkit (`seqkit rmdup -s -D`) and the corresponding reads are merged into one file.

The consensus sequences are polished with dorado polish, assuming the relevant polishing model has been downloaded to directory 'doradoModels'. The polished consensus sequences are then used for variant calling (minimap2 and freebayes). If high quality SNPs (score higher than 10) are detected, these will be used together with the read mapping file (`*.bam`) and the reference consensus sequence file to phase reads with devider. devider is able to detect even single SNP variants of very low frequency (default minimum is 0.25%) from high sequencing depth (>10000X) without the need to specify the expected number of variants. The limitation is that if there are indel-only variants (i.e. no substitutions), these will not be detected, but the indels will be taken into account in the following steps if in addition there is any substitution associated with the indel(s). The original reads corresponding to the detected variants are used to create consensus sequences separately for each variant and polished as before (`abpoa` and `dorado polish`).


The final deduplicated polished consensus sequences are found in the bcCons folder and the corresponding reads in the specimen_reads folder. Read counts for each variant consensus sequence are given in the space delimited file readCountID.txt in bcCons folder.

### Installation

To install dependencies, you can use [conda](https://docs.anaconda.com/miniconda/) or [mamba](https://mamba.readthedocs.io/en/latest/):

```bash
## Create conda environment with required dependencies
mamba create -n Minovar \
  -c conda-forge -c bioconda \
  seqkit==2.12.0 \
  seqtk==1.5 \
  cutadapt==5.2 \
  vsearch==2.30.4 \
  samtools==1.23 \
  minimap2==2.30 \
  freebayes==1.3.10 \
  abpoa==1.5.6 \
  devider==0.0.1 \
  python==3.12

## Add Dorado to the conda environment
conda activate Minovar
mkdir -p "$CONDA_PREFIX/opt"
cd "$CONDA_PREFIX/opt"
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.3.1-linux-x64.tar.gz
tar -xzf dorado-1.3.1-linux-x64.tar.gz
rm dorado-1.3.1-linux-x64.tar.gz
ln -sf "$CONDA_PREFIX/opt/dorado-1.3.1-linux-x64/bin/dorado" "$CONDA_PREFIX/bin/dorado"
dorado --help

## Download dorado model
# cd ~/Nanopore/    # select directory where you want to download the model

doradoModel=dna_r10.4.1_e8.2_400bps_sup@v5.2.0_polish_rl_mv
mkdir -p doradoModels
dorado download --model $doradoModel --models-directory doradoModels
```

### Usage

1. Pull Minovar code to your local machine

```bash
wget https://raw.githubusercontent.com/Mycology-Microbiology-Center/fullITS-multiplatform-eval/refs/heads/main/minovar.sh
chmod +x minovar.sh
```

2. Prepare files with primer and amplicon length information.
For examples, see:
```bash
wget https://raw.githubusercontent.com/Mycology-Microbiology-Center/fullITS-multiplatform-eval/refs/heads/main/bcGenes.txt
wget https://raw.githubusercontent.com/Mycology-Microbiology-Center/fullITS-multiplatform-eval/refs/heads/main/primers.fas
```

3. Execute the workflow

To run the pipleline with the default parameters, use the following command:

```bash
conda activate Minovar
./minovar.sh
```

> [!NOTE]
> ```
> In addition to `bcGenes.txt` and `primers.fas`, the pipeline expects input data
> (`fastq.gz` files that can have move tables for Dorado polishing) in the `reads` subdirectory,
> and dorado models in the `doradoModels` subdirectory of the current working directory.
> Input `fastq.gz` file sequence identifier line has to include basecalling model for dorado polishing
> (which is for example added with --emit-fastq flag with dorado basecaller).
> 
> If basecaller model is missing, it can be added in the following form.
> 
> Original sequence identifier:
> @sequenceID
> New sequence identifier:
> @sequenceID RG:Z:-_dna_r10.4.1_e8.2_400bps_sup@v4.3.0
> 
> Sequence identifier and basecalling model can be separated by space or tab.
> ```

To use customized parameters, run the script with the desired parameters. For example:

```bash
conda activate Minovar
./minovar.sh \
  --cutadapt-error-rate 0.15 \
  --cutadapt-overlap 15 \
  --cutadapt-min-length 400 \
  --cluster-id-first 0.8 \
  --cluster-id-second 0.93 \
  --min-cluster-size 5 \
  --max-reads-consensus 100 \
  --reads-for-polishing 200 \
  --variant-quality-threshold 10 \
  --min-coverage 5 \
  --read-length-filter 0.9 \
  --reads-dir reads \
  --dorado-models-dir doradoModels
```

To see all available parameters, run the script with the `--help` option:

```bash
conda activate Minovar
./minovar.sh --help
```

