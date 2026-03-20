# Benchmarking full-length ITS metabarcoding across Illumina 2x500, PacBio, and Oxford Nanopore sequencing using mock and soil communities

## Contents

- [Minovar](#minovar)
- [PRONAME](#proname)
- [NextITS](#nextits)


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

## PRONAME

The analysis of ONT sequencing data was also performed using the [PRONAME](https://github.com/benn888/PRONAME) pipeline. It consists of four main scripts ((i) `proname_import`, (ii) `proname_filter`, (iii) `proname_refine` and (iv) `proname_taxonomy`), but only the first three were used here, as the taxonomic assignment step was performed outside PRONAME in order to standardize the procedure across datasets from all sequencing platforms.

First, demultiplexed FASTQ sequencing files, placed in the `RawData` folder, were imported into PRONAME v2.2.0:

```bash
proname_import \
  --inputpath RawData \
  --duplex no \
  --trimadapters no \
  --trimprimers no \
  --threads 48 \
  --plotformat html
```

Low-quality reads were then discarded by retaining only reads ranging from 450 to 3000 bp with a Phred score of at least 15:

```bash
proname_filter \
  --datatype simplex \
  --filtminlen 450 \
  --filtmaxlen 3000 \
  --filtminqual 15 \
  --inputpath RawData \
  --threads 48 \
  --plotformat html
```

The core of the analysis was then performed using `proname_refine` to cluster reads according to a 98% identity threshold and remove singletons, followed by an error-correction step using `dorado polish`. The `proname_refine` script was slightly modified to deactivate the chimera detection, as this step was performed at a later stage for all datasets.

```bash
proname_refine \
  --clusterid 0.98 \
  --clusterthreads 48 \
  --inputpath RawData \
  --minreadspercluster 2 \
  --polisher dorado \
  --polisherthreads 48 \
  --chimeramethod denovo \
  --qiime2import no
```

Full-length ITS regions were then extracted using ITSx and an updated hidden Markov model (HMM) profile database:

```bash
docker run --rm -v "$PWD":/data vmikk/nextits:1.1.0 \
  ITSx -i /data/rep_seqs.fasta -o /data/ITSx_extracted \
      --complement T \
      --save_regions all \
      --positions T \
      --not_found T \
      -E 0.1 \
      -t all \
      --partial 0 \
      --cpu 48 \
      --preserve T
```

The abundance table was then filtered accordingly, and ITS-extracted sequences were dereplicated:

```bash
# Filtering out from abundance table (rep_table.tsv) the sequences that were discarded during ITS extraction
grep "^>" ITSx_extracted.full.fasta | sed 's/^>//' > ids_kept.txt

awk -F'\t' 'FNR==NR { gsub(/\r/,"",$1); keep[$1]; next } FNR==1 || ($1 in keep)' \
  ids_kept.txt rep_table.tsv > ITSx_extracted_rep_table.tsv

rm ids_kept.txt

mv ITSx_extracted.full.fasta ITSx_extracted_rep_seqs.fasta


# Dereplicating ITS extracted sequences
vsearch \
  --derep_fulllength ITSx_extracted_rep_seqs.fasta \
  --strand both \
  --fasta_width 0 \
  --uc ITSx_extracted_rep_seqs.derep.uc \
  --output ITSx_extracted_rep_seqs.derep.fasta \
  --threads 48


# Dereplicating the abundance table accordingly
awk -F'\t' '
  BEGIN{ OFS="\t" }
  FNR==NR {
    gsub(/\r/,"")
    # S: centroid -> itself ; H: query -> centroid
    if ($1=="S") { q=$(NF-1); map[q]=q }
    else if ($1=="H") { q=$(NF-1); t=$NF; map[q]=t }
    next
  }
  FNR==1 { header=$0; ncol=NF; next }
  {
    id=$1
    c = (id in map) ? map[id] : id
    if (!(c in seen)) { order[++m]=c; seen[c]=1 }
    for (i=2; i<=ncol; i++) sum[c,i] += $i+0
  }
  END {
    print header
    for (k=1; k<=m; k++) {
      c = order[k]
      printf "%s", c
      for (i=2; i<=ncol; i++) printf OFS "%d", sum[c,i]
      printf "\n"
    }
  }
' ITSx_extracted_rep_seqs.derep.uc ITSx_extracted_rep_table.tsv > ITSx_extracted_rep_table.derep.tsv
```

## NextITS

[NextITS](https://github.com/vmikk/NextITS) v1.1.0 ([Zenodo](https://doi.org/10.5281/zenodo.15074882)) was run with Nextflow. It was applied to PacBio HiFi reads and to Illumina reads after demultiplexing, paired-end merging, and primer trimming, using the same NextITS settings as for PacBio. Step-1 performs QC (including LIMA demultiplexing for PacBio). Step-2 pools Step-1 outputs and clusters sequences. Full pipeline options, container profiles, and HPC notes are available in the [NextITS documentation](https://next-ITS.github.io/).

### Step-1

Standard settings (LIMA minimum score 85, full ITS, reference-based chimera checking with a EUKARYOME-compatible USEARCH UDB):

```bash
INPUT_FASTQ="/path/to/run/*.fastq.gz"
BARCODES_FASTA="/path/to/run/*_barcodes.fasta"
CHIMERA_DB="/path/to/Eukaryome_1.9.3_241222_FullITS_100-800.udb"
OUTDIR_STEP1="$(pwd)/Step1_Results/PACBvs"
WORK_DIR_STEP1="$(pwd)/Step1_work/PACBvs"
mkdir -p "$OUTDIR_STEP1" "$WORK_DIR_STEP1"

export NXF_OPTS="-Xms500M -Xmx3G"

nextflow run vmikk/NextITS -r main \
  -profile singularity \
  -resume \
  --step     "Step1" \
  --input    "$INPUT_FASTQ" \
  --barcodes "$BARCODES_FASTA" \
  --primer_forward   GTACACACCGCCCGTCG \
  --primer_reverse   CCTSCSCTTANTDATATGC \
  --lima_barcodetype "dual_symmetric" \
  --lima_minscore    85 \
  --its_region "full" \
  --chimera_db "$CHIMERA_DB" \
  --outdir     "$OUTDIR_STEP1" \
  -work-dir    "$WORK_DIR_STEP1"
```

On HPC clusters, switch to `-profile singularity,hpc`, set `NXF_SINGULARITY_CACHEDIR` if needed, and pass `-qs "$SLURM_NTASKS"` as described in the [HPC](https://next-its.github.io/HPC/) instructions.


Minimum-features Step-1 (chimera removal, tag-jump filtering, and homopolymer correction disabled):

```bash
LABEL="PACBvs"
INPUT_FASTQ="/path/to/${LABEL}/*.fastq.gz"
BARCODES_FASTA="/path/to/${LABEL}/*_barcodes.fasta"
OUTDIR_STEP1="$(pwd)/Step1_Results/${LABEL}"
WORK_DIR_STEP1="$(pwd)/Step1_work/${LABEL}"
mkdir -p "$OUTDIR_STEP1" "$WORK_DIR_STEP1"

nextflow run vmikk/NextITS -r main \
  -profile singularity \
  -resume \
  --step     "Step1" \
  --input    "$INPUT_FASTQ" \
  --barcodes "$BARCODES_FASTA" \
  --primer_forward   GTACACACCGCCCGTCG \
  --primer_reverse   CCTSCSCTTANTDATATGC \
  --lima_barcodetype "dual_symmetric" \
  --lima_minscore    85 \
  --its_region "full" \
  --chimera_methods "none" \
  --hp false \
  --tj false \
  --outdir  "$OUTDIR_STEP1" \
  -work-dir "$WORK_DIR_STEP1"
```

### Step-2 (clustering)

UNOISE pre-clustering and VSEARCH clustering at 98% identity (`--data_path` is the parent directory that contains the Step-1 output folders):

```bash
DATA_PATH="$(pwd)/Step1_Results"
OUTDIR_STEP2="$(pwd)/Step2_Results"
WORK_DIR_STEP2="$(pwd)/Step2_work"
mkdir -p "$OUTDIR_STEP2" "$WORK_DIR_STEP2"

export NXF_OPTS="-Xms500M -Xmx2G"

nextflow run vmikk/NextITS -r main \
  -profile singularity \
  -resume \
  --step "Step2" \
  --ampliconlen_min 100 \
  --preclustering "unoise" \
  --unoise_alpha 6.0 \
  --unoise_minsize 1 \
  --clustering "vsearch" \
  --otu_id 0.98 \
  --merge_replicates false \
  --max_MEEP 0.6 \
  --max_ChimeraScore 0.6 \
  --lulu false \
  --data_path "$DATA_PATH" \
  --outdir    "$OUTDIR_STEP2" \
  -work-dir   "$WORK_DIR_STEP2"
```
