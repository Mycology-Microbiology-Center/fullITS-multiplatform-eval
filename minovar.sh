#!/bin/bash

set -e  # Exit on error

## --- Default parameters (can be overridden by command-line arguments)

# Multithreading (e.g., for cutadapt, devider, minimap2, etc.)
THREADS=4

# Cutadapt parameters
CUTADAPT_ERROR_RATE=0.15
CUTADAPT_OVERLAP=15
CUTADAPT_MIN_LENGTH=400

# Clustering parameters
CLUSTER_ID_FIRST=0.8
CLUSTER_ID_SECOND=0.93
MIN_CLUSTER_SIZE=5

# Sampling parameters
MAX_READS_CONSENSUS=100
READS_FOR_POLISHING=200

# Variant calling parameters
VARIANT_QUALITY_THRESHOLD=10
MIN_COVERAGE=5
READ_LENGTH_FILTER=0.9  # 90% of consensus length

# Directory paths
READS_DIR="reads"
PRIMERS_FILE="primers.fas"
BC_GENES_FILE="bcGenes.txt"
DORADO_MODELS_DIR="doradoModels"

# Devider parameters
DEVIDER_PRESET="nanopore-r10"


## --- Usage information

usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

OPTIONS:
  Cutadapt parameters:
    --cutadapt-error-rate FLOAT       Error rate for primer matching (default: $CUTADAPT_ERROR_RATE)
    --cutadapt-overlap INT            Minimum overlap for primer detection (default: $CUTADAPT_OVERLAP)
    --cutadapt-min-length INT         Minimum length after trimming (default: $CUTADAPT_MIN_LENGTH)

  Clustering parameters:
    --cluster-id-first FLOAT          Identity threshold for first clustering (default: $CLUSTER_ID_FIRST)
    --cluster-id-second FLOAT         Identity threshold for second clustering (default: $CLUSTER_ID_SECOND)
    --min-cluster-size INT            Minimum reads per cluster (default: $MIN_CLUSTER_SIZE)

  Sampling parameters:
    --max-reads-consensus INT         Maximum reads for consensus generation (default: $MAX_READS_CONSENSUS)
    --reads-for-polishing INT         Reads to use for polishing (default: $READS_FOR_POLISHING)

  Variant calling parameters:
    --variant-quality-threshold INT   Minimum quality score for SNPs (default: $VARIANT_QUALITY_THRESHOLD)
    --min-coverage INT                Minimum coverage for variants (default: $MIN_COVERAGE)
    --read-length-filter FLOAT        Fraction of consensus length for filtering (default: $READ_LENGTH_FILTER)

  Directory and file paths:
    --reads-dir PATH                  Directory containing input reads (default: $READS_DIR)
    --primers-file PATH               Primers FASTA file (default: $PRIMERS_FILE)
    --bc-genes-file PATH              Barcode genes configuration file (default: $BC_GENES_FILE)
    --dorado-models-dir PATH          Dorado models directory (default: $DORADO_MODELS_DIR)

  Devider parameters:
    --devider-preset STRING           Devider preset (default: $DEVIDER_PRESET)

  Other:
    --threads INT                     Number of CPU threads for tools (default: $THREADS)
    -h, --help                        Display this help message and exit

EXAMPLES:
  # Run with default parameters
  ./$(basename "$0")

  # Run with customized parameters
  ./$(basename "$0") --cluster-id-first 0.85 --cluster-id-second 0.95 --min-coverage 10

EOF
    exit 1
}


## --- Parsing command-line arguments

while [[ $# -gt 0 ]]; do
    case $1 in
        # Cutadapt parameters
        --cutadapt-error-rate)
            CUTADAPT_ERROR_RATE="$2"
            shift 2 # move to the next argument
            ;;
        --cutadapt-overlap)
            CUTADAPT_OVERLAP="$2"
            shift 2
            ;;
        --cutadapt-min-length)
            CUTADAPT_MIN_LENGTH="$2"
            shift 2
            ;;
        --cluster-id-first)
            CLUSTER_ID_FIRST="$2"
            shift 2
            ;;
        --cluster-id-second)
            CLUSTER_ID_SECOND="$2"
            shift 2
            ;;
        --min-cluster-size)
            MIN_CLUSTER_SIZE="$2"
            shift 2
            ;;
        --max-reads-consensus)
            MAX_READS_CONSENSUS="$2"
            shift 2
            ;;
        --reads-for-polishing)
            READS_FOR_POLISHING="$2"
            shift 2
            ;;
        --variant-quality-threshold)
            VARIANT_QUALITY_THRESHOLD="$2"
            shift 2
            ;;
        --min-coverage)
            MIN_COVERAGE="$2"
            shift 2
            ;;
        --read-length-filter)
            READ_LENGTH_FILTER="$2"
            shift 2
            ;;
        --reads-dir)
            READS_DIR="$2"
            shift 2
            ;;
        --primers-file)
            PRIMERS_FILE="$2"
            shift 2
            ;;
        --bc-genes-file)
            BC_GENES_FILE="$2"
            shift 2
            ;;
        --dorado-models-dir)
            DORADO_MODELS_DIR="$2"
            shift 2
            ;;
        --devider-preset)
            DEVIDER_PRESET="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done



## --- Validation

echo "Validating dependencies..."

check_availability() {
  local tool="$1"
  if ! command -v "$tool" >/dev/null 2>&1; then
      echo -e "  \033[31mERROR: $tool is not installed\033[0m" >&2
      exit 1
  fi
}

check_availability "seqkit"
check_availability "cutadapt"
check_availability "vsearch"
check_availability "seqtk"
check_availability "samtools"
check_availability "minimap2"
check_availability "freebayes"
check_availability "dorado"
check_availability "abpoa"
check_availability "devider"
echo -e "  \033[32mAll dependencies are installed\033[0m" >&2


if [ ! -f $PRIMERS_FILE ]; then
  echo -e "  \033[31mERROR: $PRIMERS_FILE file not found\033[0m" >&2
  exit 1
fi
if [ ! -f $BC_GENES_FILE ]; then
  echo -e "  \033[31mERROR: $BC_GENES_FILE file not found\033[0m" >&2
  exit 1
fi
if [ ! -d "$DORADO_MODELS_DIR" ]; then
  echo -e "  \033[31mERROR: $DORADO_MODELS_DIR directory not found\033[0m" >&2
  echo -e "  \033[31mYou may specify a path to the directory with dorado models using --dorado-models-dir option\033[0m" >&2
  exit 1
fi
if [ ! -d "$READS_DIR" ]; then
  echo -e "  \033[31mERROR: $READS_DIR directory not found\033[0m" >&2
  exit 1
fi

# Check for at least one *.fastq.gz file in $READS_DIR and that each has a matching .bam
shopt -s nullglob
fastq_files=("$READS_DIR"/*.fastq.gz)
if [ ${#fastq_files[@]} -eq 0 ]; then
  echo -e "  \033[31mERROR: No .fastq.gz files found in $READS_DIR\033[0m" >&2
  exit 1
fi
echo -e "  \033[32mFound ${#fastq_files[@]} .fastq.gz file(s) in $READS_DIR\033[0m" >&2
missing_bam=()
for f in "${fastq_files[@]}"; do
  base=$(basename "$f" .fastq.gz)
  if [ ! -f "$READS_DIR/$base.bam" ]; then
    missing_bam+=("$base.bam")
  fi
done
if [ ${#missing_bam[@]} -gt 0 ]; then
  echo -e "  \033[31mERROR: Each .fastq.gz must have a matching .bam in $READS_DIR. Missing:\033[0m" >&2
  printf '    %s\n' "${missing_bam[@]}" >&2
  exit 1
fi
shopt -u nullglob



## --- Show specified parameters

echo "Configuration:"
echo "  CPU threads: $THREADS"
echo "  Cutadapt:"
echo "    Error rate: $CUTADAPT_ERROR_RATE"
echo "    Overlap:    $CUTADAPT_OVERLAP"
echo "    Min length: $CUTADAPT_MIN_LENGTH"
echo "  Clustering:"
echo "    First round ID:   $CLUSTER_ID_FIRST"
echo "    Second round ID:  $CLUSTER_ID_SECOND"
echo "    Min cluster size: $MIN_CLUSTER_SIZE"
echo "  Sampling:"
echo "    Max reads for consensus: $MAX_READS_CONSENSUS"
echo "    Reads for polishing:     $READS_FOR_POLISHING"
echo "  Variant calling:"
echo "    Quality threshold:  $VARIANT_QUALITY_THRESHOLD"
echo "    Min coverage:       $MIN_COVERAGE"
echo "    Read length filter: $READ_LENGTH_FILTER"
echo "  Paths:"
echo "    Reads directory: $READS_DIR"
echo "    Primers file:    $PRIMERS_FILE"
echo "    BC genes file:   $BC_GENES_FILE"
echo "    Dorado models:   $DORADO_MODELS_DIR"
echo ""


## --- Dir setup

mkdir -p vsearch
mkdir -p specimen_reads
mkdir -p cons_withPrimers
mkdir -p mixed
mkdir -p nonTargetAmp
mkdir -p clusters
mkdir -p clusters_small
mkdir -p bcCons
mkdir -p cons


## --- Main pipeline

seqkit seq -n $PRIMERS_FILE > gIDs.txt
find $READS_DIR/ -name "*.fastq.gz" | sed 's/.fastq.gz//' | sed "s|$READS_DIR/||" > bcSpecimens_all.txt

for b in $(find $READS_DIR/ -name "*.fastq.gz" | sed 's/.fastq.gz//' | sed "s|$READS_DIR/||")
do echo "Filtering reads of $b according to minimum and maximum length specified in the file $BC_GENES_FILE and sorting by average read quality in descending order"
  seqkit seq --quiet -m $(awk 'min>$2 || NR==1{min=$2} END{print min}' $BC_GENES_FILE) -M $(awk 'max<$3 || NR==1{max=$3} END{print max}' $BC_GENES_FILE) $READS_DIR/$b.fastq.gz | 
  seqkit replace -p "\s.+" > np_out.fq
  find vsearch/ -type f -delete
  cd vsearch
  echo "Classifying reads according to primers"
  cutadapt \
  -g file:../$PRIMERS_FILE \
  -e $CUTADAPT_ERROR_RATE \
  --revcomp \
  --action trim \
  --untrimmed-output ../nonTargetAmp/$b.fq \
  --overlap $CUTADAPT_OVERLAP \
  --minimum-length $CUTADAPT_MIN_LENGTH \
  --cores $THREADS \
  -o {name}.fq \
  ../np_out.fq
  find -name "*.fq" -type 'f' -empty -delete
if (($(find -name "*.fq" | grep -f ../gIDs.txt | wc -l) < 1))
  then
    echo "No target reads for $b"
  else
    for k in $(find -name "*.fq" | grep -o -f ../gIDs.txt)
    do seqkit replace -p "\s.+" $k.fq | seqkit fx2tab -q | sort -nrk4 | cut -f1,2,3 | seqkit tab2fx > np_out.fq
       mv np_out.fq $k.fq
       vsearch --cluster_smallmem $k.fq --usersort --clusters vcluster --id $CLUSTER_ID_FIRST --iddef 1 --threads $THREADS  # clustering reads # clustering reads without sorting reads in the input file
      if (($(find -name "vcluster*" 2> /dev/null | wc -l) < 1))
        then
        echo "no clusters for $b"
      else
        for j in $(find -name "vcluster*" | sed 's/.\///')
        do
          if (($(grep -c ">" $j) < $MIN_CLUSTER_SIZE)) # selecting clusters with at least MIN_CLUSTER_SIZE reads for consensus sequence creation
          then cat $j >> $k.fas
          else cp $j cluster.fas
          seqtk sample -s$(echo $RANDOM) cluster.fas $MAX_READS_CONSENSUS > inputc.fas # creating consensus sequence based on maximum of MAX_READS_CONSENSUS random reads of the cluster
          abpoa inputc.fas > temp.fas
          seqkit replace -p Consensus_sequence -r $b.$k.$j temp.fas >> $k.fas # collecting consensus sequences of the sample in the same file
          seqkit seq -n $j > $b.$k.$j.txt
          fi
        done
        find . -type f -name "vcluster*" -delete
        echo "Second round of clustering of $b $k clusters"
        seqkit seq --quiet -m $(grep $k -w ../$BC_GENES_FILE | cut -f2) $k.fas > filt.fas # filtering out sequences too short for the gene
        if (($(grep -c ">" filt.fas) < 1))
        then echo "No sequences of $k of sufficient length"
        else mv filt.fas $k.fas
             vsearch --cluster_fast $k.fas --clusters $k.vcluster --id $CLUSTER_ID_SECOND --iddef 1 --threads $THREADS  # second round of clustering
             find -name "$k.vcluster*" | sed 's/.\///' > IDs.txt
             find . -type f -name "readCount.txt" -delete
             for j in $(find -name "$k.vcluster*" | sed 's/.\///')
             do seqkit grep -rp vcluster $j | seqkit seq -n | wc -l >> readCount.txt # counting sequence names containing "vcluster" in clusters
             done
             paste IDs.txt readCount.txt | awk '$2 == 0' | awk '{print $1}' > temp.txt
             xargs -I{} rm -r "{}" < temp.txt # removing clusters not containing any sequence names containing "vcluster"
             if (($(find -name "$k.*" | wc -l) < 3))
             then echo "No clusters for $k"
             else
               for j in $(find -name "$k.vcluster*" | sed 's/.\///')
               do
                 if (($(grep -c ">" $j) < $MIN_CLUSTER_SIZE)) # selecting clusters with at least MIN_CLUSTER_SIZE reads for consensus sequence creation
                 then seqkit grep -rvp vcluster $j | seqkit seq -n > IDs.txt
                      seqkit grep -rp vcluster $j | seqkit seq -n | sed s/$/.txt/ > temp.txt
                      for v in $(cat temp.txt)
                      do cat $v | sort | uniq >> IDs.txt
                      done
                      samtools view -N IDs.txt ../$READS_DIR/$b.bam -O BAM -o ../specimen_reads/$b.$j.reads.bam
                      samtools fastq ../specimen_reads/$b.$j.reads.bam > ../specimen_reads/$b.$j.reads.fq
                      seqkit grep -f IDs.txt $k.fq > seq.fq
                      seqtk sample -s$(echo $RANDOM) seq.fq $MAX_READS_CONSENSUS > inputc.fas
                      abpoa inputc.fas > temp.fas
                      seqkit replace -p Consensus_sequence -r $b.$j temp.fas > ../specimen_reads/$b.$j.fas
                 else
                      seqtk sample -s$(echo $RANDOM) $j $MAX_READS_CONSENSUS > inputc.fas # creating consensus sequence based on maximum of MAX_READS_CONSENSUS random reads of the cluster
                      abpoa inputc.fas > temp.fas
                      seqkit replace -p Consensus_sequence -r $b.$j temp.fas > ../specimen_reads/$b.$j.fas
                      seqkit grep -rvp vcluster $j | seqkit seq -n > IDs.txt
                      seqkit grep -rp vcluster $j | seqkit seq -n | sed s/$/.txt/ > temp.txt
                      for v in $(cat temp.txt)
                      do cat $v | sort | uniq >> IDs.txt
                      done
                      samtools view -N IDs.txt ../$READS_DIR/$b.bam -O BAM -o ../specimen_reads/$b.$j.reads.bam
                      samtools fastq ../specimen_reads/$b.$j.reads.bam > ../specimen_reads/$b.$j.reads.fq
                 fi
                 cat ../specimen_reads/$b.$k.*.fas | seqkit rmdup -s -D ../clusters/$b.$k.IDs.txt > ../clusters/$b.$k.classified.fas # removing duplicate identical sequences and getting IDs of duplicate sequences
               done
             fi
        fi
      fi
    done
fi
cd ..
done
#
  echo "Removing duplicate identical sequences and merging corresponding *.reads.bam and *.reads.fq files"
  if (($(find clusters/ -name "*.IDs.txt" | wc -l) < 1))
  then echo "No duplicate sequences found"
       else for j in $(find clusters/ -name "*.IDs.txt" | sed 's/.IDs.txt//' | sed 's/clusters\///')
            do while read line
               do echo $line | sed 's/ /\t/' | cut -f2 | sed 's/, /\n/g' | awk '{print "mv specimen_reads/"$1".reads.fq cons/"$1".reads.fq"}' > mv.sh
                  echo $line | sed 's/ /\t/' | cut -f2 | sed 's/, /\n/g' | awk '{print "rm specimen_reads/"$1".reads.bam"}' >> mv.sh
                  echo $line | sed 's/ /\t/' | cut -f2 | sed 's/, /\n/g' | tail -n+2 | awk '{print "rm specimen_reads/"$1".fas"}' >> mv.sh
                  chmod +x mv.sh
                  ./mv.sh
                  toMerge=$(echo $line | sed 's/ /\t/' | cut -f2 | sed 's/, /\n/g' | head -n1)
                  cat cons/*.fq > specimen_reads/$toMerge.reads.fq
                  seqkit seq -n specimen_reads/$toMerge.reads.fq > IDs.txt
                  samtools view -N IDs.txt $READS_DIR/$(echo $j | grep -wo -f bcSpecimens_all.txt).bam -O BAM -o specimen_reads/$toMerge.reads.bam
                  rm cons/*.fq
                  rm -f ./mv.sh
               done < clusters/$j.IDs.txt
            done
  fi
#
  echo "Moving small clusters to the folder clusters_small"
  rm -f ./mv.sh
  cd specimen_reads
  for i in $(find -name "*.fq" | grep -f ../bcSpecimens_all.txt -o | sort | uniq)
  do seqkit stats $i.*.fq -T | awk '$4 < 5' | awk '{print $1}' | awk '{print "mv specimen_reads/"$1" clusters_small/"$1}' >> ../mv.sh
     seqkit stats $i.*.fq -T | awk '$4 < 5' | awk '{print $1}' | sed s/reads.fq/fas/ | awk '{print "mv specimen_reads/"$1" clusters_small/"$1}' >> ../mv.sh
     seqkit stats $i.*.fq -T | awk '$4 < 5' | awk '{print $1}' | sed s/reads.fq/reads.bam/ | awk '{print "mv specimen_reads/"$1" clusters_small/"$1}' >> ../mv.sh
  done
  cd ..
  chmod +x mv.sh
  ./mv.sh
#
echo "Polishing and variant calling"
  for j in $(find specimen_reads/ -name "*.fas" | sed 's/.fas//' | sed 's/specimen_reads\///')
  do echo "Consensus sequence polishing of $j"
     seqtk sample -s$(echo $RANDOM) specimen_reads/$j.reads.fq $READS_FOR_POLISHING | seqkit seq -n > IDs.txt # taking READS_FOR_POLISHING random reads for consensus polishing
     samtools view -N IDs.txt specimen_reads/$j.reads.bam -O BAM -o seq.bam
     find . -type f -name "*.bai" -delete
     find . -type f -name "*.fai" -delete
     cp specimen_reads/$j.fas cons.fa # dorado polish accepts only *.fasta or *.fa extensions, not *.fas
     samtools faidx cons.fa
     dorado aligner cons.fa seq.bam | samtools sort > alignment.bam
     samtools index alignment.bam
     dorado polish alignment.bam cons.fa --ignore-read-groups --models-directory $DORADO_MODELS_DIR > consmed.fas # dorado polish alignment.bam cons.fa --ignore-read-groups > consmed.fas
     seqkit replace -p $(seqkit seq -n consmed.fas) -r $j consmed.fas > specimen_reads/$j.fas
     echo "Variant calling and consensus sequence creation and polishing of $j"
     if (($(samtools view -c seq.bam) < 5))
     then echo "Fewer than 5 reads for $j"
     else rlen=$(seqkit stats -T specimen_reads/$j.fas | cut -f5 | tail -n+2 | awk -v filter="$READ_LENGTH_FILTER" '{print int($1*filter)}') # getting READ_LENGTH_FILTER of consensus sequence length to filter out too short reads after mapping
         minimap2 -a --sam-hit-only -x map-ont --secondary=no -t $THREADS specimen_reads/$j.fas specimen_reads/$j.reads.fq | samtools sort | samtools view -e "rlen>=$rlen" -O BAM > alignment.bam # $rlen works only with double quotes in samtools view -e
          samtools index alignment.bam
          samtools faidx specimen_reads/$j.fas
          freebayes -i --haplotype-length -1 -f specimen_reads/$j.fas alignment.bam > var.vcf # variant calling
          grep '#' var.vcf > varf.vcf # creating a separate file for detected variable positions or SNPs having high quality, first adding only the header containing comments of the original file
          grep -v '#' var.vcf | awk -v threshold="$VARIANT_QUALITY_THRESHOLD" '$6>threshold' >> varf.vcf # adding only detected SNPs with high quality (score higher than VARIANT_QUALITY_THRESHOLD)
        if (($(grep -v '#' varf.vcf | wc -l) < 1)) # if no high quality SNPs were detected, do not process the consensus sequence further
        then echo "No SNPs for $j detected"
             cp alignment.bam specimen_reads/$j.bam
             samtools index specimen_reads/$j.bam
             else devider -b alignment.bam -v varf.vcf -r specimen_reads/$j.fas -o devider_output -t $THREADS --preset $DEVIDER_PRESET --min-cov $MIN_COVERAGE -O # keeping only variants supported by at least MIN_COVERAGE reads
             if [ ! -f devider_output/majority_vote_haplotypes.fasta ]
             then echo "No variants with more than 4 reads detected"
                  mv specimen_reads/$j.reads.fq mixed/
                  mv specimen_reads/$j.fas mixed/
             else seqkit seq -n devider_output/majority_vote_haplotypes.fasta | cut -d ',' -f3 > temp.txt # haplotype IDs
                  samtools view -f 16 alignment.bam | cut -f1 | sort | uniq > rev.txt # reverse complemented read IDs
                  samtools view alignment.bam | awk '$2 == 0' | cut -f1 | sort | uniq > for.txt  # forward read IDs
                  seqkit grep -f for.txt specimen_reads/$j.reads.fq | seqkit fq2fa > for.fas
                  seqkit grep -f rev.txt specimen_reads/$j.reads.fq | seqkit seq -p -r -v -t DNA | seqkit fq2fa >> for.fas # all reads in the same orientation
               for v in $(cat temp.txt)
               do grep $v devider_output/ids.txt | cut -f 4- | sed 's/\t/\n/g' > IDs.txt
                  seqkit grep -f IDs.txt specimen_reads/$j.reads.fq > specimen_reads/$j.$(echo $v | sed 's/Haplotype://').reads.fq # copying fastq reads of the variant to folder specimen_reads
                  seqkit grep -f IDs.txt for.fas > inputvar.fas
                  seqtk sample -s$(echo $RANDOM) inputvar.fas $MAX_READS_CONSENSUS > inputc.fas # MAX_READS_CONSENSUS random reads for consensus sequence computation
                  abpoa inputc.fas > constemp.fas # computing consensus sequence or a variant
                  seqkit replace -p $(seqkit seq -n constemp.fas) -r $j.$(echo $v | sed 's/Haplotype://') constemp.fas > cons.fa # renaming consensus sequence of the variant # dorado polish accepts only *.fasta or *.fa extensions, not *.fas
                  seqtk sample -s$(echo $RANDOM) specimen_reads/$j.$(echo $v | sed 's/Haplotype://').reads.fq $READS_FOR_POLISHING | seqkit seq -n > IDs.txt
                  samtools view -N IDs.txt specimen_reads/$j.reads.bam -O BAM -o seq.bam
                  find . -type f -name "*.bai" -delete
                  find . -type f -name "*.fai" -delete
                  samtools faidx cons.fa
                  dorado aligner cons.fa seq.bam | samtools sort > alignment.bam
                  samtools index alignment.bam
                  dorado polish alignment.bam cons.fa --ignore-read-groups --models-directory $DORADO_MODELS_DIR > cons_withPrimers/$j.$(echo $v | sed 's/Haplotype://').fas
               done
                 mv specimen_reads/$j.reads.fq mixed/
                 mv specimen_reads/$j.fas mixed/
             fi
        fi
     fi
  done
#
  echo "Collecting same amplicon single variant sequences into one file in the folder bcCons"
  cd specimen_reads
  find . -type f -name "*.fai" -delete
  if (($(find -name "*.fas" | grep -f ../gIDs.txt -o | wc -l) < 1))
  then echo "No single variant sequences"
  else 
    for j in $(find -name "*.fas" | sed 's/.\///' | sed 's/.vcluster.*$//' | sort | uniq)
    do cat $j*.fas > ../bcCons/$j.new.fas
    done
  fi
  cd ..
#
  echo "Removing primers from multi variant sequences, creating mapping files for each variant, and adding the variants to corresponding amplicon files in the folder bcCons"
  cd cons_withPrimers
  if (($(find -name "*.fas" | grep -f ../gIDs.txt -o | wc -l) < 1))
  then echo "No sequences with more than one variant"
  else 
     for j in $(find -name "*.fas" | sed 's/.\///' | sed 's/.vcluster.*$//' | sort | uniq)
     do echo "Cutting primers"
     cat $j*.fas > withPrimers.fasta
     find . -type f -name "primersNotCut.fasta" -delete
     cutadapt \
       -g file:../$PRIMERS_FILE \
       -e $CUTADAPT_ERROR_RATE \
       --action trim \
       --untrimmed-output primersNotCut.fasta \
       --overlap $CUTADAPT_OVERLAP \
       --minimum-length $CUTADAPT_MIN_LENGTH \
       --cores $THREADS \
       -o noPrimers.fasta \
       withPrimers.fasta
       find -name "*.fasta" -type 'f' -empty -delete
       cat noPrimers.fasta >> ../bcCons/$j.new.fas
       if [ -f primersNotCut.fasta ]; then
         cat primersNotCut.fasta >> ../bcCons/primersNotCut.fas
       fi
     done
  fi
  cd ../bcCons
  for i in $(find -name "*.new.fas" | sed 's/.\///' | sed 's/.new.fas//')
  do echo "Removing duplicate consensus sequences of $i"
     seqkit rmdup -s -D $i.IDs.txt $i.new.fas > $i.deduplicated.fas
  done
#
  if (($(find -name "*.IDs.txt" | wc -l) < 1))
  then echo "No duplicate sequences found"
    else for i in $(find -name "*.IDs.txt" | sed 's/.\///' | sed 's/.IDs.txt//')
         do echo "Merging read files of identical consensus sequences of $i"
         while read line
           do echo $line | sed 's/ /\t/' | cut -f2 | sed 's/, /\n/g' | awk '{print "mv ../specimen_reads/"$1".reads.fq ../cons/"$1".reads.fq"}' > mv.sh
              chmod +x mv.sh
              ./mv.sh
              toMerge=$(echo $line | sed 's/ /\t/' | cut -f2 | sed 's/, /\n/g' | head -n1)
              cat ../cons/*.fq > ../specimen_reads/$toMerge.reads.fq
              find ../cons -type f -name "*.fq" -delete
              rm -f ./mv.sh
           done < $i.IDs.txt
         done
  fi
#
  cd ..
  cd specimen_reads
  for i in $(find -name "*reads.fq" | grep -f ../bcSpecimens_all.txt -o | sort | uniq)
  do seqkit stats $i.*reads.fq -T | tail -n +2 | awk '{print $1, $4}' | sed 's/.reads.fq//' | sed 's/.vcluster/ vcluster/' | sort -k1,1 -k3,3nr | awk '{print $1"."$2, $0}' >> ../bcCons/readCountID.txt
  done
  cd ..
#
