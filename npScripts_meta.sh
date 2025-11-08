seqkit seq -n primers.fas > gIDs.txt
ls reads/*.fastq.gz | sed 's/.fastq.gz//' | sed 's/reads\///' > bcSpecimens_all.txt
#doradoModel=dna_r10.4.1_e8.2_400bps_sup@v5.2.0_polish_rl_mv
#mkdir doradoModels
#dorado download --model $doradoModel --models-directory doradoModels
for b in $(ls reads/*.fastq.gz | sed 's/.fastq.gz//' | sed 's/reads\///')
do echo "Filtering reads of $b according to minimum and maximum length specified in the file bcGenes.txt and sorting by average read quality in descending order"
  seqkit seq -m $(awk 'min>$2 || NR==1{min=$2} END{print min}' bcGenes.txt) -M $(awk 'max<$3 || NR==1{max=$3} END{print max}' bcGenes.txt) reads/$b.fastq.gz | 
  seqkit replace -p "\s.+" > np_out.fq
  cd vsearch
  rm *
  echo "Classifying reads according to primers"
  cutadapt \
  -g file:../primers.fas \
  -e 0.15 \
  --revcomp \
  --action trim \
  --untrimmed-output ../nonTargetAmp/$b.fq \
  --overlap 15 \
  --minimum-length 400 \
  --cores 3 \
  -o {name}.fq \
  ../np_out.fq
  find -name "*.fq" -type 'f' -empty -delete
if (($(ls *.fq | grep -f ../gIDs.txt | wc -l) < 1))
  then
    echo "No target reads for $b"
  else
    for k in $(ls *.fq | grep -o -f ../gIDs.txt)
    do seqkit replace -p "\s.+" $k.fq | seqkit fx2tab -q | sort -nrk4 | cut -f1,2,3 | seqkit tab2fx > np_out.fq
       mv np_out.fq $k.fq
       vsearch --cluster_smallmem $k.fq --usersort --clusters vcluster --id 0.8 --iddef 1 # clustering reads # clustering reads without sorting reads in the input file
      if (($(ls vcluster* 2> /dev/null | wc -l) < 1))
        then
        echo "no clusters for $b"
      else
        for j in $(ls vcluster*)
        do
          if (($(grep -c ">" $j) < 5)) # selecting clusters with at least 5 reads for consensus sequence creation
          then cat $j >> $k.fas
          else cp $j cluster.fas
          seqtk sample -s$(echo $RANDOM) cluster.fas 100 > inputc.fas # creating consensus sequence based on maximum of 100 random reads of the cluster
          abpoa inputc.fas > temp.fas
          seqkit replace -p Consensus_sequence -r $b.$k.$j temp.fas >> $k.fas # collecting consensus sequences of the sample in the same file
          seqkit seq -n $j > $b.$k.$j.txt
          fi
        done
        rm vcluster*
        echo "Second round of clustering of $b $k clusters"
        seqkit seq -m $(grep $k -w ../bcGenes.txt | cut -f2) $k.fas > filt.fas # filtering out sequences too short for the gene
        if (($(grep -c ">" filt.fas) < 1))
        then echo "No sequences of $k of sufficient length"
        else mv filt.fas $k.fas
             vsearch --cluster_fast $k.fas --clusters $k.vcluster --id 0.93 --iddef 1 # second round of clustering
             ls $k.vcluster* > IDs.txt
             rm readCount.txt
             for j in $(ls $k.vcluster*)
             do seqkit grep -rp vcluster $j | seqkit seq -n | wc -l >> readCount.txt # counting sequence names containing "vcluster" in clusters
             done
             paste IDs.txt readCount.txt | awk '$2 == 0' | awk '{print $1}' > temp.txt
             xargs -I{} rm -r "{}" < temp.txt # removing clusters not containing any sequence names containing "vcluster"
             if (($(ls $k.* | wc -l) < 3))
             then echo "No clusters for $k"
             else
               for j in $(ls $k.vcluster*)
               do
                 if (($(grep -c ">" $j) < 5)) # selecting clusters with at least 5 reads for consensus sequence creation
                 then seqkit grep -rvp vcluster $j | seqkit seq -n > IDs.txt
                      seqkit grep -rp vcluster $j | seqkit seq -n | sed s/$/.txt/ > temp.txt
                      for v in $(cat temp.txt)
                      do cat $v | sort | uniq >> IDs.txt
                      done
                      samtools view -N IDs.txt ../reads/$b.bam -O BAM -o ../specimen_reads/$b.$j.reads.bam
                      samtools fastq ../specimen_reads/$b.$j.reads.bam > ../specimen_reads/$b.$j.reads.fq
                      seqkit grep -f IDs.txt $k.fq > seq.fq
                      seqtk sample -s$(echo $RANDOM) seq.fq 100 > inputc.fas
                      abpoa inputc.fas > temp.fas
                      seqkit replace -p Consensus_sequence -r $b.$j temp.fas > ../specimen_reads/$b.$j.fas
                 else
                      seqtk sample -s$(echo $RANDOM) $j 100 > inputc.fas # creating consensus sequence based on maximum of 100 random reads of the cluster
                      abpoa inputc.fas > temp.fas
                      seqkit replace -p Consensus_sequence -r $b.$j temp.fas > ../specimen_reads/$b.$j.fas
                      seqkit grep -rvp vcluster $j | seqkit seq -n > IDs.txt
                      seqkit grep -rp vcluster $j | seqkit seq -n | sed s/$/.txt/ > temp.txt
                      for v in $(cat temp.txt)
                      do cat $v | sort | uniq >> IDs.txt
                      done
                      samtools view -N IDs.txt ../reads/$b.bam -O BAM -o ../specimen_reads/$b.$j.reads.bam
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
  if (($(ls clusters/*.IDs.txt | wc -l) < 1))
  then echo "No duplicate sequences found"
       else for j in $(ls clusters/*.IDs.txt | sed 's/.IDs.txt//' | sed 's/clusters\///')
            do while read line
               do echo $line | sed 's/ /\t/' | cut -f2 | sed 's/, /\n/g' | awk '{print "mv specimen_reads/"$1".reads.fq cons/"$1".reads.fq"}' > mv.sh
                  echo $line | sed 's/ /\t/' | cut -f2 | sed 's/, /\n/g' | awk '{print "rm specimen_reads/"$1".reads.bam"}' >> mv.sh
                  echo $line | sed 's/ /\t/' | cut -f2 | sed 's/, /\n/g' | tail -n+2 | awk '{print "rm specimen_reads/"$1".fas"}' >> mv.sh
                  chmod +x mv.sh
                  ./mv.sh
                  toMerge=$(echo $line | sed 's/ /\t/' | cut -f2 | sed 's/, /\n/g' | head -n1)
                  cat cons/*.fq > specimen_reads/$toMerge.reads.fq
                  seqkit seq -n specimen_reads/$toMerge.reads.fq > IDs.txt
                  samtools view -N IDs.txt reads/$(echo $j | grep -wo -f bcSpecimens_all.txt).bam -O BAM -o specimen_reads/$toMerge.reads.bam
                  rm cons/*.fq
               done < clusters/$j.IDs.txt
            done
  fi
#
  echo "Moving small clusters to the folder clusters_small"
  cd specimen_reads
  seqkit stats *.fq -T | awk '$4 < 5' | awk '{print $1}' | awk '{print "mv specimen_reads/"$1" clusters_small/"$1}' > ../mv.sh
  seqkit stats *.fq -T | awk '$4 < 5' | awk '{print $1}' | sed s/reads.fq/fas/ | awk '{print "mv specimen_reads/"$1" clusters_small/"$1}' >> ../mv.sh
  seqkit stats *.fq -T | awk '$4 < 5' | awk '{print $1}' | sed s/reads.fq/reads.bam/ | awk '{print "mv specimen_reads/"$1" clusters_small/"$1}' >> ../mv.sh
  cd ..
  chmod +x mv.sh
  ./mv.sh
#
echo "Polishing and variant calling"
  for j in $(ls specimen_reads/*.fas | sed 's/.fas//' | sed 's/specimen_reads\///')
  do echo "Consensus sequence polishing of $j"
     seqtk sample -s$(echo $RANDOM) specimen_reads/$j.reads.fq 200 | seqkit seq -n > IDs.txt # taking 200 random reads for consensus polishing
     samtools view -N IDs.txt specimen_reads/$j.reads.bam -O BAM -o seq.bam
     rm *.bai
     rm *.fai
     cp specimen_reads/$j.fas cons.fa # dorado polish accepts only *.fasta or *.fa extensions, not *.fas
     samtools faidx cons.fa
     dorado aligner cons.fa seq.bam | samtools sort > alignment.bam
     samtools index alignment.bam
     dorado polish alignment.bam cons.fa --ignore-read-groups --models-directory doradoModels > consmed.fas # dorado polish alignment.bam cons.fa --ignore-read-groups > consmed.fas
     seqkit replace -p $(seqkit seq -n consmed.fas) -r $j consmed.fas > specimen_reads/$j.fas
     echo "Variant calling and consensus sequence creation and polishing of $j"
     if (($(samtools view -c seq.bam) < 5))
     then echo "Fewer than 5 reads for $j"
     else rlen=$(seqkit stats -T specimen_reads/$j.fas | cut -f5 | tail -n+2 | awk '{print int($1*0.9)}') # getting 90% of consensus sequence length to filter out too short reads after mapping
          minimap2 -a --sam-hit-only -x map-ont --secondary=no -t 1 specimen_reads/$j.fas specimen_reads/$j.reads.fq | samtools sort | samtools view -e "rlen>=$rlen" -O BAM > alignment.bam # $rlen works only with double quotes in samtools view -e
          samtools index alignment.bam
          samtools faidx specimen_reads/$j.fas
          freebayes -i --haplotype-length -1 -f specimen_reads/$j.fas alignment.bam > var.vcf # variant calling
          grep '#' var.vcf > varf.vcf # creating a separate file for detected variable positions or SNPs having high quality, first adding only the header containing comments of the original file
          grep -v '#' var.vcf | awk '$6>10' >> varf.vcf # adding only detected SNPs with high quality (score higher than 10)
        if (($(grep -v '#' varf.vcf | wc -l) < 1)) # if no high quality SNPs were detected, do not process the consensus sequence further
        then echo "No SNPs for $j detected"
             cp alignment.bam specimen_reads/$j.bam
             samtools index specimen_reads/$j.bam
        else $HOME/devider -b alignment.bam -v varf.vcf -r specimen_reads/$j.fas -o devider_output -t 4 --preset nanopore-r10 --min-cov 5 -O # keeping only variants supported by at least 5 reads
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
                  seqtk sample -s$(echo $RANDOM) inputvar.fas 100 > inputc.fas # 100 random reads for consensus sequence computation
                  abpoa inputc.fas > constemp.fas # computing consensus sequence or a variant
                  seqkit replace -p $(seqkit seq -n constemp.fas) -r $j.$(echo $v | sed 's/Haplotype://') constemp.fas > cons.fa # renaming consensus sequence of the variant # dorado polish accepts only *.fasta or *.fa extensions, not *.fas
                  seqtk sample -s$(echo $RANDOM) specimen_reads/$j.$(echo $v | sed 's/Haplotype://').reads.fq 200 | seqkit seq -n > IDs.txt
                  samtools view -N IDs.txt specimen_reads/$j.reads.bam -O BAM -o seq.bam
                  rm *.bai
                  rm *.fai
                  samtools faidx cons.fa
                  dorado aligner cons.fa seq.bam | samtools sort > alignment.bam
                  samtools index alignment.bam
                  dorado polish alignment.bam cons.fa --ignore-read-groups --models-directory doradoModels > cons_withPrimers/$j.$(echo $v | sed 's/Haplotype://').fas
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
  rm *.fai
  if (($(ls *.fas | grep -f ../gIDs.txt -o | wc -l) < 1))
  then echo "No single variant sequences"
  else 
    for j in $(ls *.fas | sed 's/.vcluster.*$//' | sort | uniq)
    do cat $j*.fas > ../bcCons/$j.new.fas
    done
  fi
  cd ..
#
  echo "Removing primers from multi variant sequences, creating mapping files for each variant, and adding the variants to corresponding amplicon files in the folder bcCons"
  cd cons_withPrimers
  if (($(ls *.fas | grep -f ../gIDs.txt -o | wc -l) < 1))
  then echo "No sequences with more than one variant"
  else 
     for j in $(ls *.fas | sed 's/.vcluster.*$//' | sort | uniq)
     do echo "Cutting primers"
     cat $j*.fas > withPrimers.fasta
     cutadapt \
       -g file:../primers.fas \
       -e 0.15 \
       --action trim \
       --untrimmed-output primersNotCut.fasta \
       --overlap 15 \
       --minimum-length 400 \
       --cores 3 \
       -o noPrimers.fasta \
       withPrimers.fasta
       find -name "*.fasta" -type 'f' -empty -delete
       cat noPrimers.fasta >> ../bcCons/$j.new.fas
       cat primersNotCut.fasta >> ../bcCons/primersNotCut.fas
     done
  fi
  cd ../bcCons
  for i in $(ls *.new.fas | sed 's/.new.fas//')
  do echo "Removing duplicate consensus sequences of $i"
     seqkit rmdup -s -D $i.IDs.txt $i.new.fas > $i.deduplicated.fas
  done
#
  if (($(ls *.IDs.txt | wc -l) < 1))
  then echo "No duplicate sequences found"
    else for i in $(ls *.IDs.txt | sed 's/.IDs.txt//')
         do echo "Merging read files of identical consensus sequences of $i"
         while read line
           do echo $line | sed 's/ /\t/' | cut -f2 | sed 's/, /\n/g' | awk '{print "mv ../specimen_reads/"$1".reads.fq ../cons/"$1".reads.fq"}' > mv.sh
              chmod +x mv.sh
              ./mv.sh
              toMerge=$(echo $line | sed 's/ /\t/' | cut -f2 | sed 's/, /\n/g' | head -n1)
              cat ../cons/*.fq > ../specimen_reads/$toMerge.reads.fq
              rm ../cons/*.fq
           done < $i.IDs.txt
         done
  fi
#
  cd ..
  cd specimen_reads
  seqkit stats *reads.fq -T | tail -n +2 | awk '{print $1, $4}' | sed 's/.reads.fq//' | sed 's/.vcluster/ vcluster/' | sort -k1,1 -k3,3nr | awk '{print $1"."$2, $0}' > ../bcCons/readCountID.txt
  cd ..
#
