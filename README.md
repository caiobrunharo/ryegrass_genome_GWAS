# This repository is associated with the genome assembly, annotation, and mapping of paraquat resistance in _Lolium multiflorum_
This repository stores scripts and raw data associated with the ryegrass genome assembly and GWAS of paraquat resistance 

## Assembly of Lolium multiflorum genome with hifiasm
hifiasm -o ryegrass_l1.asm -t36 --h1 HiC_1.fq.gz --h2 HiC_2.fq.gz cell1.fastq.gz cell2.fastq.gz  cell3.fastq.gz  cell4.fastq.gz cell5.fastq.gz cell6.fastq.gz cell7.fastq.gz cell8.fastq.gz  cell9.fastq.gz  cell10.fastq.gz  cell11.fastq.gz --primary -l1 -D 10 

### purge_dups was run on the initial assembly
for i in cell1.fastq.gz cell2.fastq.gz  cell3.fastq.gz  cell4.fastq.gz cell5.fastq.gz cell6.fastq.gz  cell7.fastq.gz cell8.fastq.gz  cell9.fastq.gz  cell10.fastq.gz  cell11.fastq.gz; do
minimap2 -xasm20 -t8 -I3g ryegrass_l1.asm.hic.hap1.p_ctg.fa  ${i} | gzip -c - > ${i}.paf.gz
done

pbcstat ryegrass_cell*.paf.gz

split_fa ryegrass_l1.asm.hic.hap1.p_ctg.fa  > ryegrass_l1.asm.fa.split

minimap2 -xasm5 -I3g -DP ryegrass_l1.asm.fa.split ryegrass_l1.asm.fa.split -t20 | gzip -c - > ryegrass_l1.asm.fa.split.self.paf.gz

calcuts PB.stat > cutoffs 2>calcults.log

python3 hist_plot.py -c cutoffs PB.stat PB.cov.png

purge_dups -a 50 -2 -T cutoffs_manual -c PB.base.cov ryegrass_l1.asm.fa.split.self.paf.gz >  dups_mancut.bed 2> purge_dups_mancut.log

get_seqs dups_mancut.bed ryegrass_l1.asm.hic.hap1.p_ctg.fa  -p cutoffs 

### HiC folllowing the Arima pipeline. https://github.com/ArimaGenomics/mapping_pipeline

bwa mem -t 36 purged_assembly.fa fixed1.fq | samtools view -@ 36 -Sb - > alignment_R1.bam

bwa mem -t 36 purged_assembly.fa fixed2.fq | samtools view -@ 36 -Sb - > alignment_R2.bam

samtools view -h alignment_R1.bam | perl filter_five_end.pl | samtools view -Sb - > alignment_filtered_R1.bam

samtools view -h alignment_R2.bam | perl filter_five_end.pl | samtools view -Sb - > alignment_filtered_R2.bam

perl two_read_bam_combiner.pl alignment_filtered_R1.bam alignment_filtered_R2.bam samtools 10 | samtools view -bS -t purged_assembly.fa.fai - | samtools sort -@ 36 -o paired.bam -

java -Xmx50G -Djava.io.tmpdir=./ -jar ./picard.jar AddOrReplaceReadGroups INPUT=paired.bam OUTPUT=paired_RG.bam ID=CKDL210022428 LB=CKDL210022428 SM=hic PL=ILLUMINA PU=none

java -Xmx50G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=./ -jar ./picard.jar MarkDuplicates INPUT=paired_RG.bam OUTPUT=paired_RG_duplicates.bam METRICS_FILE=metrics.txt TMP_DIR=./ ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

### Scaffolding with YAHS
yahs purged_assembly.fa  paired_RG_duplicates_sorted.bam

juicer pre -a -o out_JBAT ./yahs.out.bin ./yahs.out_scaffolds_final.agp ./purged_assembly.fa.fai  > out_JBAT.log 2>&1

juicer pre ./yahs.out.bin yahs.out_scaffolds_final.agp ./purged_assembly.fa.fai   | sort -k2,2d -k6,6d -T ./ --parallel=16 -S500G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

(java -jar -Xmx50G juicer_tools.1.9.9_jcuda.0.8.jar  pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)
#### Final assembly was denominated LOLMU.fa

## QTL-seq
### Trimming raw reads with HTStream
hts_Stats -N 'RawReads stats' -L ./${sample}_htsStats.json \
	-1 ${sample}_R1.fastq.gz \
	-2 ${sample}_R2.fastq.gz | \
      hts_SeqScreener -N 'screen phix' -A ./${sample}_htsStats.json | \
      hts_AdapterTrimmer -N 'trim adapters' -A ./${sample}_htsStats.json | \
      hts_QWindowTrim -N 'trim low qulity bases from ends of reads' -A ./${sample}_htsStats.json | \
      hts_NTrimmer -N 'remove any remanining N characters' -A ./${sample}_htsStats.json | \
      hts_LengthFilter -N 'remove reads < 50bp' -m 50 -n -A ./${sample}_htsStats.json | \
      hts_Stats -N 'final stats' -A ./${sample}_htsStats.json -f ./${sample}.htstream

### mapping with minimap2
minimap2 -ax sr LOLMU.mmi -t10 -I3g \
 ${sample}.htstream_R1.fastq.gz \
 ${sample}.htstream_R2.fastq.gz -R '@RG\tID:'${sample}'\tSM:'${sample} | \
samblaster --removeDups | \
 samtools view -bh -@ 10 - | \
 samtools sort -@ 10 -O bam -o ${sample}.sorted.bam - 

### Merging mapped bam files among susceptible or resistant individuals
samtools merge -@ 20  S10_100.sorted.bam  S10_114.sorted.bam  S10_147.sorted.bam  S10_184.sorted.bam  S10_202.sorted.bam  S10_227.sorted.bam  S10_31.sorted.bam   S10_359.sorted.bam S10_54.sorted.bam S10_80.sorted.bam S10_103.sorted.bam  S10_118.sorted.bam  S10_157.sorted.bam  S10_187.sorted.bam  S10_207.sorted.bam  S10_23.sorted.bam   S10_35.sorted.bam   S10_37.sorted.bam S10_64.sorted.bam  S10_90.sorted.bam S10_110.sorted.bam  S10_128.sorted.bam  S10_161.sorted.bam  S10_191.sorted.bam  S10_213.sorted.bam  S10_230.sorted.bam  S10_352.sorted.bam  S10_4.sorted.bam    S10_69.sorted.bam  S10_95.sorted.bam S10_111.sorted.bam  S10_132.sorted.bam  S10_166.sorted.bam  S10_196.sorted.bam  S10_217.sorted.bam S10_240.sorted.bam S10_354.sorted.bam  S10_42.sorted.bam   S10_75.sorted.bam S10_113.sorted.bam  S10_136.sorted.bam  S10_178.sorted.bam  S10_20.sorted.bam   S10_222.sorted.bam S10_258.sorted.bam S10_356.sorted.bam  S10_50.sorted.bam   S10_8.sorted.bam -o PQ_S.bam

samtools merge -@ 20   R10_119.sorted.bam  R10_156.sorted.bam  R10_188.sorted.bam  R10_220.sorted.bam   R10_29.sorted.bam   R10_60.sorted.bam  R10_78.sorted.bam  R10_94.sorted.bam R10_101.sorted.bam  R10_127.sorted.bam  R10_162.sorted.bam  R10_193.sorted.bam  R10_225.sorted.bam   R10_355.sorted.bam  R10_61.sorted.bam  R10_81.sorted.bam R10_102.sorted.bam  R10_13.sorted.bam   R10_169.sorted.bam  R10_2.sorted.bam    R10_2351.sorted.bam  R10_357.sorted.bam  R10_65.sorted.bam  R10_87.sorted.bam R10_112.sorted.bam  R10_131.sorted.bam  R10_17.sorted.bam   R10_206.sorted.bam  R10_239.sorted.bam   R10_46.sorted.bam   R10_73.sorted.bam  R10_88.sorted.bam R10_117.sorted.bam  R10_154.sorted.bam  R10_177.sorted.bam  R10_219.sorted.bam  R10_28.sorted.bam    R10_47.sorted.bam   R10_76.sorted.bam  R10_93.sorted.bam  -o PQ_R.bam

### Adding read groups, sorting, and indexing R or S bam files
samtools addreplacerg -r ID:PQ_S -r SM:PQ_S -o PQ_S_RG.bam -@ 20 PQ_S.bam

samtools addreplacerg -r ID:PQ_R -r SM:PQ_R -o PQ_R_RG.bam -@ 20 PQ_R.bam

samtools sort -o PQ_S_sorted.bam PQ_S_RG.bam -@ 20

samtools sort -o PQ_R_sorted.bam PQ_R_RG.bam -@ 20

samtools index PQ_S_sorted.bam -@ 20

samtools index PQ_R_sorted.bam -@ 20

### Calling SNPs between R and S groups
freebayes-parallel <(fasta_generate_regions.py LOLMU.fa.fai 1000000) 48 --fasta-reference LOLMU.fa --bam-list bam_list.txt --use-best-n-alleles 4 --populations populations.txt -C 5 -m 10 --min-coverage 5  > output.vcf

### creating a table file that can be input into qtlseqR
gatk CreateSequenceDictionary -R LOLMU.fa

gatk VariantsToTable -V output.vcf -R LOLMU.fa -F CHROM -F POS -F REF -F ALT -F DP -F GQ -GF DP -GF AD -GF GQ -GF PL -GF RO -GF AO -O PW_QTL_RvsS.table

### qtlseqR was performed following tutorials available on https://github.com/bmansfeld/QTLseqr.


## Genome-wide association study
### Trimming raw reads with HTSTream
hts_Stats -N 'RawReads stats' -L ./${sample}_htsStats.json \
	-1 ${sample}_R1.fastq.gz \
	-2 ${sample}_R2.fastq.gz | \
      hts_SeqScreener -N 'screen phix' -A ./${sample}_htsStats.json | \
      hts_AdapterTrimmer -N 'trim adapters' -A ./${sample}_htsStats.json | \
      hts_QWindowTrim -N 'trim low qulity bases from ends of reads' -A ./${sample}_htsStats.json | \
      hts_NTrimmer -N 'remove any remanining N characters' -A ./${sample}_htsStats.json | \
      hts_LengthFilter -N 'remove reads < 50bp' -m 50 -n -A ./${sample}_htsStats.json | \
      hts_Stats -N 'final stats' -A ./${sample}_htsStats.json -f ./${sample}.htstream

### mapping with minimap2
minimap2 -ax sr ref.mmi -t18 -I3g \
 ${sample}.htstream_R1.fastq.gz \
 ${sample}.htstream_R2.fastq.gz -R '@RG\tID:'${sample}'\tSM:'${sample} | \
 samtools view -bh -@ 18 - | \
 samtools sort -@ 18 -O bam -o ${sample}.sorted.bam - 

### removing duplicates with MarkDuplicates
java -Xmx400G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=./ -jar picard.jar MarkDuplicates INPUT=${sample}.sorted.bam OUTPUT=${sample}_dups.bam METRICS_FILE=metrics.txt TMP_DIR=./ ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

### SNP calling with freebayes
freebayes-parallel <(fasta_generate_regions.py yahs.out_scaffolds_final.fa 1000000) 20 --fasta-reference yahs.out_scaffolds_final.fa --bam-list bam.fofn --use-best-n-alleles 4 --populations population.txt -C 5 -m 10 --min-coverage 5  > output.vcf

bcftools view -i 'F_MISSING < 0.25 || (INFO/DP >= 10 && INFO/DP <= 250) || QUAL >= 20' --threads 48 output.vcf.gz -Oz -o output_075_dp5_qual20.vcf.gz

### The final SNP dataset was used for GWAS with GAPIT's implementation of the ECMLM model in R
myGAPIT_ECMLM <- GAPIT( Y=myY, 
                        G=myGenotype, 
                        PCA.total=4,
                        Multiple_analysis = F,
                        model = "ECMLM")

## RNA-seq analysis
### Raw reads were processed with bbduk
bbduk.sh in=${sample}_L003_R1_001.fastq.gz out=${sample}_trimmed.fastq.gz ref=trueseq.fa literal=AAAAAAAAAAAAAAAAAA k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 stats=${sample}_stats

### mapping with STAR
STAR \
--runThreadN 48 \
--runMode genomeGenerate \
--sjdbGTFtagExonParentTranscript Parent \
--genomeDir ./ \
--genomeFastaFiles LOLMU.fa \
--sjdbGTFfile LOLMU.gff \
--sjdbOverhang 99

STAR \
--runThreadN 20 \
--genomeDir ./genome \
--readFilesIn ${trim}_trimmed.fastq.gz \
--readFilesCommand zcat \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${trim}_aligned \
--quantMode GeneCounts

### generation of a read count table 
for sample in $(cut -f1 read_list.txt); do cat ${sample} | tail -n +5  | cut -f3 > ./counts/tmp/${sample}.count ; done

tail -n +5 sampleA_alignedReadsPerGene.out.tab | cut -f1 > counts/tmp/geneids.txt

paste counts/tmp/geneids.txt counts/tmp/*.count > counts/tmp/tmp.out

cat <(cat read_list.txt  | sort | paste -s) counts/tmp/tmp.out > counts/rnaseq_count_table.txt

### RNA-seq analysis followed tutorial available at https://ucdavis-bioinformatics-training.github.io/2020-mRNA_Seq_Workshop/







