#1.reformat UMIs
for i in `find . -maxdepth 1 -name "*R1.fastq"`; do bsub -q short -W 240 -n 4 -o $i.out -e $i.err "python RNAseq_reformat_umi_fastq.py -l $i -r ${i/R1/R2} -L $i.reformated  -R ${i/R1/R2}.reformated"; done

#2.map and deduplicate ERCC reads
module load python/2.7.9
module load gcc/5.1.0
module load htslib/1.3
module load libcurl/7.37.0
module load samtools/0.1.19
module load bowtie/1.0.0
pip install --user pysam
bowtie-build ercc.fa ercc
for i in `find . -maxdepth 1 -name "*R1.fastq.reformated"`; do /project/umw_phil_zamore/common/pipelines/ms_short_64G "bowtie -S -I 100 -X 600 --fr ercc -1 $i -2 ${i/R1/R2} > $i.ercc.sam; samtools view -Sbt ercc.fa.fai $i.ercc.sam >  $i.ercc.bam; samtools sort $i.ercc.bam $i.ercc.sorted; python RNAseq_umi_mark_duplicates.py -f $i.ercc.sorted.bam -p 4; samtools view -b -F 0x400 $i.ercc.deumi.sorted.bam > $i.ercc.deumi.sorted.bam.dedup; samtools view $i.ercc.deumi.sorted.bam.dedup > $i.ercc.deumi.sorted.bam.dedup.sam; rm $i.ercc.deumi.sorted.bam.dedup $i.ercc.deumi.sorted.bam.bai $i.ercc.deumi.sorted.bam $i.ercc.sorted.bam.bai $i.ercc.sorted.bam $i.ercc.bam $i.ercc.sam"; done
for i in `find . -maxdepth 1 -name "*.ercc.deumi.sorted.bam.dedup.sam"`;  do echo $i; cut -f3 $i | sort | uniq -c |awk '{print $2"\t"$1}' > temp.tmp; awk -v f="$i" '{a+=$2/2}END{print f"\t"a}' temp.tmp >> ercc.reads; done; rm temp.tmp

#3.align with STAR with piPipes
for i in `find . -maxdepth 1 -name "*R1.fastq.reformated"`; do /project/umw_phil_zamore/common/pipelines/ms "/project/umw_phil_zamore/common/pipelines/piPipes/piPipes rna -l $i -r ${i/R1/R2} -g mm10 -o $i.pipipes -c 8"; done

#4a.in bam files from piPipes mark duplicates
module load python/2.7.9
module load gcc/5.1.0
module load htslib/1.3
module load libcurl/7.37.0
module load samtools/0.1.19
pip install --user pysam
for i in $(find . -maxdepth 1 -name "*.mm10.sorted.bam"); do /project/umw_phil_zamore/common/pipelines/ms  "python RNAseq_umi_mark_duplicates.py -f $i -p 8"; done

#4b.deduplicate bams
for i in $(find . -maxdepth 1 -name "*.deumi.sorted.bam"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "samtools view -b -F 0x400 $i > $i.dedup"; done

#4c.sort deduplicated bams by chrom pos
for i in $(find . -maxdepth 1 -name "*.deumi.sorted.bam.dedup"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "samtools sort $i $i.srt"; done

#5.calculate the number of all mapped reads from .bam
module load bedtools/2.26.0
for i in $(find . -maxdepth 1 -name "*.bam"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "bedtools bamtobed -bed12 -tag NH -i $i | awk -v f=$i '{a[\$4]++}END{for (j in a){total++};print f\"\\t\"total/2}' > $i.all.mapped.reads.txt"; done

#6a.bedgraph: seq.depth normalized coverage
module load bedtools/2.26.0
for i in $(find . -name "*.deumi.sorted.bam.dedup.srt.bam"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "bedtools bamtobed -bed12 -tag NH -i $i | awk 'BEGIN{FS=OFS=\"\\t\"}(\$5==1){if(substr(\$4,length(\$4))==1){\$6=(\$6==\"+\"?\"-\":\"+\")}; print \$0}' > $i.unique.bed12; normscale=\$(awk '{a[\$4]++}END{for (j in a){total++};print 2000000/total}' $i.unique.bed12); bedtools genomecov -scale \$normscale -split -bg -strand + -i $i.unique.bed12 -g mm10.genome > $i.plus; bedtools genomecov -scale \$normscale -split -bg -strand - -i $i.unique.bed12 -g mm10.genome > $i.temp; awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t-\"\$4}' $i.temp > $i.minus ;cat $i.plus $i.minus | sort -k1,1 -k2,2n > $i.bedgraph; rm $i.minus $i.plus $i.temp $i.unique.bed12"; done

#6b.convert bed graphs to bigwigs
for i in $(find . -maxdepth 1 -name "*.bedgraph"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "awk '(\$4>0)' $i > $i.plus; awk '(\$4<0)' $i > $i.minus; /project/umw_phil_zamore/common/pipelines/piPipes/bin/bedGraphToBigWig $i.plus mm10.genome $i.plus.bigWig; /project/umw_phil_zamore/common/pipelines/piPipes/bin/bedGraphToBigWig $i.minus mm10.genome $i.minus.bigWig; rm $i.plus $i.minus"; done

#7.STRINGTIE tpm, fpkm to create the following ctab files for 5P-seq_1_initial_processing.sh
WT1_t_data.ctab
WT2_t_data.ctab
WT3_t_data.ctab
WT4_t_data.ctab
WT5_t_data.ctab
WT6_t_data.ctab

module load stringtie/1.3.4
for i in $(find . -maxdepth 1 -name "*.deumi.sorted.bam.dedup.srt.bam"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "stringtie $i  --rf -o $i.stringtie_mm10 -b $i.ctab_mm10 -p 20 -G mm10.38.92.gtf -A $i.gene_abund.tab_mm10 -e"; done