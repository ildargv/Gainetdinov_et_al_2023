#0a.create bed file for merged exons for each gene
awk 'BEGIN{FS=OFS="\t"}($3=="exon"){split($9,a,"\""); print $1"\t"$4-1"\t"$5"\t"a[6]"\t"a[2]"\t"$7}' mm10.38.92.gtf > mm10.38.92.gtf.exons.bed
awk '{print > $4".gene"}' mm10.38.92.gtf.exons.bed
module load bedtools/2.26.0
for i in $(find . -name "*.gene"); do bedtools merge -s -c 4 -o distinct -i <(sort -k1,1 -k2,2n $i) > $i.merged; rm $i; done
cat *.merged > mm10.38.92.gtf.genes.bed
rm *.merged
awk '{print $1"\t"$2"\t"$3"\t"$5"\tna\t"$4}' mm10.38.92.gtf.genes.bed > mm10.38.92.gtf.exons.genes.bed

#0b.get genic nonexons (i.e., introns)
num=$(cat mm10.38.92.gtf.exons.genes.bed | wc -l); sort -k4,4 -k2,2n mm10.38.92.gtf.exons.genes.bed | awk -v num="$num" 'BEGIN{FS=OFS="\t"}(FNR==1){cur_name=$4;cur_st=$3;exons=1}(FNR>1){ while (cur_name==$4){print $1"\t"cur_st"\t"$2"\t"$4"\tna\t"$6; cur_st=$3; getline; if((FNR==num)&&(cur_name==$4)){print $1"\t"cur_st"\t"$2"\t"$4"\tna\t"$6}}; cur_name=$4;cur_st=$3;exons=1}' > mm10.38.92.gtf.nonexons.genes.bed

#0c.remove rmsk from nonexons
module load bedtools/2.26.0
bedtools subtract -a mm10.38.92.gtf.nonexons.genes.bed -b mm10.rmsk.bed > mm10.38.92.gtf.nonexons.genes.rmsk_subtracted.bed
bedtools getfasta -name -s -fi mm10.fa -bed mm10.38.92.gtf.nonexons.genes.rmsk_subtracted.bed > mm10.38.92.gtf.nonexons.genes.rmsk_subtracted.fa

#0d.prepare repBase ref files by merging each entry into one line
for i in `find . -name "*.ref"`; do echo $i; num=$(cat $i | wc -l); awk -v num="$num" '(FNR==1){print}(FNR>1){full=""; while (index($1,">")==0){full=full $1; getline; if(FNR==num){full=full $1; print toupper(full);exit}}; print toupper(full); print}' $i > $i.fa; done

#1a.extract all "Mus musculus" from RepBase
for i in `find . -name "*.ref.fa"`; do echo $i; awk '($0~/Mus musculus/){f=substr($1,2)".fasta";print $0 > f; getline; print > f}' $i; done

#1b.get index of kmers in TEs
for l in {6..30}; do for f in `ls *.fasta`; do echo -e $l"\t"$f; awk -v l="$l" 'BEGIN{FS=OFS="\t"}(FNR%2==0){for (i=1;i<=(length($1)-l+1);i++){ind[substr($1,i,l)]++}}END{for (i in ind){print i"\t"ind[i]}}' $f > $f.$l.ind; done; done

#2a.count kmer overlap for exons
for f in `ls *.ind`; do /project/umw_phil_zamore/common/pipelines/ms_1G_1core "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){ind[\$1]=\$2}(FNR<NR)&&(FILENAME==\"mm10.38.92.gtf.t_g_n_type\"){gene[\$1]=\$2;type[\$1]=\$4}(FNR<NR)&&(FILENAME==\"mm10.38.92.gtf.fasta.fa\")&&(FNR%2==1){tr=substr(\$0,2,index(\$1,\"gene=\")-3);gene_total[gene[tr]]=1;getline;for(kmer in ind){if(index(\$1,kmer)>0){gene_hit[gene[tr]]=1;if(gene_type[gene[tr]]!=\"protein_coding\"){gene_type[gene[tr]]=type[tr]}}}}END{for (g in gene_total){total++;if(gene_hit[g]==1){hit++; type_hit[gene_type[g]]++}};print FILENAME\"\\t\"hit\"\\t\"total\"\\t\"hit/total;for (t in type_hit){print t\"\\t\"type_hit[t]}}' $f mm10.38.92.gtf.t_g_n_type mm10.38.92.gtf.fasta.fa > $f.exons"; done

#2b.count kmer overlap for non-exons
for f in `ls *.ind`; do /project/umw_phil_zamore/common/pipelines/ms_1G_1core "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){ind[\$1]=\$2}(FNR<NR)&&(FILENAME==\"mm10.38.92.gtf.t_g_n_type\"){type[\$2]=\$4}(FNR<NR)&&(FILENAME==\"mm10.38.92.gtf.nonexons.genes.rmsk_subtracted.fa\")&&(FNR%2==1){gene=substr(\$0,2,index(\$1,\"::\")-2);gene_total[gene]=1;getline;for(kmer in ind){if(index(\$1,kmer)>0){gene_hit[gene]=1;if(gene_type[gene]!=\"protein_coding\"){gene_type[gene]=type[gene]}}}}END{for (g in gene_total){total++;if(gene_hit[g]==1){hit++; type_hit[gene_type[g]]++}};print FILENAME\"\\t\"hit\"\\t\"total\"\\t\"hit/total;for (t in type_hit){print t\"\\t\"type_hit[t]}}' $f mm10.38.92.gtf.t_g_n_type mm10.38.92.gtf.nonexons.genes.rmsk_subtracted.fa > $f.nonexons"; done

#3.merge kmer counts
awk 'BEGIN{FS=OFS="\t"}(FNR==1){n=split(FILENAME,f,".");dat[f[1],f[n-2]]=$4;dat_all[f[1]]=1}END{for (i in dat_all){print i"\t"dat[i,6]"\t"dat[i,7]"\t"dat[i,8]"\t"dat[i,9]"\t"dat[i,10]"\t"dat[i,11]"\t"dat[i,12]"\t"dat[i,13]"\t"dat[i,14]"\t"dat[i,15]"\t"dat[i,16]"\t"dat[i,17]"\t"dat[i,18]"\t"dat[i,19]"\t"dat[i,20]"\t"dat[i,21]"\t"dat[i,22]"\t"dat[i,23]"\t"dat[i,24]"\t"dat[i,25]"\t"dat[i,26]"\t"dat[i,27]"\t"dat[i,28]"\t"dat[i,29]"\t"dat[i,30]}}' *.exons > kmers.exons.txt
awk 'BEGIN{FS=OFS="\t"}(FNR==1){n=split(FILENAME,f,".");dat[f[1],f[n-2]]=$4;dat_all[f[1]]=1}END{for (i in dat_all){print i"\t"dat[i,6]"\t"dat[i,7]"\t"dat[i,8]"\t"dat[i,9]"\t"dat[i,10]"\t"dat[i,11]"\t"dat[i,12]"\t"dat[i,13]"\t"dat[i,14]"\t"dat[i,15]"\t"dat[i,16]"\t"dat[i,17]"\t"dat[i,18]"\t"dat[i,19]"\t"dat[i,20]"\t"dat[i,21]"\t"dat[i,22]"\t"dat[i,23]"\t"dat[i,24]"\t"dat[i,25]"\t"dat[i,26]"\t"dat[i,27]"\t"dat[i,28]"\t"dat[i,29]"\t"dat[i,30]}}' *.nonexons > kmers.nonexons.txt
