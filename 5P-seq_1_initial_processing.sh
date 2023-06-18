#1.find and remove 5' UMI part
for i in `find . -maxdepth 1 -name "*.R1.fastq"`; do /project/umw_phil_zamore/common/pipelines/ms_480G_10days "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR%4==1)&&(FNR==NR){name=\$1;getline;if((index(\$1,\"N\")==0)&&(((substr(\$1,4,3)==\"ATC\")&&(substr(\$1,10,3)==\"AGT\"))||((substr(\$1,4,3)==\"CGA\")&&(substr(\$1,10,3)==\"TAC\")))){a[FNR-1]=1; print name > FILENAME\".rfmtd.fq\"; 
print substr(\$1,16)> FILENAME\".rfmtd.fq\"; getline;
print > FILENAME\".rfmtd.fq\"; getline
print substr(\$1,16)> FILENAME\".rfmtd.fq\"}}(FNR<NR)&&(a[FNR]==1)&&(FNR%4==1){print > FILENAME\".rfmtd.fq\" ; getline;print > FILENAME\".rfmtd.fq\";getline;print > FILENAME\".rfmtd.fq\";getline;print > FILENAME\".rfmtd.fq\"}' $i ${i/.R1./.R2.}"; done

#2.trim the sequencing adapters
module load cutadapt/2.9
for i in $(ls *.R1.fastq.rfmtd.fq); do cutadapt -a TGGAATTCTCGGGTGCCAAGG -A GATCGTCGGACTGTAGAACTC -m 20 --pair-filter=any -j 40 --overlap 15 -o ${i/.R1.fastq.rfmtd.fq/.R1.rfmt.trmd.fq} -p ${i/.R1.fastq.rfmtd.fq/.R2.rfmt.trmd.fq} $i ${i/.R1./.R2.} 1> $i.txt; done

#3.to harmonize nextseq and novaseq read length, make R1+R2 64nt+79nt
for i in `find . -name "*.R1.rfmt.trmd.fq"`; do /project/umw_phil_zamore/common/pipelines/ms " awk 'BEGIN{FS=OFS=\"\\t\"}(FNR%4==1){print;getline;print substr(\$1,1,64);getline;print; getline;print substr(\$1,1,64)}' $i > ${i/.R1.rfmt.trmd.fq/.R1.rfmt.trmd.6479.fq} "; done
for i in `find . -name "*.R2.rfmt.trmd.fq"`; do /project/umw_phil_zamore/common/pipelines/ms " awk 'BEGIN{FS=OFS=\"\\t\"}(FNR%4==1){print;getline;print substr(\$1,1,79);getline;print; getline;print substr(\$1,1,79)}' $i > ${i/.R2.rfmt.trmd.fq/.R2.rfmt.trmd.6479.fq} "; done

#4.remove rRNA reads and align onto mm10 with STAR in piPipes
for i in `find . -name "*.R1.rfmt.trmd.6479.fq"`; do /project/umw_phil_zamore/common/pipelines/ms_large "/project/umw_phil_zamore/common/pipelines/piPipes/piPipes deg -l $i -r ${i/.R1./.R2.} -g mm10 -o $i.pipipes.deg -c 20"; done

#5.collapse 5ends for all uniquely mapped reads
for i in $(find . -maxdepth 1  -name "*.unique.bed12"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "awk '{if (\$6==\"+\"){aplus[\$1,\$2]+=\$5}else{aminus[\$1,\$3]+=\$5}}END{for (ij in aplus) {split(ij,indices,SUBSEP);
 i=indices[1];
 j=indices[2]; print i\"\\t\"j\"\\t\"j+1\"\\t\"aplus[i,j]\"\\tna\\t+\"};for (ij in aminus) {split(ij,indices,SUBSEP);
 i=indices[1];
 j=indices[2]; print i\"\\t\"j-1\"\\t\"j\"\\t\"aminus[i,j]\"\\tna\\t-\"}}' $i > $i.1"; done

#6.normalize to seq.depth: make rpm 
for i in $(find . -name "*.1"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "awk '{a[FNR]=\$0; total+=\$4; num=FNR}END{for (i=1;i<=num;i++){split (a[i],b,\"\\t\");print b[1]\"\\t\"b[2]\"\\t\"b[3]\"\\t\"1000000*b[4]/total\"\\t\"b[4]\"\\t\"b[6]}}' $i > $i.rpm"; done

#7.keep only reads outside pachytene clusters
module load bedtools/2.26.0
for i in $(find . -name "*.rpm"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "bedtools intersect -wa -v -a $i -b piRNA.cluster.pachytene.merged.bed6 > $i.out.pch"; done

#8.create all 16 combinations of 4 WT vs 4 p2p9p17 reps
for i in $(find . -maxdepth 1 -name "*WT*.out.pch"); do si=${i:2}; si=${si%%.x_rRNA*}; si=${si%%.unique*}; for j in $(find . -maxdepth 1 -name "*p7*.out.pch" -o -name "*p2p9p17*.out.pch"); do sj=${j:2}; sj=${sj%%.x_rRNA*}; /project/umw_phil_zamore/common/pipelines/ms_short_64G "awk '(FNR==NR){mutrpm[\$1,\$6,\$2]=\$4;mutcount[\$1,\$6,\$2]=\$5;mutcoord[\$1,\$6,\$2]=\$1\"\\t\"\$2\"\\t\"\$3;mutstrand[\$1,\$6,\$2]=\$6}(FNR<NR){if(mutrpm[\$1,\$6,\$2]==\"\"){print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\",0\\t\"\$5\",0\\t\"\$6}else{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\",\"mutrpm[\$1,\$6,\$2]\"\\t\"\$5\",\"mutcount[\$1,\$6,\$2]\"\\t\"\$6}}' $j $i > $si.$sj.unique.out.m"; done; done

#9.filter at 5 reads in WT to speed up downstream steps
for i in $(find . -name "*.m"); do si=${i//unique.bed12.1.rpm.in.pch./}; si=${si//unique.bed12.1.rpm.out.pch./}; /project/umw_phil_zamore/common/pipelines/ms_short_8G "awk '(FS=OFS=\"\\t\"){split (\$5,a,\",\");if((a[1])>=5){split (\$4,a,\",\");if(a[1]==0){\$4=15}else{if(a[2]==0){\$4=-15}else {\$4=log(a[2]/a[1])/log(2)}}; \$5=\$5\",\"a[1]\",\"a[2]; print}}' $i > $si.5.d"; done

#10a. intersect with denovo assembled genes, rmsk, and mm10.38.92.genesNtype
module load bedtools/2.26.0
for i in $(find . -name "*.d"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "bedtools intersect -s -wao -a $i -b 4N.2samples.merged.stringtie.gtf.gene_merged_exons.bed | awk '{a[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6]=\$10}END{for (i in a){print i\"\\t\"a[i]}}' > $i.gdn"; done
#10b
for i in $(find . -name "*.gdn"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "bedtools intersect -wao -a $i -b mm10.rmsk.bed | awk '{if(\$6==\$13){a[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7]=\$11\"\\t+\"}else{a[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7]=\$11\"\\t-\"}}END{for (i in a){print i\"\\t\"a[i]}}' > $i.r"; done
#10c
for i in $(find . -name "*.r"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "bedtools intersect -s -wao -a $i -b mm10.38.92.gtf.exons.geneNnameNtype.bed | awk '{if (\$14==-1){b[1]=\".\";b[2]=\".\"}else{split (\$14,b,\":\")}; a[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9]=\$13\";\"b[1]\";\"b[2]}END{for (i in a){print i\"\\t\"a[i]}}' > $i.g"; done


#11.get GENOMIC -25+15 sequence for REVERSED and -60+60 sequence for not REVERSED
for i in $(find . -name "*.d.gdn.r.g"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "awk '(\$6==\"+\"){if((\$2-60)>-1){
print \$1\"\\t\"\$2-25\"\\t\"\$2+15\"\\t\"\$1\";\"\$2\";\"\$3\";\"\$4\";\"\$5\";\"\$6\";\"\$7\";\"\$8\";\"\$9\";\"\$10\"\\t\"\$5\"\\t-\";
print \$1\"\\t\"\$2-60\"\\t\"\$2+60\"\\t\"\$1\";\"\$2\";\"\$3\";\"\$4\";\"\$5\";\"\$6\";\"\$7\";\"\$8\";\"\$9\";\"\$10\"\\t\"\$5\"\\t+\"}}(\$6==\"-\"){if((\$3-60)>-1){
print \$1\"\\t\"\$3-15\"\\t\"\$3+25\"\\t\"\$1\";\"\$2\";\"\$3\";\"\$4\";\"\$5\";\"\$6\";\"\$7\";\"\$8\";\"\$9\";\"\$10\"\\t\"\$5\"\\t+\";
print \$1\"\\t\"\$3-60\"\\t\"\$3+60\"\\t\"\$1\";\"\$2\";\"\$3\";\"\$4\";\"\$5\";\"\$6\";\"\$7\";\"\$8\";\"\$9\";\"\$10\"\\t\"\$5\"\\t-\"}}' $i > $i.200" ; done
module load bedtools/2.26.0
for i in $(find . -name "*.d.gdn.r.g.200"); do si=${i//DEG_UMI_on_/}; si=${si//.d.gdn.r.g.200/}; /project/umw_phil_zamore/common/pipelines/ms_short_8G "bedtools getfasta -s -name -fi mm10.fa -bed $i -fo $si.f"; done
for i in $(find . -name "*.f"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "awk '(FNR%4==1){print; getline; print;getline; getline; print}' $i > $i.both"; done

#OUTPUT FORMAT: line1->deg_info::ccordinates; line2->revcomp -25+15; line3->sense -60+60

#12.get context sequence from transcripts in addition to genome: use 6reps of WT RNAseq from RNAseq_whole_pipeline.sh
#12a.get mm10.38.92.gtf reference for all transcripts
awk 'BEGIN{FS=OFS="\t"}
(FILENAME=="WT1_t_data.ctab")&&(FNR>1){ab1[$6]=$12}
(FILENAME=="WT2_t_data.ctab")&&(FNR>1){ab2[$6]=$12}
(FILENAME=="WT3_t_data.ctab")&&(FNR>1){ab3[$6]=$12}
(FILENAME=="WT4_t_data.ctab")&&(FNR>1){ab4[$6]=$12}
(FILENAME=="WT5_t_data.ctab")&&(FNR>1){ab5[$6]=$12}
(FILENAME=="WT6_t_data.ctab")&&(FNR>1){ab6[$6]=$12}
(FILENAME=="mm10.38.92.gtf.exons.geneNnameNtype"){type[$1]=$3","$2}
(FILENAME=="mm10.38.92.gtf.exon_coord.fasta.fa.as.form")&&(FNR%2==1){split(substr($1,2),gname,";");name=$1;getline; print name";"(ab1[gname[1]]+ab2[gname[1]]+ab3[gname[1]]+ab4[gname[1]]+ab5[gname[1]]+ab6[gname[1]])/6";"type[gname[2]]";"ab1[gname[1]]":"ab2[gname[1]]":"ab3[gname[1]]":"ab4[gname[1]]":"ab5[gname[1]]":"ab6[gname[1]];print $0}' *.ctab mm10.38.92.gtf.exons.geneNnameNtype mm10.38.92.gtf.exon_coord.fasta.fa.as.form > mm10.38.92.gtf.exon_coord.fasta.fa.as.form.WT_0fpkm

#OUTPUT by COLUMN
1-ENSMUST
2-ENSMUSG
3-refname
4-chr
5-strand
6-exons (1-2,3-4,5-6) 1-based on genomic strand
7-segs (1-2,3-4,5-6) 1-based on transcribed strand (reverse order for - genomic strand)
8-cds start 1-based on transcribed strand (reverse order for - genomic strand)
9-cds end 1-based on transcribed strand (reverse order for - genomic strand)
10-ENSMUST MEAN fpkm abundance in WT37/55/56/61/62/63 by stringtie
11-refname,type
12-ENSMUST fpkm abundance in each WT37/55/56/61/62/63 by stringtie

#$12b. Add sense sequence to transcripts
awk 'BEGIN{FS=OFS="\ "}(FNR==NR)&&(FNR%2==1){name=substr($1,2);getline;a[name]=$1}(FNR<NR)&&(FNR%2==1){name=substr($1,2,index($1,";")-2);print;getline;print;print a[name]}' mm10.38.92.gtf.exon_coord.fasta.fa mm10.38.92.gtf.exon_coord.fasta.fa.as.form.WT_0fpkm > mm10.38.92.gtf.exon_coord.fasta.fa.asNs.form.WT_0fpkm


#12c.split "both" files from 11 into 60,000 lines per file for parallel computing
for i in $(find . -name "*.m.5.f.both"); do echo $i; split -a 4 -d -l 60000 $i $i.xx. ; done 

#12d.match deg reads with transcripts
for i in $(find . -maxdepth 1 -name "*.xx.*"); do si=${i//.xx./.xy.}; /project/umw_phil_zamore/common/pipelines/ms_8G_1core_12h "awk  'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(FNR%3==1){split(substr(\$1,2),dat,\";\");
 split(dat[6],exons,\",\");
 m=split(dat[7],segs,\",\");
 for (i=1;i<=m;i++){split(exons[i],b,\"-\");
 excoor[dat[1],i,1]=b[1]-1;
 excoor[dat[1],i,2]=b[2];
 split(segs[i],b,\"-\");
 segcoor[dat[1],i,1]=b[1]-1;
 segcoor[dat[1],i,2]=b[2]};
 chr[dat[1]]=dat[4];
 strand[dat[1]]=dat[5];
 cdsst[dat[1]]=dat[8]-1;
 cdsend[dat[1]]=dat[9];
 segnum[dat[1]]=m;
 st[dat[1]]=excoor[dat[1],1,1];
 end[dat[1]]= excoor[dat[1],m,2];
 fpkm[dat[1]]= dat[10];
 getline;seq[dat[1]]=toupper(\$1);
 getline;seqse[dat[1]]=toupper(\$1)}(FNR<NR)&&(FNR%3==1){tr=\"\";ndeg=split(substr(\$1,2,index(\$1,\"::\")-2),deg,\";\");
 
 if(deg[6]==\"+\"){for (i in seq){if ((chr[i]==deg[1])&&(deg[6]==strand[i])&&(st[i]+25<=deg[2])&&(deg[2]<=end[i]-15)){ seg=0;for(j=1;j<=segnum[i];j++){if((excoor[i,j,1]<=deg[2])&&(deg[2]<excoor[i,j,2])){seg=j}};if(seg>0){utr=\"\";loc=(segcoor[i,seg,1]+deg[2]-excoor[i,seg,1]);if(cdsend[i]>0){if(loc<cdsst[i]){utr=\"5utr\"};if(cdsend[i]<loc){utr=\"3utr\"};if((loc>=cdsst[i])&&(cdsend[i]>=loc)){utr=\"cds\"}}else{utr=\"nc\"}; if(loc<100){
 tr=tr\";\"i\";\"loc\";\"utr\";\"fpkm[i]\";\"substr(seq[i],segcoor[i,segnum[i],2]-loc-15+1,40)\";\"substr(seqse[i],1,loc+100)\";\"loc\";\"cdsst[i]\";\"cdsend[i]\";\"segcoor[i,segnum[i],2]}else{
 tr=tr\";\"i\";\"loc\";\"utr\";\"fpkm[i]\";\"substr(seq[i],segcoor[i,segnum[i],2]-loc-15+1,40)\";\"substr(seqse[i],loc-100+1,200)\";\"100\";\"cdsst[i]\";\"cdsend[i]\";\"segcoor[i,segnum[i],2]}}}}};
 
 if(deg[6]==\"-\"){for (i in seq){if ((chr[i]==deg[1])&&(deg[6]==strand[i])&&(st[i]+15<=deg[3])&&(deg[3]<=end[i]-25)){
 seg=0;for(j=1;j<=segnum[i];j++){if((excoor[i,j,1]<deg[3])&&(deg[3]<=excoor[i,j,2])){seg=j}};if(seg>0){utr=\"\";loc=(segcoor[i,segnum[i]-seg+1,1]+excoor[i,seg,2]-deg[3]);if(cdsend[i]>0){if(loc<cdsst[i]){utr=\"5utr\"};if(cdsend[i]<loc){utr=\"3utr\"};if((loc>=cdsst[i])&&(cdsend[i]>=loc)){utr=\"cds\"};}else{utr=\"nc\"}; if(loc<100){
 tr=tr\";\"i\";\"loc\";\"utr\";\"fpkm[i]\";\"substr(seq[i],segcoor[i,segnum[i],2]-loc-15+1,40)\";\"substr(seqse[i],1,loc+100)\";\"loc\";\"cdsst[i]\";\"cdsend[i]\";\"segcoor[i,segnum[i],2]}else{
 tr=tr\";\"i\";\"loc\";\"utr\";\"fpkm[i]\";\"substr(seq[i],segcoor[i,segnum[i],2]-loc-15+1,40)\";\"substr(seqse[i],loc-100+1,200)\";\"100\";\"cdsst[i]\";\"cdsend[i]\";\"segcoor[i,segnum[i],2]}}}}};
 
 printf (deg[1]);for (k=2;k<=ndeg;k++){printf (\"\\t\"deg[k])}; 
 getline; printf (\"\\t\"toupper(\$1)); 
 getline;printf (\"\\t\"toupper(\$1)); printf (\"\\t\"tr); printf (\"\\n\")}' mm10.38.92.gtf.exon_coord.fasta.fa.asNs.form.WT_0fpkm $i > $si.tr"; done

#OUTPUT by COLUMN
 1-6 deg fields from .d
 7-denovo gene
 8-RMSK
 9-RMSK_strand
 10-mm10_gene
 11-mm10_gene_type
 12-mm10_gene_name 
 13-genome revcomp for -25+15 
 14-genome sense for -100+100, 
fields of 15: 
1-empty
2-ENSMUST
3-loc of cleavage site in ENSMUST
4-utr_type
5-mean fpkm
6-revcomp of -25+15
7-sense of -100+100
8-loc of cleaveage site in 7
9-start of CDS
10-end of CDS
11-length of ENSMUST

#12e.merge data after parallel computing
for i in $(find . -maxdepth 1 -name "*.xy.0000.tr"); do echo $i; first=${i%.xy.0000.*}; last=${i#*.xy.0000.}; cat $first.*.$last >  $first.merged.$last; done

#12f.keep context only for most abundantly expressed transcript (or if only transcript is found it has to be at 1+fpkm; 1fpkm corresponds to ~3 mRNAs per primary spermatocytes)
for i in $(find . -maxdepth 1 -name "*.tr"); do /project/umw_phil_zamore/common/pipelines/ms_8G_1core_12h "awk  'BEGIN{FS=OFS=\"\\t\"}{n=split (\$15,tr,\";\");if (n>1){max=0;maxn=0;for (i=1;i<=(n-1)/10;i++){if((tr[((i-1)*10)+5]>=1)&&(tr[((i-1)*10)+5]>max)){max=tr[((i-1)*10)+5];maxn=i}}; if(maxn>0){ 
print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$12\"\\t\" tr[((maxn-1)*10)+2]\"\\t\" tr[((maxn-1)*10)+3]\"\\t\"tr[((maxn-1)*10)+4]\"\\t\"tr[((maxn-1)*10)+5]\"\\t\"tr[((maxn-1)*10)+6]\"\\t\"tr[((maxn-1)*10)+7]\"\\t\"tr[((maxn-1)*10)+8]\"\\t\"tr[((maxn-1)*10)+9]\"\\t\"tr[((maxn-1)*10)+10]\"\\t\"tr[((maxn-1)*10)+11]\"\\t\"\$13\"\\t\"\$14\"\\ttranscript\"}else{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$12\"\\t\" 0\"\\t\" 0\"\\t\"0\"\\t\"0\"\\t\"\$13\"\\t\"\$14\"\\t\"0\"\\t\"0\"\\t\"0\"\\t\"0\"\\t\"0\"\\t\"0\"\\tnonexpressedloc\"} }else{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$12\"\\t\" 0\"\\t\" 0\"\\t\"0\"\\t\"0\"\\t\"\$13\"\\t\"\$14\"\\t\"0\"\\t\"0\"\\t\"0\"\\t\"0\"\\t\"0\"\\t\"0\"\\tgenome\"}}' $i > $i.1tr"; done

#OUTPUT by COLUMN
 1-chr
 2-start
 3-end
 4-log2change
 5-wt_count,mut_count,wt_ppm,mut_ppm
 6-strand
 7-denovo gene
 8-RMSK
 9-RMSK_strand
 10-mm10_gene
 11-mm10_gene_type
 12-mm10_gene_name 
 13-ENSMUST,
 14-loc_position in ENSMUST, 
 15-utr_type, 
 16-fpkm, 
!17-revcomp transcriptseq/or/genomeseq for -25+15, 
 18-sense transcriptseq/or/genomeseq for -100+100, 
 19-locpos of cleavage site in sense transcriptseq/or/genomeseq for -100+100
 20-start of CDS
 21-end of CDS
 22-length of ENSMUST
 23-genome revcomp for -25+15 
 24-genome sense for -100+100, 
 25-annotation of the site (genome/transcript/nonexpressedloc)
