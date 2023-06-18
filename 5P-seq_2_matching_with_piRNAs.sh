#USE "1tr" files for finding 5P-seq reads explained by various patterns of complementarity with piRNAs from 12.inter.2555.1ppm.realbed.piRNAs from SmallRNAseq_2_prepare_data_for_matching_with_5P-seq.sh
#to save time and disk space, these matches were found for only 4 out of 16 combinations of WT/pi2pi9pi17 datasets (WT1/mut1, WT2/mut1, WT3/mut1, WT4/mut1), and /mut2, /mut3, /mut4 combinations were added later (see 5P-seq_3_calculating_fraction_cleaved.sh)

#1.CONTIGUOUS COMPLEMENTARITY from g2
#1a.split 12.inter.2555.1ppm.realbed.piRNAs into 3,000 piRNAs per chunk for parallel computing
for i in $(find . -name "*.piRNAs"); do echo $i; split -a 4 -d -l 3000 $i $i.xx. ; done 

#1b.find matches with piRNAs
for j in $(find . -maxdepth 1 -name "*.1tr"); do for i in $(find . -maxdepth 1 -name "*.xx.*"); do si=${i:2}; si=${si//.xx./.xy.};  for n in {7..24}; do /project/umw_phil_zamore/common/pipelines/ms_1G_1core_24h "awk -v n=$n 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){pirnan++;pirna[pirnan]=\$0;pirnaseq[pirnan]=\$4}
(FNR<NR){for (i=1;i<=pirnan;i++){if(substr(pirnaseq[i],2,n-1)==substr(\$17,7,n-1)){print \$0\"\\t\"pirna[i]}}}' $i $j > $j.$si.$n.ms"; done; done; done

#1c.merge data after parallel computing
for i in $(find . -maxdepth 1 -name "*.xy.0000.*.ms"); do echo $i; first=${i%.xy.0000.*}; last=${i#*.xy.0000.}; sfirst=${first//.unique.out.m.5.f.both.merged.tr.1tr.12.inter.2555.1ppm.realbed./.unique.out.12wt.2555.}; sfirst=${sfirst//.unique.in.m.5.f.both.merged.tr.1tr.12.inter.2555.1ppm.realbed./.unique.in.12wt.2555.}; cat $first.*.$last >  $sfirst.merged.$last; done
 
#1d.get non-overlapping sets by requiring 1-3 nt after the contiguous match to be non-matching 
for i in $(find . -maxdepth 1 -name "*.ms"); do /project/umw_phil_zamore/common/pipelines/ms_short_1G "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-1]+1,1)!=substr(\$17,5+a[n-1]+1,1))&&(substr(\$29,a[n-1]+2,1)!=substr(\$17,5+a[n-1]+2,1))&&(substr(\$29,a[n-1]+3,1)!=substr(\$17,5+a[n-1]+3,1)){print}' $i > $i.3nt"; done

for i in $(find . -maxdepth 1 -name "*.ms"); do /project/umw_phil_zamore/common/pipelines/ms_short_1G "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-1]+1,1)!=substr(\$17,5+a[n-1]+1,1))&&(substr(\$29,a[n-1]+2,1)!=substr(\$17,5+a[n-1]+2,1)){print}' $i > $i.2nt"; done

for i in $(find . -maxdepth 1 -name "*.ms"); do /project/umw_phil_zamore/common/pipelines/ms_short_1G "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-1]+1,1)!=substr(\$17,5+a[n-1]+1,1)){print}' $i > $i.1nt"; done

#OUTPUT by COLUMN (1-25: 5' P site, 26-36: piRNA)
 1-5P site _chr
 2-5P site _start
 3-5P site _end
 4-5P site _log2change
 5-wt_count,mut_count,wt_ppm,mut_ppm
 6-5P site _strand
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
 26-chr piRNA
 27-start piRNA
 28-end piRNA
!29-piRNA sequence
 30-12wt abundance
 31-piRNA strand
 32-cluster
 33-piRNA gene
 34-piRNA RMSK
 35-piRNA RMSK strand
 36-piRNA 9reps of triple data abundance




#2.Single mismatch in a contiguous COMPLEMENTARITY from g2
#2a.split 12.inter.2555.1ppm.realbed.piRNAs into 50 piRNA per chunk for parallel computing
for i in $(find . -name "*.piRNAs"); do echo $i; split -a 4 -d -l 50 $i $i.xx. ; done 

#2b.find matches with piRNAs
for j in $(find . -maxdepth 1 -name "*.1tr"); do sj=${j//.unique.out.m.5.f.both.merged.tr.1tr/.out}; for i in $(find . -maxdepth 1 -name "*.xx.*"); do si=${i:2}; si=${si//.xx./.xy.}; si=${si//12.inter.2555.1ppm.realbed.piRNAs/12.inter.2555}; for n in {11..21}; do /project/umw_phil_zamore/common/pipelines/ms_1G_1core_24h "awk -v n=$n -v name=$sj.$si 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){pirnan++;pirna[pirnan]=\$0;pirnaseq[pirnan]=\$4}
(FNR<NR){for (i=1;i<=pirnan;i++){if(substr(pirnaseq[i],2,n-1)==substr(\$17,7,n-1)){print \$0\"\\t\"pirna[i]> name\".\"n\".ms\"}else{for (k=2;k<=n;k++){if(substr(pirnaseq[i],2,k-2)\".\"substr(pirnaseq[i],k+1,n-k)==substr(\$17,7,k-2)\".\"substr(\$17,6+k,n-k)){print \$0\"\\t\"pirna[i] > name\".\"n\".ms.\"k\".mm\"}}}}}' $i $j"; done; done; done

#2c.merge data after parallel computing
for i in $(find . -maxdepth 1 -name "*.xx.0000"); do si=${i:2}; si=${si//.xx./.xy.}; si=${si//12.inter.2555.1ppm.realbed.piRNAs/12.inter.2555}; for j in $(find . -maxdepth 1 -name "*.1tr"); do sj=${j//.unique.out.m.5.f.both.merged.tr.1tr/.out}; for n in {11..21}; do fname=$sj.$si.$n.ms; echo $fname; first=${fname%.xy.0000.*}; last=${fname#*.xy.0000.}; cat $first.*.$last > $first.merged.$last; for m in $(seq 2 $n); do fname=$sj.$si.$n.ms.$m.mm; echo $fname; 
first=${fname%.xy.0000.*}; last=${fname#*.xy.0000.}; cat $first.*.$last > $first.merged.$last; done; done; done; done

#2d.get non-overlapping sets by requiring 1-3 nt after the contiguous match to be non-matching
for i in $(find . -maxdepth 1 -name "*.mm"); do /project/umw_phil_zamore/common/pipelines/ms_short_1G "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-3]+1,1)!=substr(\$17,5+a[n-3]+1,1))&&(substr(\$29,a[n-3]+2,1)!=substr(\$17,5+a[n-3]+2,1))&&(substr(\$29,a[n-3]+3,1)!=substr(\$17,5+a[n-3]+3,1)){print}' $i > $i.3nt"; done
for i in $(find . -maxdepth 1 -name "*.mm"); do /project/umw_phil_zamore/common/pipelines/ms_short_1G "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-3]+1,1)!=substr(\$17,5+a[n-3]+1,1))&&(substr(\$29,a[n-3]+2,1)!=substr(\$17,5+a[n-3]+2,1))' $i > $i.2nt"; done
for i in $(find . -maxdepth 1 -name "*.mm"); do /project/umw_phil_zamore/common/pipelines/ms_short_1G "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-3]+1,1)!=substr(\$17,5+a[n-3]+1,1)){print}' $i > $i.1nt"; done

#OUTPUT by COLUMN (1-25: 5' P site, 26-36: piRNA)
 1-5P site _chr
 2-5P site _start
 3-5P site _end
 4-5P site _log2change
 5-wt_count,mut_count,wt_ppm,mut_ppm
 6-5P site _strand
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
 26-chr piRNA
 27-start piRNA
 28-end piRNA
!29-piRNA sequence
 30-12wt abundance
 31-piRNA strand
 32-cluster
 33-piRNA gene
 34-piRNA RMSK
 35-piRNA RMSK strand
 36-piRNA 9reps of triple data abundance



#3.Mononucleotide bulges: i.e., target insertions and deletions for contiguous COMPLEMENTARITY from g2
#3a.split 12.inter.2555.1ppm.realbed.piRNAs into 50 piRNAs per chunk for parallel computing
for i in $(find . -name "*.piRNAs"); do echo $i; split -a 4 -d -l 50 $i $i.xx. ; done 

#3b.find matches with piRNAs with MONObuldgesINguide (=target deletion)
for j in $(find . -maxdepth 1 -name "*.1tr"); do sj=${j//.unique.out.m.5.f.both.merged.tr.1tr/.out}; for i in $(find . -maxdepth 1 -name "*.xx.*"); do si=${i:2}; si=${si//.xx./.xy.}; si=${si//12.inter.2555.1ppm.realbed.piRNAs/12.inter.2555}; for n in {11..20}; do /project/umw_phil_zamore/common/pipelines/ms_1G_1core_24h "awk -v n=$n -v name=$sj.$si 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){pirnan++;pirna[pirnan]=\$0;pirnaseq[pirnan]=\$4}
(FNR<NR){for (i=1;i<=pirnan;i++){if(substr(pirnaseq[i],2,n-1)!=substr(\$17,7,n-1)){for (k=2;k<n;k++){if(substr(pirnaseq[i],2,k-1)\"\"substr(pirnaseq[i],k+2,n-k)==substr(\$17,7,n-1)){print \$0\"\\t\"pirna[i] > name\".\"n\".ms.\"k\".gb\"}}}}}' $i $j"; done; done; done
#3c.merge tb data after parallel computing
for i in $(find . -maxdepth 1 -name "*.xx.0000"); do si=${i:2}; si=${si//.xx./.xy.}; si=${si//12.inter.2555.1ppm.realbed.piRNAs/12.inter.2555}; for j in $(find . -maxdepth 1 -name "*.1tr"); do sj=${j//.unique.out.m.5.f.both.merged.tr.1tr/.out}; for n in {11..20}; do  let "nn = $n - 1"; for m in $(seq 2 $nn); do fname=$sj.$si.$n.ms.$m.gb; echo $fname; 
first=${fname%.xy.0000.*}; last=${fname#*.xy.0000.}; cat $first.*.$last > $first.merged.$last; done; done; done; done

#3d.find matches with piRNAs with MONObuldgesINtarget (=target insertion)
for j in $(find . -maxdepth 1 -name "*.1tr"); do sj=${j//.unique.out.m.5.f.both.merged.tr.1tr/.out}; for i in $(find . -maxdepth 1 -name "*.xx.*"); do si=${i:2}; si=${si//.xx./.xy.}; si=${si//12.inter.2555.1ppm.realbed.piRNAs/12.inter.2555}; for n in {11..20}; do /project/umw_phil_zamore/common/pipelines/ms_1G_1core_24h "awk -v n=$n -v name=$sj.$si 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){pirnan++;pirna[pirnan]=\$0;pirnaseq[pirnan]=\$4}
(FNR<NR){for (i=1;i<=pirnan;i++){if(substr(pirnaseq[i],2,n-1)!=substr(\$17,7,n-1)){for (k=2;k<n;k++){if(substr(pirnaseq[i],2,n-1)==substr(\$17,7,k-1)\"\"substr(\$17,7+k,n-k)){print \$0\"\\t\"pirna[i] > name\".\"n\".ms.\"k\".tb\"}}}}}' $i $j"; done; done; done
#3e.merge tb data after parallel computing
for i in $(find . -maxdepth 1 -name "*.xx.0000"); do si=${i:2}; si=${si//.xx./.xy.}; si=${si//12.inter.2555.1ppm.realbed.piRNAs/12.inter.2555}; for j in $(find . -maxdepth 1 -name "*.1tr"); do sj=${j//.unique.out.m.5.f.both.merged.tr.1tr/.out}; for n in {11..20}; do  let "nn = $n - 1"; for m in $(seq 2 $nn); do fname=$sj.$si.$n.ms.$m.tb; echo $fname; 
first=${fname%.xy.0000.*}; last=${fname#*.xy.0000.}; cat $first.*.$last > $first.merged.$last; done; done; done; done

#3f.get non-overlapping sets by requiring 1 nt or 3 nt after the contiguous match to be non-matching
for i in $(find . -maxdepth 1 -name "*.gb"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-3]+2,1)!=substr(\$17,5+a[n-3]+1,1))&&(substr(\$29,a[n-3]+3,1)!=substr(\$17,5+a[n-3]+2,1))&&(substr(\$29,a[n-3]+4,1)!=substr(\$17,5+a[n-3]+3,1)){print}' $i > $i.3nt"; done
for i in $(find . -maxdepth 1 -name "*.gb"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-3]+2,1)!=substr(\$17,5+a[n-3]+1,1)){print}' $i > $i.1nt"; done

for i in $(find . -maxdepth 1 -name "*.tb"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-3]+1,1)!=substr(\$17,5+a[n-3]+2,1))&&(substr(\$29,a[n-3]+2,1)!=substr(\$17,5+a[n-3]+3,1))&&(substr(\$29,a[n-3]+3,1)!=substr(\$17,5+a[n-3]+4,1)){print}' $i > $i.3nt"; done
for i in $(find . -maxdepth 1 -name "*.tb"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-3]+1,1)!=substr(\$17,5+a[n-3]+2,1)){print}' $i > $i.1nt"; done

#OUTPUT by COLUMN (1-25: 5' P site, 26-36: piRNA)
 1-5P site _chr
 2-5P site _start
 3-5P site _end
 4-5P site _log2change
 5-wt_count,mut_count,wt_ppm,mut_ppm
 6-5P site _strand
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
 26-chr piRNA
 27-start piRNA
 28-end piRNA
!29-piRNA sequence
 30-12wt abundance
 31-piRNA strand
 32-cluster
 33-piRNA gene
 34-piRNA RMSK
 35-piRNA RMSK strand
 36-piRNA 9reps of triple data abundance




#4.Contigous matches from g3,g4,g5,g6,..., g14, g15.
#4a.split 12.inter.2555.1ppm.realbed.piRNAs into 500 piRNAs per chunk for parallel computing
for i in $(find . -name "*.piRNAs"); do echo $i; split -a 4 -d -l 500 $i $i.xx. ; done 

#4b.find matches with piRNAs
for j in $(find . -maxdepth 1 -name "*.1tr"); do sj=${j//.unique.out.m.5.f.both.merged.tr.1tr/};  for i in $(find . -maxdepth 1 -name "*.xx.*"); do si=${i:2}; si=${si//.xx./.xy.};  for s in {2..15}; do  let "nn = 24 - $s + 1"; for (( n=6; n<=nn; n++)); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core_24h "awk -v s=$s -v n=$n 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){pirnan++;pirna[pirnan]=\$0;pirnaseq[pirnan]=\$4}
(FNR<NR){for (i=1;i<=pirnan;i++){if(substr(pirnaseq[i],s,n)==substr(\$17,s+5,n)){print \$0\"\\t\"pirna[i]}}}' $i $j > $sj.$si.$s.$n.ms"; done; done; done; done

#4c.merge data after parallel computing
for i in $(find . -maxdepth 1 -name "*.xy.0000.*.ms"); do echo $i; first=${i%.xy.0000.*}; last=${i#*.xy.0000.}; sfirst=${first//.12.inter.2555.1ppm.realbed.piRNAs/}; cat $first.*.$last >  $sfirst.merged.$last; rm $first.xy.*.$last; done

#4d.get non-overlapping sets for g5,g6,g7,g8-... by requiring 1nt/3nts to be not paired before and after the contiguous match
for i in $(find . -maxdepth 1 -name "*.ms"); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core_24h "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-2]+a[n-1],1)!=substr(\$17,5+a[n-2]+a[n-1],1))&&(substr(\$29,a[n-2]+a[n-1]+1,1)!=substr(\$17,5+a[n-2]+a[n-1]+1,1))&&(substr(\$29,a[n-2]+a[n-1]+2,1)!=substr(\$17,5+a[n-2]+a[n-1]+2,1))&&(substr(\$29,a[n-2]-1,1)!=substr(\$17,5+a[n-2]-1,1))&&(substr(\$29,a[n-2]-2,1)!=substr(\$17,5+a[n-2]-2,1))&&(substr(\$29,a[n-2]-3,1)!=substr(\$17,5+a[n-2]-3,1)){print}'  $i > $i.3nt"; done
for i in $(find . -maxdepth 1 -name "*.??.ms"); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core_24h "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-2]+a[n-1],1)!=substr(\$17,5+a[n-2]+a[n-1],1))&&(substr(\$29,a[n-2]-1,1)!=substr(\$17,5+a[n-2]-1,1)){print}'  $i > $i.1nt"; done

#4e.get non-overlapping sets for g3-... by requiring 1nt before and 3nts after the contiguous stretch not to be paired
for i in $(find . -maxdepth 1 -name "*.merged.3.*.ms"); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core_24h "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-2]+a[n-1],1)!=substr(\$17,5+a[n-2]+a[n-1],1))&&(substr(\$29,a[n-2]+a[n-1]+1,1)!=substr(\$17,5+a[n-2]+a[n-1]+1,1))&&(substr(\$29,a[n-2]+a[n-1]+2,1)!=substr(\$17,5+a[n-2]+a[n-1]+2,1))&&(substr(\$29,a[n-2]-1,1)!=substr(\$17,5+a[n-2]-1,1)){print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$29}'  $i > $i.3nt"; done


#4f.get non-overlapping sets for g4-... by requiring 2nt before and 3nts after the contiguous stretch not to be paired
#get non-overlapping sets for g4contig KEEP ONLY 7 fields and require 2nt before and 3nts after the contiguous stretch not paired
for i in $(find . -maxdepth 1 -name "*.merged.4.*.ms"); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core_24h "awk -v f=\"$i\" 'BEGIN{FS=OFS=\"\\t\";n=split (f,a,\".\")}
(substr(\$29,a[n-2]+a[n-1],1)!=substr(\$17,5+a[n-2]+a[n-1],1))&&(substr(\$29,a[n-2]+a[n-1]+1,1)!=substr(\$17,5+a[n-2]+a[n-1]+1,1))&&(substr(\$29,a[n-2]+a[n-1]+2,1)!=substr(\$17,5+a[n-2]+a[n-1]+2,1))&&(substr(\$29,a[n-2]-1,1)!=substr(\$17,5+a[n-2]-1,1))&&(substr(\$29,a[n-2]-2,1)!=substr(\$17,5+a[n-2]-2,1)){print}'  $i > $i.3nt"; done

#OUTPUT by COLUMN (1-25: 5' P site, 26-36: piRNA)
 1-5P site _chr
 2-5P site _start
 3-5P site _end
 4-5P site _log2change
 5-wt_count,mut_count,wt_ppm,mut_ppm
 6-5P site _strand
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
 26-chr piRNA
 27-start piRNA
 28-end piRNA
!29-piRNA sequence
 30-12wt abundance
 31-piRNA strand
 32-cluster
 33-piRNA gene
 34-piRNA RMSK
 35-piRNA RMSK strand
 36-piRNA 9reps of triple data abundance
