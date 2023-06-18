#1.to save time and disk space, first, use 4 wt "1tr" files from 5P-seq_1_initial_processing.sh to merge all 5P-5'ends in all WT in all 4 reps at >=0.1ppm, run all to all comparisons of piRNAs and 5P'end sequences, then create all 16+16 permutations of four wt vs four p2p9p17 or p7 mutant datasets
#1a.merge
/project/umw_phil_zamore/common/pipelines/ms_short_256G_8cores "awk 'BEGIN{FS=OFS=\"\\t\"}{split(\$5,a,\",\");if(a[3]>=0.1){comb[\$1\"\\t\"\$2\"\\t\"\$6]=1}}END{for (i in comb){print i}}' *.1tr > all_deg_sites_01ppm_from_1tr"
#1b.add info to merged data
/project/umw_phil_zamore/common/pipelines/ms_short_256G_8cores "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){site[\$1\"\\t\"\$2\"\\t\"\$3]=1}(FNR<NR)&&(site[\$1\"\\t\"\$2\"\\t\"\$6]==1){print;site[\$1\"\\t\"\$2\"\\t\"\$6]++}' all_deg_sites_01ppm_from_1tr *.1tr > all_deg_sites_01ppm_from_1tr.info"
#OUTPUT by COLUMN
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
 
#2.add mean piRNA abundance from 3 reps of p7 mutants and calculate mean  piRNA abundance from 12 reps of wt animals
awk 'BEGIN{FS=OFS="\t"}(FILENAME~/realbed2/){split($4,a,";");pirna[substr($4,1,25)]+=a[2]}(FILENAME!~/realbed2/){print $0"\t"pirna[substr($4,1,25)]/3}' *p7*.2555.new.realbed2 12.inter.2555.1ppm.realbed.piRNAs > 12.inter.2555.1ppm.realbed.piRNAs.p7
awk 'BEGIN{FS=OFS="\t"}{split($5,a,";");total=0;for(i=2;i<=13;i++){total+=a[i]};print $0"\t"total/12}' 12.inter.2555.1ppm.realbed.piRNAs.p7 > 12.inter.2555.1ppm.realbed.piRNAs.p7_WT
#OUTPUT by COLUMN
 1-chr piRNA
 2-start piRNA
 3-end piRNA
!4-piRNA sequence
 5-12reps of wt piRNA abundance for each rep separately
 6-piRNA strand
 7-piRNA cluster
 8-piRNA gene
 9-piRNA RMSK
 10-piRNA RMSK strand
 11-mean piRNA abundance in 9reps of triple pi2pi9pi17 mutants
 12-mean piRNA abundance in 3reps of pi7 mutants
 13-mean piRNA abundance in 12reps of wt

#3.find all possible pairs of piRNAs and 5'P ends with at least 12 complementary nucleotides
#3a.split piRNAs into 30 piRNA chunks form arallel computing
for i in 12.inter.2555.1ppm.realbed.piRNAs.p7_WT; do echo $i; split -a 4 -d -l 30 $i $i.xx. ; done 

#3b.find pairs
for j in $(find . -maxdepth 1 -name "all_deg_sites_01ppm_from_1tr.info"); do for i in $(find . -maxdepth 1 -name "*.xx.*"); do si=${i:2}; si=${si//.xx./.xy.};  /project/umw_phil_zamore/common/pipelines/ms_1G_1core_24h "awk -v n=$n 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){pirnan++;pirna[pirnan]=\$4\"\\t\"\$11\"\\t\"\$12\"\\t\"\$13;pirnaseq[pirnan]=\$4}(FNR<NR){for (i=1;i<=pirnan;i++){num=0;j=2;if(substr(pirnaseq[i],j,1)==substr(\$17,5+j,1)){num++;pattern=\"1\"}else{pattern=\"0\"};for (j=3;j<=25;j++){if(substr(pirnaseq[i],j,1)==substr(\$17,5+j,1)){num++;pattern=pattern\",1\"}else{pattern=pattern\",0\"}};if(num>=12){print \$1\":\"\$2\":\"\$6\"\\t\"\$4\"\\t\"pirna[i]\"\\t\"pattern\"\\t\"num}}}' $i $j > $j.$si.alltoall"; done; done

#3c.merge data after parallel computing
cat *.xy.* > all_deg_sites_01ppm_from_1tr.info.12.inter.2555.1ppm.realbed.piRNAs.p7_WT.merged.alltoall

#OUTPUT by COLUMN (1-2:5'P site, 3-6:piRNA, 7,8:pairing pattern)
1-chr:start:end
2-5P site_log2_fc mut/wt
3-piRNA sequence
4-mean piRNA abundance in 9reps of triple pi2pi9pi17 mutants
5-mean piRNA abundance in 3reps of pi7 mutants
6-mean piRNA abundance in 12reps of wt
7-matching pattern between g2g25 (0-unpaired,1-paired)
8-total number of matches between g2g25

#4.make 32 permutations of 4 WT and 8 p2p9p17 and p7 datasets, use ".unique.out.m.5.d" files from 5P-seq_1_initial_processing.sh
for i in $(find . -maxdepth 1 -name "*.d"); do /project/umw_phil_zamore/common/pipelines/ms_64G_1days "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){split (\$5,a,\",\");if(a[3]>0.1){data[\$1\":\"\$2\":\"\$6]=\$4}}(FNR<NR)&&(data[\$1]!=\"\"){print \$1\"\\t\"data[\$1]\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10}' $i all_deg_sites_01ppm_from_1tr.info.12.inter.2555.1ppm.realbed.piRNAs.p7_WT.merged.alltoall > $i.all"; done

#OUTPUT by COLUMN (1-2:5'P site, 3-6:piRNA, 7,8:pairing pattern)
1-chr:start:end
2-5P site_log2_fc mut/wt
3-piRNA sequence
4-mean piRNA abundance in 9reps of triple pi2pi9pi17 mutants
5-mean piRNA abundance in 3reps of pi7 mutants
6-mean piRNA abundance in 12reps of wt
7-matching pattern between g2g25 (0-unpaired,1-paired)
8-total number of matches between g2g25


#5.keep only unique 5'P site-piRNA pairs with 12..20 match thresholds
for i in $(find . -maxdepth 1 -name "*p2p9p17*.all"); do for n in {12..20}; do /project/umw_phil_zamore/common/pipelines/ms_16G_1core_24h "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(\$8>=$n)&&(\$4<0.1){count[\$1]++}(FNR<NR)&&(\$8>=$n)&&(\$4<0.1)&&(count[\$1]==1){string=\$1\"\\t\"\$3\"\\t\"\$6\"\\t\"\$8;split(\$7,a,\",\");for(i=1;i<=24;i++){string=string\"\\t\"a[i]}; if(\$2<=-3){string=string\"\\t1\"}else{string=string\"\\t0\"};print string}' $i $i > $i.01ppm.$n.dat"; done; done
for i in $(find . -maxdepth 1 -name "*p7*.all"); do for n in {12..20}; do /project/umw_phil_zamore/common/pipelines/ms_16G_1core_24h "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(\$8>=$n)&&(\$5<0.1){count[\$1]++}(FNR<NR)&&(\$8>=$n)&&(\$5<0.1)&&(count[\$1]==1){string=\$1\"\\t\"\$3\"\\t\"\$6\"\\t\"\$8;split(\$7,a,\",\");for(i=1;i<=24;i++){string=string\"\\t\"a[i]}; if(\$2<=-3){string=string\"\\t1\"}else{string=string\"\\t0\"};print string}' $i $i > $i.01ppm.$n.dat"; done; done

#OUTPUT by COLUMN (1:5'P site, 2-4:piRNA, 5-29:pairing pattern)
1-5'P site _site coord
2-piRNA sequence
3-mean piRNA abundance in 12reps of wt
4-total number of matches between g2g25
5-28-matching pattern between g2-g25
29-cleaved / not cleaved; target site is cleaved if 5'P site read abundance log2(mut/wt) <= -3, otherwise - uncleaved

#6.get dG for piRNA:site pairing with RNAplex and info about target site location
module load ViennaRNA/2.1.6h
for i in $(ls *.dat); do /project/umw_phil_zamore/common/pipelines/ms_short_1G "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR){deg_pirna[\$1]=\$2}(FNR<NR)&&(deg_pirna[\$1\":\"\$2\":\"\$6]!=\"\"){print \">\"\$1\":\"\$2\":\"\$6\":\"deg_pirna[\$1\":\"\$2\":\"\$6]\":\"\$13\":\"\$14\":\"\$15\":\"\$17\":\"\$8\":\"\$9; guide = substr(deg_pirna[\$1\":\"\$2\":\"\$6],2,24); gsub (\"T\",\"U\",guide); print guide; target=\"\";rev_target=substr(\$17,7,24);for (j=length(rev_target);j>=1;j--){if(substr(rev_target,j,1)==\"A\"){target=target\"U\"};if(substr(rev_target,j,1)==\"T\"){target=target\"A\"};if(substr(rev_target,j,1)==\"G\"){target=target\"C\"};if(substr(rev_target,j,1)==\"C\"){target=target\"G\"}}; print target }' $i all_deg_sites_01ppm_from_1tr.info | RNAplex -T 33 | awk 'BEGIN{FS=OFS=\"\\t\"}(FNR%3==1){split(substr(\$1,2),dat,\":\");getline;n=split(\$1,en,\"(\");print dat[1]\":\"dat[2]\":\"dat[3]\"\\t\"dat[4]\"\\t\"dat[5]\"\\t\"dat[6]\"\\t\"dat[7]\"\\t\"dat[8]\"\\t\"substr(en[n],1,length(en[n])-1)}' > $i.moreinfo"; done

#OUTPUT by COLUMN
1-5P site _chr:deg_start:deg_strand
2-piRNA sequence
3-ENSMUST name
4-loc_position in ENSMUST
5-utr_type
6-revcomp transcriptseq/or/genomeseq for -25+15, 
7-dG with RNAplex


#7.use "dat" files from step 5 and "moreinfo" files from step 6 to assemble tables for ML
for i in $(ls *.dat); do /project/umw_phil_zamore/common/pipelines/ms_short_32G "awk 'BEGIN{FS=OFS=\"\\t\"}
(index(FILENAME,\"moreinfo\")>0){utr[\$1\":\"\$2]=\$5;revcomp[\$1\":\"\$2]=\$6;dG[\$1\":\"\$2]=-1*\$7}
 (index(FILENAME,\"dat\")>0)&&(index(FILENAME,\"mrnf\")==0)&&(index(revcomp[\$1\":\"\$2],\"N\")==0){
 if(substr(revcomp[\$1\":\"\$2],6,1)==\"T\"){t1A=1}else{t1A=0};
 if(substr(revcomp[\$1\":\"\$2],6,1)==\"A\"){t1U=1}else{t1U=0};
 if(substr(revcomp[\$1\":\"\$2],6,1)==\"G\"){t1C=1}else{t1C=0};
 if(substr(revcomp[\$1\":\"\$2],6,1)==\"C\"){t1G=1}else{t1G=0};
 if(utr[\$1\":\"\$2]==\"5utr\"){utr5=1}else{utr5=0};
 if(utr[\$1\":\"\$2]==\"cds\"){cds=1}else{cds=0};
 if(utr[\$1\":\"\$2]==\"3utr\"){utr3=1}else{utr3=0};
 if(utr[\$1\":\"\$2]==\"nc\"){nc=1}else{nc=0};
 print \$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11\"\\t\"\$12\"\\t\"\$13\"\\t\"\$14\"\\t\"\$15\"\\t\"\$16\"\\t\"\$17\"\\t\"\$18\"\\t\"\$19\"\\t\"\$20\"\\t\"\$21\"\\t\"\$22\"\\t\"\$23\"\\t\"\$24\"\\t\"\$25\"\\t\"\$26\"\\t\"\$27\"\\t\"\$28\"\\t\"t1A\"\\t\"t1U\"\\t\"t1C\"\\t\"t1G\"\\t\"utr5\"\\t\"cds\"\\t\"utr3\"\\t\"nc\"\\t\"dG[\$1\":\"\$2]\"\\t\"\$29}' $i.moreinfo $i > $i.36col"; done

1   - mean piRNA abundance in 12reps of wt
2   - total number of matches between g2g25
3-26- pairing pattern between g2g25
27  - t1A
28  - t1U
29  - t1C
30  - t1G
31  - 5utr
32  - cds
33  - 3utr
34  - ncRNA
35  - dG for piRNA:site pairing with RNAplex
36  - cleaved / not cleaved; target site is cleaved if 5P site read abundance log2(mut/wt) <= -3, otherwise - uncleaved

#7.run 5-fold 5x repeated cross-validated logistic regression classifier fits
for i in $(find . -maxdepth 1 -name "*p2p9p17*.dat.36col"); do /project/umw_phil_zamore/common/pipelines/ms_16G_1core_24h "python3 ./5P-seq_5_logistic_regression_classifier.py $i > $i.lg35"; done
