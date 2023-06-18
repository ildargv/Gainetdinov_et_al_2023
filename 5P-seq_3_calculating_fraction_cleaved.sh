#1.use "1nt/2nt/3nt" files from5P-seq_2_matching_with_piRNAs.sh and  "unique.out.m.5.d" files from 5P-seq_1_initial_processing.sh to get /mut2, mut3, /mut4 combinations to get to 16 total WT/mut combinations
#1a.contiguous complementarity from g2, g3, g4, g5, g6, etc.
for f in $(find . -maxdepth 1 -name "*.ms.?nt"); do sf=${f:2}; for d in $(find . -name ${sf%%.*}."*".d);do  sd=${d:2};/project/umw_phil_zamore/common/pipelines/ms_64G_1day "awk -v f=\"$sf\" -v d=\"$sd\" 'BEGIN{FS=OFS=\"\\t\";fn=split(f,fa,\".\");dn=split(d,da,\".\")}(FILENAME==\"12.inter.2555.1ppm.realbed.piRNAs\"){pirn[substr(\$4,2,fa[fn-2]-1)]++; split(\$5,pirarray,\";\");for(i=2;i<=13;i++){pirwt[substr(\$4,2,fa[fn-2]-1)]+=(pirarray[i]/12)};pirmut[substr(\$4,2,fa[fn-2]-1)]+=\$11}(FILENAME~/m.5.d/){split(\$5,degarray,\",\");deg4[\$1\":\"\$2\":\"\$6]=\$4;deg5[\$1\":\"\$2\":\"\$6]=degarray[1]\"\\t\"degarray[2]\"\\t\"degarray[3]\"\\t\"degarray[4]}(FILENAME~/merged/){\$30=\"na\"; \$4=deg4[\$1\":\"\$2\":\"\$6];\$5=\"na\";print \$0\"\\t\"deg5[\$1\":\"\$2\":\"\$6]\"\\t\" pirwt[substr(\$29,2,fa[fn-2]-1)]\"\\t\" pirmut[substr(\$29,2,fa[fn-2]-1)]\"\\t\" pirn[substr(\$29,2,fa[fn-2]-1)] > fa[1]\".\"da[2]\".\"da[4]\".\"fa[fn-2]\".\"fa[fn-1]\".\"fa[fn]\".fd\"}' 12.inter.2555.1ppm.realbed.piRNAs $d $f"; done; done

#1b.mononucleotide mismatch in contiguous complementarity from g2
for f in $(find . -maxdepth 1 -name "*.ms.*.mm.?nt"); do sf=${f:2}; for d in $(find . -name ${sf%%.*}."*".d);do  sd=${d:2};/project/umw_phil_zamore/common/pipelines/ms_64G_1day "awk -v f=\"$sf\" -v d=\"$sd\" 'BEGIN{FS=OFS=\"\\t\";fn=split(f,fa,\".\");dn=split(d,da,\".\")}(FILENAME==\"12.inter.2555.1ppm.realbed.piRNAs\"){pirn[substr(\$4,2,fa[fn-4]-1)]++; split(\$5,pirarray,\";\");for(i=2;i<=13;i++){pirwt[substr(\$4,2,fa[fn-4]-1)]+=(pirarray[i]/12)};pirmut[substr(\$4,2,fa[fn-4]-1)]+=\$11}(FILENAME~/m.5.d/){split(\$5,degarray,\",\");deg4[\$1\":\"\$2\":\"\$6]=\$4;deg5[\$1\":\"\$2\":\"\$6]=degarray[1]\"\\t\"degarray[2]\"\\t\"degarray[3]\"\\t\"degarray[4]}(FILENAME~/merged/){\$30=\"na\"; \$4=deg4[\$1\":\"\$2\":\"\$6];\$5=\"na\";print \$0\"\\t\"deg5[\$1\":\"\$2\":\"\$6]\"\\t\" pirwt[substr(\$29,2,fa[fn-4]-1)]\"\\t\" pirmut[substr(\$29,2,fa[fn-4]-1)]\"\\t\" pirn[substr(\$29,2,fa[fn-4]-1)] > fa[1]\".\"da[2]\".\"da[4]\".\"fa[fn-4]\".\"fa[fn-2]\".\"fa[fn-1]\".\"fa[fn]\".fd\"}' 12.inter.2555.1ppm.realbed.piRNAs $d $f"; done; done

#1b.mononucleotide bulge (target insertion or deletion) in contiguous complementarity from g2
for f in $(find . -maxdepth 1 -name "*.ms.*.?b.?nt"); do sf=${f:2}; for d in $(find . -name ${sf%%.*}."*".d);do  sd=${d:2};/project/umw_phil_zamore/common/pipelines/ms_64G_1day "awk -v f=\"$sf\" -v d=\"$sd\" 'BEGIN{FS=OFS=\"\\t\";fn=split(f,fa,\".\");dn=split(d,da,\".\")}(FILENAME==\"12.inter.2555.1ppm.realbed.piRNAs\"){pirn[substr(\$4,2,fa[fn-4]-1)]++; split(\$5,pirarray,\";\");for(i=2;i<=13;i++){pirwt[substr(\$4,2,fa[fn-4]-1)]+=(pirarray[i]/12)};pirmut[substr(\$4,2,fa[fn-4]-1)]+=\$11}(FILENAME~/m.5.d/){split(\$5,degarray,\",\");deg4[\$1\":\"\$2\":\"\$6]=\$4;deg5[\$1\":\"\$2\":\"\$6]=degarray[1]\"\\t\"degarray[2]\"\\t\"degarray[3]\"\\t\"degarray[4]}(FILENAME~/merged/){\$30=\"na\"; \$4=deg4[\$1\":\"\$2\":\"\$6];\$5=\"na\";print \$0\"\\t\"deg5[\$1\":\"\$2\":\"\$6]\"\\t\" pirwt[substr(\$29,2,fa[fn-4]-1)]\"\\t\" pirmut[substr(\$29,2,fa[fn-4]-1)]\"\\t\" pirn[substr(\$29,2,fa[fn-4]-1)] > fa[1]\".\"da[2]\".\"da[4]\".\"fa[fn-4]\".\"fa[fn-2]\".\"fa[fn-1]\".\"fa[fn]\".fd\"}' 12.inter.2555.1ppm.realbed.piRNAs $d $f"; done; done


#OUTPUT by COLUMN (1-25: 5' P site, 26-36: piRNA)
 1-5 P site_chr
 2-5 P site_start
 3-5 P site_end
 4-5 P site_log2change
 5-na
 6-5 P site_strand
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
 30-na
 31-piRNA strand
 32-cluster
 33-piRNA gene
 34-piRNA RMSK
 35-piRNA RMSK strand
 36-piRNA 9reps of triple data abundance
 37-wt_count
 38-mut_count
 39-wt_ppm
 40-mut_ppm
 41-wt cumulative abundance of all piRNAs with given pairing pattern
 42-mut cumulative abundance of all piRNAs with given pairing pattern
 43-mut cumulative number of all piRNAs with given pairing pattern


#2.Calculate fraction cleaved fractions of 5P-seq reads decreased 2-,4-,8-,inf-fold in mut compared to wt for REMOVED piRNAs(<0.2ppm in p2p9p17) minus CONTROL/NOT REMOVED piRNAs (same or higher piRNA abundance in mut compared to wt)
#2a.for all piRNAs
for i in $(find . -maxdepth 1 -name "*nt.fd"); do /project/umw_phil_zamore/common/pipelines/ms_short_32G " awk 'BEGIN{FS=OFS=\"\\t\"}(\$39>=0.05)&&(\$42<0.1){total[\$1\".\"\$2\".\"\$6]=1}(\$39>=0.05)&&(\$42<0.1)&&(\$4<=-1){fr1[\$1\".\"\$2\".\"\$6]=1}(\$39>=0.05)&&(\$42<0.1)&&(\$4<=-2){fr2[\$1\".\"\$2\".\"\$6]=1}(\$39>=0.05)&&(\$42<0.1)&&(\$4<=-3){fr3[\$1\".\"\$2\".\"\$6]=1}(\$39>=0.05)&&(\$42<0.1)&&(\$4<=-15){fr0[\$1\".\"\$2\".\"\$6]=1}(\$39>=0.05)&&(\$42>=\$41){totalC[\$1\".\"\$2\".\"\$6]=1}(\$39>=0.05)&&(\$42>=\$41)&&(\$4<=-1){fr1C[\$1\".\"\$2\".\"\$6]=1}(\$39>=0.05)&&(\$42>=\$41)&&(\$4<=-2){fr2C[\$1\".\"\$2\".\"\$6]=1}(\$39>=0.05)&&(\$42>=\$41)&&(\$4<=-3){fr3C[\$1\".\"\$2\".\"\$6]=1}(\$39>=0.05)&&(\$42>=\$41)&&(\$4<=-15){fr0C[\$1\".\"\$2\".\"\$6]=1}END{for(i in total){totalt++};for(i in fr1){fr1t++};for(i in fr2){fr2t++};for(i in fr3){fr3t++};for(i in fr0){fr0t++};for(i in totalC){totaltC++};for(i in fr1C){fr1tC++};for(i in fr2C){fr2tC++};for(i in fr3C){fr3tC++};for(i in fr0C){fr0tC++};if((totalt!=0)&&(totaltC!=0)){print FILENAME\"\t\"((fr1t/totalt)-(fr1tC/totaltC))\"\t\"((fr2t/totalt)-(fr2tC/totaltC))\"\t\"((fr3t/totalt)-(fr3tC/totaltC))\"\t\"((fr0t/totalt)-(fr0tC/totaltC))\"\t\"totalt\"\t\"totaltC}else{totalt=0;print FILENAME\"\t\"totalt\"\t\"totalt\"\t\"totalt\"\t\"totalt\"\t\"totalt\"\t\"totalt}}' $i > $i.005.rNr"; done

#2b.binned by wt piRNA abundance: columns: all piRNAs, 0-3ppm, 3-5ppm, 5-10ppm, 10-50ppm, >50ppm; 1ppm ~ 10 piRNAs per primary spermatocyte (~10pM intracellular concentration) 
for i in $(find . -maxdepth 1 -name "*nt.fd"); do /project/umw_phil_zamore/common/pipelines/ms_short_32G " awk 'BEGIN{FS=OFS=\"\\t\"}(\$39>=0.05)&&(\$42<0.1){total[\$1\".\"\$2\".\"\$6]=1;ttype[\$1\".\"\$2\".\"\$6]=\$41}(\$39>=0.05)&&(\$42<0.1)&&(\$4<=-3){fr[\$1\".\"\$2\".\"\$6]=1;frtype[\$1\".\"\$2\".\"\$6]=\$41}(\$39>=0.05)&&(\$42>=\$41){totalC[\$1\".\"\$2\".\"\$6]=1;ttypeC[\$1\".\"\$2\".\"\$6]=\$41}(\$39>=0.05)&&(\$42>=\$41)&&(\$4<=-3){frC[\$1\".\"\$2\".\"\$6]=1;frtypeC[\$1\".\"\$2\".\"\$6]=\$41}END{for(i in total){t++; 
if(ttype[i]>50){t0++};
if((ttype[i]>10)&&(ttype[i]<=50)){t1++};
if((ttype[i]>5)&&(ttype[i]<=10)){t2++};
if((ttype[i]>3)&&(ttype[i]<=5)){t3++};
if(ttype[i]<=3){t4++}};

for(i in totalC){tC++; 
if(ttypeC[i]>50){t0C++};
if((ttypeC[i]>10)&&(ttypeC[i]<=50)){t1C++};
if((ttypeC[i]>5)&&(ttypeC[i]<=10)){t2C++};
if((ttypeC[i]>3)&&(ttypeC[i]<=5)){t3C++};
if(ttypeC[i]<=3){t4C++}};


for(i in fr){frt++;
if(frtype[i]>50){fr0++};
if((frtype[i]>10)&&(frtype[i]<=50)){fr1++};
if((frtype[i]>5)&&(frtype[i]<=10)){fr2++};
if((frtype[i]>3)&&(frtype[i]<=5)){fr3++};
if(frtype[i]<=3){fr4++}};

for(i in frC){frtC++;
if(frtypeC[i]>50){fr0C++};
if((frtypeC[i]>10)&&(frtypeC[i]<=50)){fr1C++};
if((frtypeC[i]>5)&&(frtypeC[i]<=10)){fr2C++};
if((frtypeC[i]>3)&&(frtypeC[i]<=5)){fr3C++};
if(frtypeC[i]<=3){fr4C++}};

if((t==0)||(tC==0)){ratt=0}else{ratt=((frt/t)-(frtC/tC))};
if((t4==0)||(t4C==0)){ratt4=0}else{ratt4=((fr4/t4)-(fr4C/t4C))};
if((t3==0)||(t3C==0)){ratt3=0}else{ratt3=((fr3/t3)-(fr3C/t3C))};
if((t2==0)||(t2C==0)){ratt2=0}else{ratt2=((fr2/t2)-(fr2C/t2C))};
if((t1==0)||(t1C==0)){ratt1=0}else{ratt1=((fr1/t1)-(fr1C/t1C))};
if((t0==0)||(t0C==0)){ratt0=0}else{ratt0=((fr0/t0)-(fr0C/t0C))};

print FILENAME\"\t\"ratt\"\\t\"ratt4\"\\t\"ratt3\"\\t\"ratt2\"\\t\"ratt1\"\\t\"ratt0}' $i > $i.005.ab.rNr"; done






















