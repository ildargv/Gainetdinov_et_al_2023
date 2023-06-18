#1.merge reads with the same prefix.2555 (output is longest sequence followed by nt18-nt55,t18-t55,untailed,tailed,total)
for i in $(find . -maxdepth 1 -name "*.p20.bed2.rpm"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "sort -k3,3 $i > $i.sorted; gawk 'BEGIN{ch=1;started=0}(length (\$3)>24)&&(length (\$3)<56)&&(started==0){for (i=18;i<=55;i++){tempnt[i]=0;
tempt[i]=0};prefix=\$3;started=1;tempnt[length(prefix)]+=\$4;
tempt[length(prefix)]+=\$5;}
(length (\$3)>24)&&(length (\$3)<56)&&(started==1){while ((length (\$3)<56)&&(length (\$3)>24)&&(substr(\$3,1,length(prefix))==prefix)&&(ch!=0)){prefix=\$3;tempnt[length(prefix)]+=\$4;tempt[length(prefix)]+=\$5;ch=getline};
if((length(prefix)>24)&&(length(prefix)<56)){for (i=18;i<=55;i++){nt[prefix,i]=tempnt[i];
t[prefix,i]=tempt[i]};
all[prefix]=1};
for (i=18;
i<=55;
i++){tempnt[i]=0;
tempt[i]=0};
tempnt[length(\$3)]+=\$4;
tempt[length(\$3)]+=\$5;
prefix=\$3}END{for (r in all){totalnt=0;totalt=0;printf (r);
for (i=18;
i<=55;
i++){printf (\"\\t\");
printf ((nt[r,i]-t[r,i])); totalnt+=(nt[r,i]-t[r,i])};
for (i=18;
i<=55;
i++){printf (\"\\t\");
printf (t[r,i]); totalt+=t[r,i]};
printf (\"\\t\");
printf (totalnt);
printf (\"\\t\");
printf (totalt);
printf (\"\\t\");
printf (totalt+totalnt);
printf (\"\\n\")}}' $i.sorted > $i.2555.new; rm $i.sorted"; done

#2.get bed2 files with all piRNAs (1ppm+ corresponds to 10 molecules per cell) requiring g1g25 to be uniquely mapped
module load bowtie/1.0.0
for i in $(find . -name "*.bed2.rpm.2555.new"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "awk '(length(\$1)>23)&&(\$80>0){max1=0; max1n=0; for (i=8;i<=39;i++){if ((\$i+\$(i+38))>max1){max1=(\$i+\$(i+38)); max1n=i+16}};len[\$1]=max1n; ab[\$1]=\$80; 
ab7[substr(\$1,2,7)]+=\$80; 
ab8[substr(\$1,2,8)]+=\$80; 
ab9[substr(\$1,2,9)]+=\$80; 
ab10[substr(\$1,2,10)]+=\$80; 
ab11[substr(\$1,2,11)]+=\$80; 
ab12[substr(\$1,2,12)]+=\$80; 
ab13[substr(\$1,2,13)]+=\$80; 
ab14[substr(\$1,2,14)]+=\$80; 
ab15[substr(\$1,2,15)]+=\$80; 
ab16[substr(\$1,2,16)]+=\$80; 
ab17[substr(\$1,2,17)]+=\$80; 
ab18[substr(\$1,2,18)]+=\$80; 
ab19[substr(\$1,2,19)]+=\$80; 
ab20[substr(\$1,2,20)]+=\$80; 
ab21[substr(\$1,2,21)]+=\$80; 
ab22[substr(\$1,2,22)]+=\$80; 
ab23[substr(\$1,2,23)]+=\$80; 
ab24[substr(\$1,2,24)]+=\$80; 
ab31[substr(\$1,2,31)]+=\$80; 
num7[substr(\$1,2,7)]++; 
num8[substr(\$1,2,8)]++; 
num9[substr(\$1,2,9)]++;
num10[substr(\$1,2,10)]++; 
num11[substr(\$1,2,11)]++; 
num12[substr(\$1,2,12)]++; 
num13[substr(\$1,2,13)]++; 
num14[substr(\$1,2,14)]++; 
num15[substr(\$1,2,15)]++; 
num16[substr(\$1,2,16)]++; 
num17[substr(\$1,2,17)]++; 
num18[substr(\$1,2,18)]++; 
num19[substr(\$1,2,19)]++; 
num20[substr(\$1,2,20)]++; 
num21[substr(\$1,2,21)]++; 
num22[substr(\$1,2,22)]++; 
num23[substr(\$1,2,23)]++; 
num24[substr(\$1,2,24)]++; 
num31[substr(\$1,2,31)]++}END{for (i in ab){print \">\"i\";\"ab[i]\";\"len[i]\";\"ab7[substr(i,2,7)]\";\"ab8[substr(i,2,8)]\";\"ab9[substr(i,2,9)]\";\"ab10[substr(i,2,10)]\";\"ab11[substr(i,2,11)]\";\"ab12[substr(i,2,12)]\";\"ab13[substr(i,2,13)]\";\"ab14[substr(i,2,14)]\";\"ab15[substr(i,2,15)]\";\"ab16[substr(i,2,16)]\";\"ab17[substr(i,2,17)]\";\"ab18[substr(i,2,18)]\";\"ab19[substr(i,2,19)]\";\"ab20[substr(i,2,20)]\";\"ab21[substr(i,2,21)]\";\"ab22[substr(i,2,22)]\";\"ab23[substr(i,2,23)]\";\"ab24[substr(i,2,24)]\";\"ab31[substr(i,2,31)]\";\"num7[substr(i,2,7)]\";\"num8[substr(i,2,8)]\";\"num9[substr(i,2,9)]\";\"num10[substr(i,2,10)]\";\"num11[substr(i,2,11)]\";\"num12[substr(i,2,12)]\";\"num13[substr(i,2,13)]\";\"num14[substr(i,2,14)]\";\"num15[substr(i,2,15)]\";\"num16[substr(i,2,16)]\";\"num17[substr(i,2,17)]\";\"num18[substr(i,2,18)]\";\"num19[substr(i,2,19)]\";\"num20[substr(i,2,20)]\";\"num21[substr(i,2,21)]\";\"num22[substr(i,2,22)]\";\"num23[substr(i,2,23)]\";\"num24[substr(i,2,24)]\";\"num31[substr(i,2,31)]; print substr(i,1,25)}}' $i > $i.fa; bowtie -f -v 0 -a mm10 $i.fa | awk '{print \$3\"\\t\"\$4\"\\t\"(\$4+length(\$5))\"\\t\"\$1\"\\t\"\$7+1\"\\t\"\$2}' > $i.prebed2;  awk '(FNR==NR){a[\$4]++}(FNR<NR)&&(a[\$4]==1){print}' $i.prebed2 $i.prebed2 > $i.realbed2; rm $i.fa"; done

#OUTPUT by COLUMN
1 sequence
2 abundance
3 peak nt
4-22 total abundance g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21,g22,g23,g24,g31,
23-41 redundancy g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21,g22,g23,g24,g31,
42 mapping redundancy of g1g25

#3.merge 12 WT replicates to make realbed of the intersection of the 12 WT animals
awk '{split ($4,a,";");coord[substr(a[1],1,25)]=$1"\t"$2"\t"$3;strand[substr(a[1],1,25)]=$6;seq[substr(a[1],1,25)]=substr(a[1],1,29);abund[substr(a[1],1,25)]=abund[substr(a[1],1,25)]";"a[2]; if(a[2]>1){above1[substr(a[1],1,25)]++}else{above1[substr(a[1],1,25)]+=0}}END{for (i in above1){if(above1[i]==12){print coord[i]"\t"seq[i]"\t"abund[i]"\t"strand[i]}}}' *WT*.2555.new.realbed2 > 12.inter.2555.1ppm.realbed

#4.add cluster, gene, TE info
#4a.add cluster info
module load bedtools/2.26.0
for i in 12.inter.2555.1ppm.realbed; do /project/umw_phil_zamore/common/pipelines/ms_short_8G "bedtools intersect -s -wao -a $i -b piRNA.cluster.pachytene.bed6 | awk '{a[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6]=substr(\$10,1,length(\$10)-2)}END{for (i in a){if(a[i]!=\"\"){print i\"\\t\"a[i]}else{print i\"\\t0\"}}}' > $i.cl; "; done

#4b.add gene info
module load bedtools/2.26.0
for i in 12.inter.2555.1ppm.realbed.cl; do /project/umw_phil_zamore/common/pipelines/ms_short_8G "bedtools intersect -s -wao -a $i -b mm10.38.92.gtf.exons.geneNnameNtype.bed | awk '{a[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7]=\$12}END{for (i in a){if(a[i]!=\"\"){print i\"\\t\"a[i]}else{print i\"\\t0\"}}}' > $i.g; "; done

#4c.add TE info
module load bedtools/2.26.0
for i in 12.inter.2555.1ppm.realbed.cl.g; do /project/umw_phil_zamore/common/pipelines/ms_short_8G "bedtools intersect -wao -a $i -b mm10.rmsk.bed | awk '{if(\$6==\$14){a[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8]=\$12\"\\t+\"}else{a[\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8]=\$12\"\\t-\"}}END{for (i in a){print i\"\\t\"a[i]}}' > $i.r; "; done

#5.add mean p2p9p17 triple mutant piRNA abundance from 9 reps
for i in 12.inter.2555.1ppm.realbed.cl.g.r; do awk '(FILENAME~/realbed2/){split($4,a,";");pirna[substr($4,1,25)]+=a[2]}(FILENAME!~/realbed2/){print $0"\t"pirna[substr($4,1,25)]/9}' *p2p9p17*.2555.new.realbed2 $i > $i.p2p9p17; done

#6.remove piRNAs from genomic strand that is opposite to the annotated clusters
module load bedtools/2.26.0
for i in 12.inter.2555.1ppm.realbed.cl.g.r.p2p9p17; do bedtools intersect -wao -a $i -b piRNA.cluster.pachytene.bed6 | awk 'BEGIN{FS=OFS="\t"}($NF==0){a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11]=1}($NF!=0)&&((($6=="+")&&($17=="+"))||(($6=="-")&&($17=="-"))){a[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11]=1}END{for (i in a){print i}}' > 12.inter.2555.1ppm.realbed.piRNAs; done

#OUTPUT by COLUMN
 1-chr piRNA
 2-start piRNA
 3-end piRNA
!4-piRNA sequence
 5-12reps of wt abundance for each
 6-piRNA strand
 7-cluster
 8-piRNA gene
 9-piRNA RMSK
 10-piRNA RMSK strand
 11-mean piRNA abundance in 9reps of triple pi2pi9pi17 mutants
