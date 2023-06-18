!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for i in $(find . -name "*.fastq"); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core "split -a 4 -d -l 4000000 $i $i.xx." ; done 

#let-7a: get site types for Kd gN_Kmer matches: parallel computing
for f in $(ls *let-7a_*.fastq.xx.????); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(FNR%2==0){flt[FNR]=\$1}(FNR%4==2)&&(FNR<NR){if((length(\$1)==40)&&(index(\$1,\"N\")==0)&&(substr(\$1,38,3)==\"TGT\")){full=substr(\$1,1,37);if((index(flt[2],full)==0)&&(index(flt[4],full)==0)&&(index(full,flt[6])==0)&&(index(full,flt[8])==0)){sites=0;site_loc=\"\";
for(k=6;k<=11;k++){for (g=2;g<=(23-k);g++){num=split( full, temp, substr(\"AACTATACAACCTACTACCTC\",(24-k)-g,k)); if(num>2){sites=2}; if((num==2)&&(substr(temp[1],length(temp[1]),1)!=substr(\"NAACTATACAACCTACTACCTCN\",(24-k)-g,1))&&(substr(temp[2],1,1)!=substr(\"NAACTATACAACCTACTACCTCN\",25-g,1)))
{site_loc=\"g\"g\"-\"(g+k-1);sites++}}};if(sites==1){counts[site_loc]++};if(sites==0){counts[\"nosite\"]++}}}}END{for (site in counts){print site\"\\t\"counts[site]}}' phiX_standards.fa $f > $f.let_gN_Kmer.types"; done

#miR-1: get site types for Kd gN_Kmer matches: parallel computing
for f in $(ls *miR-1_*.fastq.xx.????); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(FNR%2==0){flt[FNR]=\$1}(FNR%4==2)&&(FNR<NR){if((length(\$1)==40)&&(index(\$1,\"N\")==0)&&(substr(\$1,38,3)==\"TGT\")){full=substr(\$1,1,37);if((index(flt[2],full)==0)&&(index(flt[4],full)==0)&&(index(full,flt[6])==0)&&(index(full,flt[8])==0)){sites=0;site_loc=\"\";
for(k=6;k<=11;k++){for (g=2;g<=(23-k);g++){num=split( full, temp, substr(\"ATACATACTTCTTTACATTCC\",(24-k)-g,k)); if(num>2){sites=2}; if((num==2)&&(substr(temp[1],length(temp[1]),1)!=substr(\"NATACATACTTCTTTACATTCCN\",(24-k)-g,1))&&(substr(temp[2],1,1)!=substr(\"NATACATACTTCTTTACATTCCN\",25-g,1)))
{site_loc=\"g\"g\"-\"(g+k-1);sites++}}};if(sites==1){counts[site_loc]++};if(sites==0){counts[\"nosite\"]++}}}}END{for (site in counts){print site\"\\t\"counts[site]}}' phiX_standards.fa $f > $f.1_gN_Kmer.types"; done

#miR-124: get site types for Kd gN_Kmer matches: parallel computing
for f in $(ls *miR-124_*.fastq.xx.????); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(FNR%2==0){flt[FNR]=\$1}(FNR%4==2)&&(FNR<NR){if((length(\$1)==40)&&(index(\$1,\"N\")==0)&&(substr(\$1,38,3)==\"TGT\")){full=substr(\$1,1,37);if((index(flt[2],full)==0)&&(index(flt[4],full)==0)&&(index(full,flt[6])==0)&&(index(full,flt[8])==0)){sites=0;site_loc=\"\";
for(k=6;k<=11;k++){for (g=2;g<=(23-k);g++){num=split( full, temp, substr(\"TTGGCATTCACCGCGTGCCTT\",(24-k)-g,k)); if(num>2){sites=2}; if((num==2)&&(substr(temp[1],length(temp[1]),1)!=substr(\"NTTGGCATTCACCGCGTGCCTTN\",(24-k)-g,1))&&(substr(temp[2],1,1)!=substr(\"NTTGGCATTCACCGCGTGCCTTN\",25-g,1)))
{site_loc=\"g\"g\"-\"(g+k-1);sites++}}};if(sites==1){counts[site_loc]++};if(sites==0){counts[\"nosite\"]++}}}}END{for (site in counts){print site\"\\t\"counts[site]}}' phiX_standards.fa $f > $f.124_gN_Kmer.types"; done

#lsy-6: get site types for Kd gN_Kmer matches: parallel computing
for f in $(ls *lsy-6_*.fastq.xx.????); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(FNR%2==0){flt[FNR]=\$1}(FNR%4==2)&&(FNR<NR){if((length(\$1)==40)&&(index(\$1,\"N\")==0)&&(substr(\$1,38,3)==\"TGT\")){full=substr(\$1,1,37);if((index(flt[2],full)==0)&&(index(flt[4],full)==0)&&(index(full,flt[6])==0)&&(index(full,flt[8])==0)){sites=0;site_loc=\"\";
for(k=6;k<=11;k++){for (g=2;g<=(23-k);g++){num=split( full, temp, substr(\"TCGAAATGCGTCTCATACAAA\",(24-k)-g,k)); if(num>2){sites=2}; if((num==2)&&(substr(temp[1],length(temp[1]),1)!=substr(\"NTCGAAATGCGTCTCATACAAAN\",(24-k)-g,1))&&(substr(temp[2],1,1)!=substr(\"NTCGAAATGCGTCTCATACAAAN\",25-g,1)))
{site_loc=\"g\"g\"-\"(g+k-1);sites++}}};if(sites==1){counts[site_loc]++};if(sites==0){counts[\"nosite\"]++}}}}END{for (site in counts){print site\"\\t\"counts[site]}}' phiX_standards.fa $f > $f.lsy_gN_Kmer.types"; done

#miR-155: get site types for Kd gN_Kmer matches: parallel computing
for f in $(ls *miR-155_*.fastq.xx.????); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(FNR%2==0){flt[FNR]=\$1}(FNR%4==2)&&(FNR<NR){if((length(\$1)==40)&&(index(\$1,\"N\")==0)&&(substr(\$1,38,3)==\"TGT\")){full=substr(\$1,1,37);if((index(flt[2],full)==0)&&(index(flt[4],full)==0)&&(index(full,flt[6])==0)&&(index(full,flt[8])==0)){sites=0;site_loc=\"\";
for(k=6;k<=11;k++){for (g=2;g<=(24-k);g++){num=split( full, temp, substr(\"ACCCCTATCACGATTAGCATTA\",(25-k)-g,k)); if(num>2){sites=2}; if((num==2)&&(substr(temp[1],length(temp[1]),1)!=substr(\"NACCCCTATCACGATTAGCATTAN\",(25-k)-g,1))&&(substr(temp[2],1,1)!=substr(\"NACCCCTATCACGATTAGCATTAN\",26-g,1)))
{site_loc=\"g\"g\"-\"(g+k-1);sites++}}};if(sites==1){counts[site_loc]++};if(sites==0){counts[\"nosite\"]++}}}}END{for (site in counts){print site\"\\t\"counts[site]}}' phiX_standards.fa $f > $f.155_gN_Kmer.types"; done

#miR-7: get site types for Kd gN_Kmer matches: parallel computing
for f in $(ls *miR-7_*.fastq.xx.????); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(FNR%2==0){flt[FNR]=\$1}(FNR%4==2)&&(FNR<NR){if((length(\$1)==40)&&(index(\$1,\"N\")==0)&&(substr(\$1,38,3)==\"TGT\")){full=substr(\$1,1,37);if((index(flt[2],full)==0)&&(index(flt[4],full)==0)&&(index(full,flt[6])==0)&&(index(full,flt[8])==0)){sites=0;site_loc=\"\";
for(k=6;k<=11;k++){for (g=2;g<=(24-k);g++){num=split( full, temp, substr(\"ACAACAAAATCACTAGTCTTCC\",(25-k)-g,k)); if(num>2){sites=2}; if((num==2)&&(substr(temp[1],length(temp[1]),1)!=substr(\"NACAACAAAATCACTAGTCTTCCN\",(25-k)-g,1))&&(substr(temp[2],1,1)!=substr(\"NACAACAAAATCACTAGTCTTCCN\",26-g,1)))
{site_loc=\"g\"g\"-\"(g+k-1);sites++}}};if(sites==1){counts[site_loc]++};if(sites==0){counts[\"nosite\"]++}}}}END{for (site in counts){print site\"\\t\"counts[site]}}' phiX_standards.fa $f > $f.7_gN_Kmer.types"; done

#merge parallels
for f in $(ls *.fastq.xx.0000.*_gN_Kmer.types); do echo $f; first=${f%.xx.0000.*}; last=${f#*.xx.0000.}; awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' $first.xx.????.$last > $first.merged.$last; done


#merge runs
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' let-7a_miR-155_miR-124_lsy-6_purification1_inputlibrary2.fastq.merged.let_gN_Kmer.types let-7a_miR-155_purification1_inputlibrary.fastq.merged.let_gN_Kmer.types   > let-7a_purification1_inputlibraryM.fastq.merged.let_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' let-7a_purification1_0.4percent.fastq.merged.let_gN_Kmer.types    > let-7a_purification1_0.4percentM.fastq.merged.let_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' let-7a_purification1_0percent.fastq.merged.let_gN_Kmer.types    > let-7a_purification1_0percentM.fastq.merged.let_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' let-7a_purification1_1.265percent.fastq.merged.let_gN_Kmer.types    > let-7a_purification1_1.265percentM.fastq.merged.let_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' let-7a_purification1_12.65percent.fastq.merged.let_gN_Kmer.types    > let-7a_purification1_12.65percentM.fastq.merged.let_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' let-7a_purification1_40percent.fastq.merged.let_gN_Kmer.types    > let-7a_purification1_40percentM.fastq.merged.let_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' let-7a_purification1_4percent.fastq.merged.let_gN_Kmer.types    > let-7a_purification1_4percentM.fastq.merged.let_gN_Kmer.types

awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' let-7a_miR-155_miR-124_lsy-6_purification1_inputlibrary2.fastq.merged.lsy_gN_Kmer.types miR-124_lsy-6_purification1_inputlibrary.fastq.merged.lsy_gN_Kmer.types   > lsy-6_purification1_inputlibraryM.fastq.merged.lsy_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' lsy-6_purification1_0.4percent.fastq.merged.lsy_gN_Kmer.types    > lsy-6_purification1_0.4percentM.fastq.merged.lsy_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' lsy-6_purification1_1.265percent.fastq.merged.lsy_gN_Kmer.types    > lsy-6_purification1_1.265percentM.fastq.merged.lsy_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' lsy-6_purification1_12.65percent.fastq.merged.lsy_gN_Kmer.types    > lsy-6_purification1_12.65percentM.fastq.merged.lsy_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' lsy-6_purification1_40percent.fastq.merged.lsy_gN_Kmer.types    > lsy-6_purification1_40percentM.fastq.merged.lsy_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' lsy-6_purification1_4percent.fastq.merged.lsy_gN_Kmer.types    > lsy-6_purification1_4percentM.fastq.merged.lsy_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-124_lsy-6_purification1_0percent.fastq.merged.lsy_gN_Kmer.types    > lsy-6_purification1_0percentM.fastq.merged.lsy_gN_Kmer.types

awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-124_purification2_0.4percent_rep1.fastq.merged.124_gN_Kmer.types miR-124_purification2_0.4percent_rep2.fastq.merged.124_gN_Kmer.types   > miR-124_purification2_0.4percentM.fastq.merged.124_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-124_purification2_0percent_rep1.fastq.merged.124_gN_Kmer.types miR-124_purification2_0percent_rep2.fastq.merged.124_gN_Kmer.types   > miR-124_purification2_0percentM.fastq.merged.124_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-124_purification2_1.265percent_rep1.fastq.merged.124_gN_Kmer.types miR-124_purification2_1.265percent_rep2.fastq.merged.124_gN_Kmer.types   > miR-124_purification2_1.265percentM.fastq.merged.124_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-124_purification2_12.65percent_rep1.fastq.merged.124_gN_Kmer.types miR-124_purification2_12.65percent_rep2.fastq.merged.124_gN_Kmer.types   > miR-124_purification2_12.65percentM.fastq.merged.124_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-124_purification2_40percent_rep1.fastq.merged.124_gN_Kmer.types miR-124_purification2_40percent_rep2.fastq.merged.124_gN_Kmer.types   > miR-124_purification2_40percentM.fastq.merged.124_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-124_purification2_4percent_rep1.fastq.merged.124_gN_Kmer.types miR-124_purification2_4percent_rep2.fastq.merged.124_gN_Kmer.types   > miR-124_purification2_4percentM.fastq.merged.124_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-124_purification2_inputlibrary_rep1.fastq.merged.124_gN_Kmer.types miR-124_purification2_inputlibrary_rep2.fastq.merged.124_gN_Kmer.types   > miR-124_purification2_inputlibraryM.fastq.merged.124_gN_Kmer.types

awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-1_purification3_0.4percent.fastq.merged.1_gN_Kmer.types    > miR-1_purification3_0.4percentM.fastq.merged.1_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-1_purification3_0percent.fastq.merged.1_gN_Kmer.types    > miR-1_purification3_0percentM.fastq.merged.1_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-1_purification3_1.265percent.fastq.merged.1_gN_Kmer.types    > miR-1_purification3_1.265percentM.fastq.merged.1_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-1_purification3_12.65percent.fastq.merged.1_gN_Kmer.types    > miR-1_purification3_12.65percentM.fastq.merged.1_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-1_purification3_40percent.fastq.merged.1_gN_Kmer.types    > miR-1_purification3_40percentM.fastq.merged.1_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-1_purification3_4percent.fastq.merged.1_gN_Kmer.types    > miR-1_purification3_4percentM.fastq.merged.1_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-1_purification3_inputlibrary.fastq.merged.1_gN_Kmer.types    > miR-1_purification3_inputlibraryM.fastq.merged.1_gN_Kmer.types

awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-7_purification2_0.4percent_rep1.fastq.merged.7_gN_Kmer.types miR-7_purification2_0.4percent_rep2.fastq.merged.7_gN_Kmer.types   > miR-7_purification2_0.4percentM.fastq.merged.7_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-7_purification2_0percent_rep1.fastq.merged.7_gN_Kmer.types miR-7_purification2_0percent_rep2.fastq.merged.7_gN_Kmer.types   > miR-7_purification2_0percentM.fastq.merged.7_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-7_purification2_1.265percent_rep1.fastq.merged.7_gN_Kmer.types miR-7_purification2_1.265percent_rep2.fastq.merged.7_gN_Kmer.types   > miR-7_purification2_1.265percentM.fastq.merged.7_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-7_purification2_12.65percent_rep1.fastq.merged.7_gN_Kmer.types miR-7_purification2_12.65percent_rep2.fastq.merged.7_gN_Kmer.types   > miR-7_purification2_12.65percentM.fastq.merged.7_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-7_purification2_40percent_rep1.fastq.merged.7_gN_Kmer.types miR-7_purification2_40percent_rep2.fastq.merged.7_gN_Kmer.types miR-7_purification2_40percent_rep3.fastq.merged.7_gN_Kmer.types miR-7_purification2_40percent_rep4.fastq.merged.7_gN_Kmer.types > miR-7_purification2_40percentM.fastq.merged.7_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-7_purification2_4percent_rep1.fastq.merged.7_gN_Kmer.types miR-7_purification2_4percent_rep2.fastq.merged.7_gN_Kmer.types miR-7_purification2_4percent_rep3.fastq.merged.7_gN_Kmer.types miR-7_purification2_4percent_rep4.fastq.merged.7_gN_Kmer.types > miR-7_purification2_4percentM.fastq.merged.7_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-7_purification2_inputlibrary_rep1.fastq.merged.7_gN_Kmer.types miR-7_purification2_inputlibrary_rep2.fastq.merged.7_gN_Kmer.types   > miR-7_purification2_inputlibraryM.fastq.merged.7_gN_Kmer.types

awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' let-7a_miR-155_miR-124_lsy-6_purification1_inputlibrary2.fastq.merged.155_gN_Kmer.types let-7a_miR-155_purification1_inputlibrary.fastq.merged.155_gN_Kmer.types   > miR-155_purification1_inputlibraryM.fastq.merged.155_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-155_purification1_0.4percent.fastq.merged.155_gN_Kmer.types    > miR-155_purification1_0.4percentM.fastq.merged.155_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-155_purification1_0percent.fastq.merged.155_gN_Kmer.types    > miR-155_purification1_0percentM.fastq.merged.155_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-155_purification1_1.265percent.fastq.merged.155_gN_Kmer.types    > miR-155_purification1_1.265percentM.fastq.merged.155_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-155_purification1_12.65percent.fastq.merged.155_gN_Kmer.types    > miR-155_purification1_12.65percentM.fastq.merged.155_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-155_purification1_40percent.fastq.merged.155_gN_Kmer.types    > miR-155_purification1_40percentM.fastq.merged.155_gN_Kmer.types
awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' miR-155_purification1_4percent.fastq.merged.155_gN_Kmer.types    > miR-155_purification1_4percentM.fastq.merged.155_gN_Kmer.types

#merge concentrations
for f in `ls * -d`; do awk 'BEGIN{FS=OFS="\t"}{split(FILENAME,name,"_");sites[$1]=1;data[$1,substr(name[3],1,2)]=$2}END{print "site\tInput\t0\t0.004\t0.01265\t0.04\t0.1265\t0.40";for (site in sites){print site"\t"data[site,"in"]"\t"data[site,"0p"]"\t"data[site,"0."]"\t"data[site,"1."]"\t"data[site,"4p"]"\t"data[site,"12"]"\t"data[site,"40"]}}' $f/* > $f.Brep1.counts; done

#MLE on station
for piRISC in 0.05 0.1 0.2 0.5 1 2 5 10 20;
do for KD in 0.05 0.1 0.2 0.5 1 2 5 10 20;
do python2 ./MLE_KD.py ./ let.Brep1.counts 100 $piRISC $KD 0;
done; done


#merge all 81 iterations for all sites Kmers
for type in 1 7 124 155 let lsy; do awk 'BEGIN{FS=OFS="\t"}(NR==1){for(i=1;i<=NF;i++){$i=$i""};print $0}(FNR==2){print $0}' $type.Brep1.counts_RISC*_all > $type.Brep1.MLE.txt; done

#reorder rep1n2
for file in `ls *.Brep1.MLE.txt`; do awk 'BEGIN{FS=OFS="\t"}(FNR==1){split(FILENAME,f,".");for (i=1;i<=NF;i++){tgs[$i]=i};printf ("num");for(k=6;k<=11;k++){for (g=2;g<=(f[1]-k+1);g++){printf ("\t"$(tgs["g"g"-"(g+k-1)]))}};print ("\t"$(tgs["RISC_stock_concentration"])"\t"$(tgs["nosite"])"\t"$(tgs["background"]))}(FNR>1){printf ("num");for(k=6;k<=11;k++){for (g=2;g<=(f[1]-k+1);g++){printf ("\t"$(tgs["g"g"-"(g+k-1)]))}};print ("\t"$(tgs["RISC_stock_concentration"])"\t"$(tgs["nosite"])"\t"$(tgs["background"]))}' $file > $file.rdrd_kg; done

for file in `ls *.Brep1.MLE.txt`; do awk 'BEGIN{FS=OFS="\t"}(FNR==1){split(FILENAME,f,".");for (i=1;i<=NF;i++){tgs[$i]=i};printf ("num");for(g=2;g<=(f[1]-6+1);g++){for (k=6;k<=(f[1]-g+1);k++){if(k<=11){printf ("\t"$(tgs["g"g"-"(g+k-1)]))}}};print ("\t"$(tgs["RISC_stock_concentration"])"\t"$(tgs["nosite"])"\t"$(tgs["background"]))}(FNR>1){printf ("num");for(g=2;g<=(f[1]-6+1);g++){for (k=6;k<=(f[1]-g+1);k++){if(k<=11){printf ("\t"$(tgs["g"g"-"(g+k-1)]))}}};print ("\t"$(tgs["RISC_stock_concentration"])"\t"$(tgs["nosite"])"\t"$(tgs["background"]))}' $file > $file.rdrd_gk; done
