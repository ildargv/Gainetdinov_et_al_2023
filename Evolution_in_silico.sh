#1.trim the adapter from E16.5 data
module load fastx_toolkit/0.0.14
for i in $(ls *.fastq); do bsub -q short -W 240 -n 4  -e $i.err "fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 18 -c -v -i $i > $i.trimmed"; done

#2.filter at p20q100
module load fastx_toolkit/0.0.14
for i in $(find . -maxdepth 1 -name "*.trimmed"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "fastq_quality_filter -q 20 -p 100 -i $i -o $i.q20p100" ; done 

#3.remove 3N and keep 24+nt reads
for i in $(find . -maxdepth 1 -name "*.trimmed.q20p100"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "awk '(NR%4==1){name=\$1;getline;len=length(\$1);if(len>=27){print name;print(substr(\$1,1,len-3));getline;print;getline;print(substr(\$1,1,len-3))}}' $i > $i.3N.fq" ; done 

#4.remove SILVA.BK000964.3.chrM.polyN with 1 mismatch allowed
module load bowtie/1.0.0
for i in $(find . -name "*.fq"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "bowtie --un $i.x_filter_1mm.fastq -k 1 -v 1 SILVA.BK000964.3.chrM.polyN $i > /dev/null" ; done 

#5.send to Tailor
for i in `find . -name "*.fq.x_filter_1mm.fastq"`; do /project/umw_phil_zamore/common/pipelines/ms "/project/umw_phil_zamore/common/pipelines/Tailor/run_tailing_pipeline.sh -i $i -c 8 -g mm10.fa -o $i.Tailor" ; done

#6.convert bed2 files into rpm files
for i in $(find . -name "*.p20.bed2"); do  /project/umw_phil_zamore/common/pipelines/ms_short_8G "awk '((\$9==0)&&(dealt[\$7]!=1)){nontailed[substr (\$7, 1, length (\$7)-\$9)]+=\$4; all[substr (\$7, 1, length (\$7)-\$9)]=1;coordinates[substr (\$7, 1, length (\$7)-\$9)]=\$1\":\"\$2\"-\"\$3\"\\t\"\$6; alluniq+=\$4; dealt[\$7]=1}((\$9>0)&&(dealt[\$7]!=1)){tailed[substr (\$7, 1, length (\$7)-\$9)]+=\$4;all[substr (\$7, 1, length (\$7)-\$9)]=1 ;coordinates[substr (\$7, 1, length (\$7)-\$9)]=\$1\":\"\$2\"-\"\$3\"\\t\"\$6; alluniq+=\$4; dealt[\$7]=1}END{print \"coordinates\\tstrand\\tsequenceASis\\ttotalRPM\\ttailedRPM\"; for (r in all){print coordinates[r]\"\\t\"r\"\\t\"(tailed[r]+nontailed[r])*1000000/alluniq\"\\t\"tailed[r]*1000000/alluniq}}' $i |sort -k5 -g -r > $i.rpm"; done

#7.merge two reps into one
awk 'BEGIN{FS=OFS="\t"}($4>10)&&(length($3)>=26)&&(length($3)<=30){ago[substr($3,2,20)]=1;piwi[substr($3,2,25)]=1}END{for (guide in ago){target=""; for (j=length(guide);j>=1;j--){if(substr(guide,j,1)=="A"){target=target"T"};if(substr(guide,j,1)=="T"){target=target"A"};if(substr(guide,j,1)=="C"){target=target"G"};if(substr(guide,j,1)=="G"){target=target"C"}};print ">"guide > "16.5.g2g21.fa"; print target > "16.5.g2g21.fa"};for (guide in piwi){target=""; for (j=length(guide);j>=1;j--){if(substr(guide,j,1)=="A"){target=target"T"};if(substr(guide,j,1)=="T"){target=target"A"};if(substr(guide,j,1)=="C"){target=target"G"};if(substr(guide,j,1)=="G"){target=target"C"}};print ">"guide > "16.5.g2g26.fa"; print target > "16.5.g2g26.fa"}}' *.rpm 

#8.resolve ambiguities in consensus sequence of L1Md repeats from RepBase
for file in `ls *.fasta`; do awk '(FNR==1){print; getline;
gsub("W","A",$1);
gsub("S","C",$1);
gsub("M","A",$1);
gsub("K","G",$1);
gsub("R","A",$1);
gsub("W","A",$1);
gsub("Y","C",$1);
gsub("B","C",$1);
gsub("D","A",$1);
gsub("H","A",$1);
gsub("V","A",$1);
gsub("N","A",$1);print}' $file > $file.resolved; done

#9.create codons table to allow only synonymous mutations
awk 'BEGIN{FS=OFS="\t"}{n=split($2,a,",");for(i=1;i<=n;i++){print $1"\t"a[i]}}' codons.txt > codons_table.txt

#10.run 100 simulations of 1000 single-nucleotide substitutions
for iter in {1..100}; do for repeat in L1MdA_I L1MdGf_I L1MdTf_I; do echo $repeat; awk -v seed=$RANDOM 'BEGIN{srand(seed);FS=OFS="\t"}(FNR==NR){codons[$2]=$1}(FNR<NR)&&(FNR==1){n=split($1,coord,"_");getline;print ">0";print;len=length($1);seq=$1;for (i=1;i<=1000;i++){found=0;
while(found==0){pos=int(rand()*len);if (pos==0){pos=1};if (pos>len){pos=len};
type=rand();
if (substr(seq,pos,1)=="A"){if(type<0.189){new="T"}else{if(type>0.863){new="C"}else{new="G"}}};
if (substr(seq,pos,1)=="T"){if(type<0.189){new="A"}else{if(type>0.863){new="G"}else{new="C"}}};
if (substr(seq,pos,1)=="C"){if(type<0.148){new="A"}else{if(type>0.284){new="T"}else{new="G"}}};
if (substr(seq,pos,1)=="G"){if(type<0.148){new="T"}else{if(type>0.284){new="A"}else{new="C"}}};
if ((pos>=coord[n-3])&&(pos<=coord[n-2])){frame=(pos-coord[n-3]+1)% 3;
if (frame==1){if(codons[substr(seq,pos,3)]==codons[new substr(seq,pos+1,2)]){found=1}};
if (frame==2){if(codons[substr(seq,pos-1,3)]==codons[substr(seq,pos-1,1) new substr(seq,pos+1,1)]){found=1}};
if (frame==0){if(codons[substr(seq,pos-2,3)]==codons[substr(seq,pos-2,2) new]){found=1}}}else{
if ((pos>=coord[n-1])&&(pos<=coord[n])){frame=(pos-coord[n-1]+1)%3;
if (frame==1){if(codons[substr(seq,pos,3)]==codons[new substr(seq,pos+1,2)]){found=1}};
if (frame==2){if(codons[substr(seq,pos-1,3)]==codons[substr(seq,pos-1,1) new substr(seq,pos+1,1)]){found=1}};
if (frame==0){if(codons[substr(seq,pos-2,3)]==codons[substr(seq,pos-2,2) new]){found=1}}}else{found=1}}};
print ">"i"_"pos"_"substr(seq,pos,1)"_"new;
if((pos!=1)&&(pos!=len)){seq=substr(seq,1,pos-1) new substr(seq,pos+1)};
if(pos==1){seq=new substr(seq,pos+1)};
if(pos==len){seq=substr(seq,1,len-1) new};
print seq }}' codons_table.txt $repeat.fasta.resolved > $repeat.$iter.evo; done; done

#11a.split 100 mutated sequenced into 20 sequence chunks for parallel computing
for file in `ls *.evo`; do echo $file; split -a 3 -d -l 40 $file $file.xx.; done

#11b.count siRNAs that cleave
for file in `ls *.evo.xx.???`; do /project/umw_phil_zamore/common/pipelines/ms_1G_1core_24h "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(FNR%2==0){targets[\$1]=1}(FNR<NR)&&(FNR%2==0){totaltargets=0;for(target in targets){for (pos=0;pos<=length(\$1)-20;pos++){if((substr(\$1,pos+9,1)==substr(target,9,1))&&(substr(\$1,pos+11,1)==substr(target,11,1))&&(substr(\$1,pos+12,1)==substr(target,12,1))&&(substr(\$1,pos+13,1)==substr(target,13,1))){mm=0;for (g=14;g<=20;g++){if(substr(\$1,pos+g,1)!=substr(target,g,1)){mm++}};if(mm<=1){for (g=1;g<=8;g++){if(substr(\$1,pos+g,1)!=substr(target,g,1)){mm++}};if(substr(\$1,pos+10,1)!=substr(target,10,1)){mm++};if(mm<=5){totaltargets++;break}}}}};print FNR/2\"\\t\"totaltargets}' 16.5.g2g21.fa $file > $file.ago.txt" ; done

#11c.count piRNAs that cleave
for file in `ls *.evo.xx.???`; do /project/umw_phil_zamore/common/pipelines/ms_1G_1core_24h "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR==NR)&&(FNR%2==0){targets[\$1]=1}(FNR<NR)&&(FNR%2==0){totaltargets=0;for(target in targets){for (pos=0;pos<=length(\$1)-25;pos++){mm=0;for (g=1;g<=25;g++){if(substr(\$1,pos+g,1)!=substr(target,g,1)){mm++};if(mm>6){break}};if(mm<=6){totaltargets++;break}}};print FNR/2\"\\t\"totaltargets}' 16.5.g2g26.fa $file > $file.piwi.6.txt" ; done

#11d.merge after parallel computing
for file in `ls *.evo`; do echo $file; awk 'BEGIN{FS=OFS="\t"}{n=split(FILENAME,f,".");print (f[n-2]*21)+$1"\t"$2}' $file.xx.???.ago.txt > $file.ago.txt; done
for file in `ls *.evo`; do echo $file; awk 'BEGIN{FS=OFS="\t"}{n=split(FILENAME,f,".");print (f[n-3]*21)+$1"\t"$2}' $file.xx.???.piwi.txt > $file.piwi.txt; done

#12.put all data together 
for srs in $(ls *.evo.ago.txt *.evo.piwi.txt); do awk '{a[FILENAME]=1;
reads[FILENAME,FNR]=$2}END{printf ("name");for (j=1;j<=1001;j++){printf ("\t");printf (j)};
printf ("\n"); for (i in a){printf (i);
 for (j=1;
j<=1001;
j++){printf ("\t");
printf (reads[i,j])};
printf ("\n")}}' $srs; done > all.agoNpiwi.txt
