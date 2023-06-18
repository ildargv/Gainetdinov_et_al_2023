#1.trim the adapter
module load fastx_toolkit/0.0.14
for i in $(ls *.fastq); do bsub -q long -W 240 -n 1  -e $i.err "fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 15 -c -v -i $i > $i.trimmed"; done

#2.collapse reads
for i in $(ls *.trimmed); do  /project/umw_phil_zamore/common/pipelines/ms_1G_1core "awk '(FNR%4==2){a[toupper(\$1)]++}END{for (i in a){print \">\"a[i];print i}}' $i > $i.collapsed.fa"; done

#3.add names
for i in $(ls *luc*.collapsed.fa); do echo $i; awk 'BEGIN{FS=OFS="\t"}(FNR==NR){names[$1]=$2","names[$1]}(FNR<NR)&&(FNR%2==1){num=substr($1,2);getline;if(names[$1]!=""){a=0;t=0;g=0;c=0;for (i=1;i<=length($1);i++){if(substr($1,i,1)=="G"){g++};if(substr($1,i,1)=="C"){c++};if(substr($1,i,1)=="A"){a++};if(substr($1,i,1)=="T"){t++}};data[names[$1]]=num"\t"$1"\t"a"\t"g"\t"c"\t"t}}END{for (i in data){print i"\t"data[i]}}' lucnames.txt $i > $i.stat;done 

for i in $(ls *let*.collapsed.fa); do echo $i; awk 'BEGIN{FS=OFS="\t"}(FNR==NR){names[$1]=$2","names[$1]}(FNR<NR)&&(FNR%2==1){num=substr($1,2);getline;if(names[$1]!=""){a=0;t=0;g=0;c=0;for (i=1;i<=length($1);i++){if(substr($1,i,1)=="G"){g++};if(substr($1,i,1)=="C"){c++};if(substr($1,i,1)=="A"){a++};if(substr($1,i,1)=="T"){t++}};data[names[$1]]=num"\t"$1"\t"a"\t"g"\t"c"\t"t}}END{for (i in data){print i"\t"data[i]}}' let7names.txt $i > $i.stat;done 

for i in $(ls *kctd7*.collapsed.fa); do echo $i; awk 'BEGIN{FS=OFS="\t"}(FNR==NR){names[$1]=$2","names[$1]}(FNR<NR)&&(FNR%2==1){num=substr($1,2);getline;if(names[$1]!=""){a=0;t=0;g=0;c=0;for (i=1;i<=length($1);i++){if(substr($1,i,1)=="G"){g++};if(substr($1,i,1)=="C"){c++};if(substr($1,i,1)=="A"){a++};if(substr($1,i,1)=="T"){t++}};data[names[$1]]=num"\t"$1"\t"a"\t"g"\t"c"\t"t}}END{for (i in data){print i"\t"data[i]}}' Kctd7names.txt $i > $i.stat;done 

for i in $(ls *l1mc*.collapsed.fa); do echo $i; awk 'BEGIN{FS=OFS="\t"}(FNR==NR){names[$1]=$2","names[$1]}(FNR<NR)&&(FNR%2==1){num=substr($1,2);getline;if(names[$1]!=""){a=0;t=0;g=0;c=0;for (i=1;i<=length($1);i++){if(substr($1,i,1)=="G"){g++};if(substr($1,i,1)=="C"){c++};if(substr($1,i,1)=="A"){a++};if(substr($1,i,1)=="T"){t++}};data[names[$1]]=num"\t"$1"\t"a"\t"g"\t"c"\t"t}}END{for (i in data){print i"\t"data[i]}}' L1MCnames.txt $i > $i.stat;done 


#4.merge: put each set of *.stat files belonging to one time course into one folder, name each folder as follows: e.g., l1mc_mili_long_1 l1mc_mili_long_2 l1mc_mili_long_3 l1mc_mili_short_1 l1mc_mili_short_2 l1mc_mili_short_3 (long for 0-960min trials, short for 0-8min trials)

#4a.merge 0-960min trials
for i in `ls | grep -v short`; do echo $i; awk -v f="$i" 'BEGIN{FS=OFS="\t";
points[1]="000min";
points[2]="020min";
points[3]="060min";
points[4]="120min";
points[5]="240min";
points[6]="480min";
points[7]="960min";
print "target\t0\t20\t60\t120\t240\t480\t960\t1000000" > f".counts";
print "target\t0\t20\t60\t120\t240\t480\t960\t1000000" > f".ppm"}{point=substr(FILENAME,index(FILENAME,"min_")-3,6);targets[$1]=1;total[point]+=$2;data[$1,point]=$2}END{for (t in targets){printf (t) > f".counts"; printf (t) > f".ppm"; for (i=1;i<=7;i++){printf ("\t"data[t,points[i]]) > f".counts"; printf ("\t"data[t,points[i]]*1000000/total[points[i]]) > f".ppm"};printf("\n") > f".counts"; printf ( "\t" ((data[t,"480min"]*1000000/total["480min"])+(data[t,"960min"]*1000000/total["960min"]))/((data[t,"000min"]*1000000/total["000min"])+(data[t,"020min"]*1000000/total["020min"]))) > f".ppm"; printf("\n") > f".ppm"}}' ./$i/*.stat ; done

#4b.merge 0-8min trials
for i in `ls | grep short`; do echo $i; awk -v f="$i" 'BEGIN{FS=OFS="\t";
points[1]="000min";
points[2]="001min";
points[3]="002min";
points[4]="004min";
points[5]="008min";
print "target\t0\t1\t2\t4\t8\t1000000" > f".counts";
print "target\t0\t1\t2\t4\t8\t1000000" > f".ppm"}{point=substr(FILENAME,index(FILENAME,"min_")-3,6);targets[$1]=1;total[point]+=$2;data[$1,point]=$2}END{for (t in targets){printf (t) > f".counts"; printf (t) > f".ppm"; for (i=1;i<=5;i++){printf ("\t"data[t,points[i]]) > f".counts"; printf ("\t"data[t,points[i]]*1000000/total[points[i]]) > f".ppm"};printf("\n") > f".counts"; printf ( "\t" ((data[t,"004min"]*1000000/total["004min"])+(data[t,"008min"]*1000000/total["008min"]))/((data[t,"000min"]*1000000/total["000min"])+(data[t,"001min"]*1000000/total["001min"]))) > f".ppm"; printf("\n") > f".ppm"}}' ./$i/*.stat ; done

#5.normalize using g2:g10:2mm (7 or fewer matches to target)
for RISC in kctd7_mili kctd7_miwi l1mc_ago2 l1mc_ago2_noG l1mc_efpiwi l1mc_efpiwi_noG l1mc_mili l1mc_mili_noG l1mc_miwi l1mc_miwi_noG let7_efpiwi let7_mili let7_miwi luc_mili luc_miwi; do echo $RISC; cd $RISC; for i in `find . -name "*.ppm" | grep -v short`; do awk 'BEGIN{FS=OFS="\t"; print "target\t0\t20\t60\t120\t240\t480\t960"}(FNR>1){if(index($1,"g2:g10:2mm:")>0){for (i=1;i<=7;i++){norm[i]+=$(i+1)}};lines[$0]=1}END{for (line in lines){split(line,fields,"\t");printf (fields[1]); for (i=1;i<=7;i++){printf ("\t"fields[i+1]/norm[i])};printf("\n")}}' $i > $i.norm ; done; for i in `find . -name "*.ppm" | grep short`; do awk 'BEGIN{FS=OFS="\t"; print "target\t0\t1\t2\t4\t8"}(FNR>1){if(index($1,"g2:g10:2mm:")>0){for (i=1;i<=5;i++){norm[i]+=$(i+1)}};lines[$0]=1}END{for (line in lines){split(line,fields,"\t");printf (fields[1]); for (i=1;i<=5;i++){printf ("\t"fields[i+1]/norm[i])};printf("\n")}}' $i > $i.norm ; done; cd ..; done

#6.infer relative product abundance by deducting relative substrate abundance from 1
for i in `ls *.norm`; do echo $i; awk 'BEGIN{FS=OFS="\t"}(FNR==1){print}(FNR>1){printf ($1);for(i=1;i<=(NF-1);i++){printf("\t"($2-$(i+1))/$2)};printf ("\n")}' $i > $i.prdct; done

#7.merge replicates
for RISC in kctd7_mili kctd7_miwi l1mc_ago2 l1mc_ago2_noG l1mc_efpiwi l1mc_efpiwi_noG l1mc_mili l1mc_mili_noG l1mc_miwi l1mc_miwi_noG let7_efpiwi let7_mili let7_miwi luc_mili luc_miwi; do echo $RISC; cd $RISC; gawk 'BEGIN{FS=OFS="\t"}(NF==8){a[$1][FILENAME]=$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;a["reference"][FILENAME]="0\t0\t0\t0\t0\t0"}(NF==6){a[$1][FILENAME]=$3"\t"$4"\t"$5"\t"$6;a["reference"][FILENAME]="0\t0\t0\t0"}END{for (r in a){printf r"\t0"; for (f in a["reference"]){if (a[r][f]==""){printf ("\t0\t0\t0\t0\t0\t0\t0\t0")}else{printf "\t"a[r][f]}};printf ("\n")}}' *.norm.prdct | grep -v reference | awk 'BEGIN{FS=OFS="\t"}($1=="target"){print}($1!="target"){lines[$0]=1}END{for (line in lines){print line}}' > ../$RISC.allreps.txt; cd ..; done

#8.run fitting in python3
for i in `ls *.allreps.txt`; do /pi/phillip.zamore-umw/home/ildar.gainetdinov-umw/common/pipelines/ms_short_8G "python3 /home/ildar.gainetdinov-umw/data/cns/processed_data/CleaveNseq_burst+ss_E_05_1_nomax_xy.py $i > $i.fit"; done
