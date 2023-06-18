#gN-... contiguous pairing with all t1 nucleotides for g2-,g3-,g4-.... 317
awk 'BEGIN{temp = "AGGACCTACAGCTTGCATAGTATTTAGTTA";
for (g_st=2;g_st<=15;g_st++){for (g_end=g_st+5;g_end<=30;g_end++){target="";for (i=2;i<=30;i++){if ((i<g_st)||(i>g_end)){
if (substr(temp,31-i,1)=="A"){target="T"target};
if (substr(temp,31-i,1)=="T"){target="A"target};
if (substr(temp,31-i,1)=="G"){target="C"target};
if (substr(temp,31-i,1)=="C"){target="G"target}}else{target=substr(temp,31-i,1) target}};
if(g_st==2){print target"A\tg"g_st":g"g_end":cont:t1A";
print target"G\tg"g_st":g"g_end":cont:t1G";
print target"C\tg"g_st":g"g_end":cont:t1C";
print target"T\tg"g_st":g"g_end":cont:t1U"}else{print target"A\tg"g_st":g"g_end":cont:t1A"}}}}' > L1MC_gNcontig.txt

#g2-... one mononucleotide mismatch of all geometries 1197 (all t1A)
awk 'BEGIN{temp = "AGGACCTACAGCTTGCATAGTATTTAGTTA";
for (g_end=10;g_end<=30;g_end++){for (m1=2;m1<=g_end;m1++){target1="";target2="";target3="";for (i=2;i<=30;i++){if (i==m1){
if (substr(temp,31-i,1)=="A"){target1="T"target1;name1="UU";target2="G"target2;name2="UG";target3="C"target3;name3="UC"};
if (substr(temp,31-i,1)=="T"){target1="A"target1;name1="AA";target2="G"target2;name2="AG";target3="C"target3;name3="AC"};
if (substr(temp,31-i,1)=="G"){target1="C"target1;name1="CC";target2="A"target2;name2="CA";target3="T"target3;name3="CU"};
if (substr(temp,31-i,1)=="C"){target1="G"target1;name1="GG";target2="A"target2;name2="GA";target3="T"target3;name3="GU"}}else{if(i>g_end){
if (substr(temp,31-i,1)=="A"){target1="T"target1;target2="T"target2;target3="T"target3};
if (substr(temp,31-i,1)=="T"){target1="A"target1;target2="A"target2;target3="A"target3};
if (substr(temp,31-i,1)=="G"){target1="C"target1;target2="C"target2;target3="C"target3};
if (substr(temp,31-i,1)=="C"){target1="G"target1;target2="G"target2;target3="G"target3}}else{
target1=substr(temp,31-i,1) target1;
target2=substr(temp,31-i,1) target2;
target3=substr(temp,31-i,1) target3}}};
print target1"A\tg2:g"g_end":1mm:"m1":"name1;
print target2"A\tg2:g"g_end":1mm:"m1":"name2;
print target3"A\tg2:g"g_end":1mm:"m1":"name3}}}' > L1MC_g2mono.txt

#g2-... two mononucleotide mismatches of a single geometry 3976 (all t1A)
awk 'BEGIN{temp = "AGGACCTACAGCTTGCATAGTATTTAGTTA";
for (g_end=10;g_end<=30;g_end++){for (m1=2;m1<=g_end;m1++){for (m2=2;m2<=g_end;m2++){if(m1<m2){target="";for (i=2;i<=30;i++){if ((i==m1)||(i==m2)||(i>g_end)){
if (substr(temp,31-i,1)=="A"){target="T"target};
if (substr(temp,31-i,1)=="T"){target="A"target};
if (substr(temp,31-i,1)=="G"){target="C"target};
if (substr(temp,31-i,1)=="C"){target="G"target}}else{target=substr(temp,31-i,1) target}};
print target"A\tg2:g"g_end":2mm:"m1":"m2}}}}}'  > L1MC_g2mono_x2.txt

#g2-... one trinucleotide mismatch of a single geometry 357 (all t1A)
awk 'BEGIN{temp = "AGGACCTACAGCTTGCATAGTATTTAGTTA";
for (g_end=10;g_end<=30;g_end++){for (m1=2;m1<=g_end-2;m1++){target="";for (i=2;i<=30;i++){if ((i==m1)||(i==(m1+1))||(i==(m1+2))||(i>g_end)){
if (substr(temp,31-i,1)=="A"){target="T"target};
if (substr(temp,31-i,1)=="T"){target="A"target};
if (substr(temp,31-i,1)=="G"){target="C"target};
if (substr(temp,31-i,1)=="C"){target="G"target}}else{target=substr(temp,31-i,1) target}};
print target"A\tg2:g"g_end":3mm:"m1}}}'  > L1MC_g2tri.txt

#g2-... one bulge in guide of a single geometry 378 (all t1A)
awk 'BEGIN{temp = "AGGACCTACAGCTTGCATAGTATTTAGTTA";
for (g_end=10;g_end<=30;g_end++){for (m1=2;m1<=g_end-1;m1++){target="";for (i=2;i<=30;i++){if (i>g_end){
if (substr(temp,31-i,1)=="A"){target="T"target};
if (substr(temp,31-i,1)=="T"){target="A"target};
if (substr(temp,31-i,1)=="G"){target="C"target};
if (substr(temp,31-i,1)=="C"){target="G"target}}else{if(i!=m1){target=substr(temp,31-i,1) target}}};
print target"A\tg2:g"g_end":gb:"m1}}}'  > L1MC_g2gb.txt

#g2-... one bulge in target of all geometries 1512 (all t1A)
awk 'BEGIN{temp = "AGGACCTACAGCTTGCATAGTATTTAGTTA";
for (g_end=10;g_end<=30;g_end++){for (m1=2;m1<=g_end-1;m1++){target1="";target2="";target3="";target4="";for (i=2;i<=30;i++){if(i>g_end){
if (substr(temp,31-i,1)=="A"){target1="T"target1;target2="T"target2;target3="T"target3;target4="T"target4};
if (substr(temp,31-i,1)=="T"){target1="A"target1;target2="A"target2;target3="A"target3;target4="A"target4};
if (substr(temp,31-i,1)=="G"){target1="C"target1;target2="C"target2;target3="C"target3;target4="C"target4};
if (substr(temp,31-i,1)=="C"){target1="G"target1;target2="G"target2;target3="G"target3;target4="G"target4}}else{
target1=substr(temp,31-i,1) target1;
target2=substr(temp,31-i,1) target2;
target3=substr(temp,31-i,1) target3;
target4=substr(temp,31-i,1) target4; if(i==m1){
target1="T"target1;name1="U";target2="G"target2;name2="G";target3="C"target3;name3="C";target4="A"target4;name4="A"}}};
print target1"A\tg2:g"g_end":tb:"m1+1":"name1;
print target2"A\tg2:g"g_end":tb:"m1+1":"name2;
print target3"A\tg2:g"g_end":tb:"m1+1":"name3;
print target4"A\tg2:g"g_end":tb:"m1+1":"name4}}}' > L1MC_g2tb.txt

#g3-...g4-...g5-... g6-...: one mono/di/trinucleotide mismatch of single geometry 2622(all t1A)
awk 'BEGIN{temp = "AGGACCTACAGCTTGCATAGTATTTAGTTA";
for (g_st=3;g_st<=6;g_st++){for (g_end=g_st+16;g_end<=30;g_end++){for (m1=g_st;m1<=g_end-2;m1++){target="";for (i=2;i<=30;i++){if ((i==m1)||(i==(m1+1))||(i==(m1+2))||(i<g_st)||(i>g_end)){
if (substr(temp,31-i,1)=="A"){target="T"target};
if (substr(temp,31-i,1)=="T"){target="A"target};
if (substr(temp,31-i,1)=="G"){target="C"target};
if (substr(temp,31-i,1)=="C"){target="G"target}}else{target=substr(temp,31-i,1) target}};
print target"A\tg"g_st":g"g_end":3mm:"m1}}}}'  > L1MC_gNtri.txt
awk 'BEGIN{temp = "AGGACCTACAGCTTGCATAGTATTTAGTTA";
for (g_st=3;g_st<=6;g_st++){for (g_end=g_st+16;g_end<=30;g_end++){for (m1=g_st;m1<=g_end-1;m1++){target="";for (i=2;i<=30;i++){if ((i==m1)||(i==(m1+1))||(i<g_st)||(i>g_end)){
if (substr(temp,31-i,1)=="A"){target="T"target};
if (substr(temp,31-i,1)=="T"){target="A"target};
if (substr(temp,31-i,1)=="G"){target="C"target};
if (substr(temp,31-i,1)=="C"){target="G"target}}else{target=substr(temp,31-i,1) target}};
print target"A\tg"g_st":g"g_end":2mm:"m1}}}}'  > L1MC_gNdi.txt
awk 'BEGIN{temp = "AGGACCTACAGCTTGCATAGTATTTAGTTA";
for (g_st=3;g_st<=6;g_st++){for (g_end=g_st+16;g_end<=30;g_end++){for (m1=g_st;m1<=g_end;m1++){target="";for (i=2;i<=30;i++){if ((i==m1)||(i<g_st)||(i>g_end)){
if (substr(temp,31-i,1)=="A"){target="T"target};
if (substr(temp,31-i,1)=="T"){target="A"target};
if (substr(temp,31-i,1)=="G"){target="C"target};
if (substr(temp,31-i,1)=="C"){target="G"target}}else{target=substr(temp,31-i,1) target}};
print target"A\tg"g_st":g"g_end":1mm:"m1}}}}'  > L1MC_gNmono.txt
#add library flanks
awk 'BEGIN{FS=OFS="\t"}{print "GTTCAGAGTTCTACAGTCCGACGATC"$1"TGGAATTCTCGGGTGCCAAGG"}' L1MC_*.txt > L1MC.txt
sort L1MC.txt | uniq > L1MCuniq.txt

awk 'BEGIN{FS=OFS="\t"}{print "GTTCAGAGTTCTACAGTCCGACGATC"$1"TGGAATTCTCGGGTGCCAAGG\t"$2}' L1MC_*.txt > L1MCwithnames.txt

awk 'BEGIN{FS=OFS="\t"}{s="GTTCAGAGTTCTACAGTCCGACGATC"$1"TGGAATTCTCGGGTGCCAAGG"; gsub("T","U",s); print s"\t"$2}' L1MC_*.txt > L1MCwithnames_RNA.txt

cat L1MC_*.txt > L1MCnames.txt
