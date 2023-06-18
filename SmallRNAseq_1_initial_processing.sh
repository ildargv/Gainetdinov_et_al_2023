#1.trim the adapter
module load fastx_toolkit/0.0.14
for i in $(ls *.fastq); do bsub -q short -W 240 -n 4  -e $i.err "fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 15 -c -v -i $i > $i.trimmed"; done

#2.filter at q20p100 (Phred >=20 at all positions) 
module load fastx_toolkit/0.0.14
for i in $(find . -maxdepth 1 -name "*.trimmed"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "fastq_quality_filter -q 20 -p 100 -i $i -o $i.q20p100" ; done 

#3.UMI:remove duplicates and anything <18nt
for i in $(find . -maxdepth 1 -name "*.trimmed.q20p100"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "awk '(NR%4==1){name=\$1}(NR%4==2){total++;if((substr(\$1,length(\$1)-11,3)==\"GTC\")&&(substr(\$1,length(\$1)-5,3)==\"TAG\")&&(((substr(\$1,4,3)==\"CGA\")&&(substr(\$1,10,3)==\"TAC\"))||((substr(\$1,4,3)==\"ATC\")&&(substr(\$1,10,3)==\"AGT\")))){umis++;if(a[\$1]!=1){nondup++;a[\$1]=1;if(length(\$1)>47){longer18++; print name; print substr (\$1,16,length (\$1)-30);getline; print; getline; print substr (\$1,16,length (\$1)-30)}}}}END{print FILENAME\"\\t\"total\"\\t\"umis\"\\t\"nondup\"\\t\"longer18 > FILENAME\".dup\"}' $i > $i.deUMI.dedup.fq" ; done 

#4.remove SILVA rRNA, genbank BK000964.3, chrM, polyN with 1 mismatch allowed
module load bowtie/1.0.0
Bowtie-build SILVA.BK000964.3.chrM.polyN.fa SILVA.BK000964.3.chrM.polyN
for i in $(find . -name "*.fq"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "bowtie --un $i.x_filter_1mm -k 1 -v 1 /nl/umw_phil_zamore/ildar/mm10_set_1/SILVA.BK000964.3.chrM.polyN $i > /dev/null" ; done 

#5.remove spikeins
module load bowtie/1.0.0
for i in $(find . -name "*.x_filter_1mm"); do /project/umw_phil_zamore/common/pipelines/ms_short_8G "bowtie --norc --un $i.x_spikein.fq -v 0 /nl/umw_phil_zamore/ildar/mm10_set_1/final.set $i | awk 'BEGIN{FS=OFS=\"\\t\"}(substr(\$3,index(\$3,\"-\")+1)==length(\$5)){spike[\$3]++}END{for (i in spike){print i\"\\t\"spike[i]}}' > $i.spikein" ; done 
	
#6.send to Tailor
for i in `find . -name "*.fq.x_filter_1mm.x_spikein.fq"`; do /project/umw_phil_zamore/common/pipelines/ms "/project/umw_phil_zamore/common/pipelines/Tailor/run_tailing_pipeline.sh -i $i -c 8 -g /project/umw_phil_zamore/common/pipelines/Tailor/indexes/mm10.fa -o $i.Tailor" ; done

#7.convert bed2 files from Tailor output to rpm files
for i in $(find . -name "*.p20.bed2"); do  /project/umw_phil_zamore/common/pipelines/ms_short_64G "awk 'BEGIN{FS=OFS=\"\\t\"}((\$9==0)&&(dealt[\$7]!=1)){nontailed[substr (\$7, 1, length (\$7)-\$9)]+=\$4; all[substr (\$7, 1, length (\$7)-\$9)]=1;coordinates[substr (\$7, 1, length (\$7)-\$9)]=\$1\":\"\$2\"-\"\$3\"\\t\"\$6; alluniq+=\$4; dealt[\$7]=1}((\$9>0)&&(dealt[\$7]!=1)){tailed[substr (\$7, 1, length (\$7)-\$9)]+=\$4;all[substr (\$7, 1, length (\$7)-\$9)]=1 ;coordinates[substr (\$7, 1, length (\$7)-\$9)]=\$1\":\"\$2\"-\"\$3\"\\t\"\$6; alluniq+=\$4; dealt[\$7]=1}END{print \"coordinates\\tstrand\\tsequenceASis\\ttotalRPM\\ttailedRPM\"; for (r in all){print coordinates[r]\"\\t\"r\"\\t\"(tailed[r]+nontailed[r])*1000000/alluniq\"\\t\"tailed[r]*1000000/alluniq}}' $i |sort -k5 -g -r > $i.rpm"; done

#8.bedgraph coverage for 5ends only
#8a.align onto mm10
module load bowtie/1.0.0
for i in $(find . -maxdepth 1 -name "*rpm"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "awk '(length(\$3)>23){print \">\"\$4-\$5\"-\"\$5; print \$3}' $i > $i.fa; bowtie -f -v 0 -a mm10 $i.fa | awk '{print \$3\"\\t\"\$4\"\\t\"(\$4+length(\$5))\"\\t\"\$1\"\\t\"\$7+1\"\\t\"\$2}' > $i.v0a.bed2.24; rm $i.fa"; done

#8b.make bedgraph for unique 5'ends
module load bedtools/2.26.0
for i in $(find . -maxdepth 1 -name "*.24"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "awk '(\$6==\"+\")&&(\$5==1){reads=substr(\$4,1,index(\$4,\"-\")-1)+substr(\$4,index(\$4,\"-\")+1);aplus[\$1,\$2]+=reads/\$5}(\$6==\"-\")&&(\$5==1){reads=substr(\$4,1,index(\$4,\"-\")-1)+substr(\$4,index(\$4,\"-\")+1);aminus[\$1,\$3]+=reads/\$5}END{for (ij in aplus) {split(ij,indices,SUBSEP); i=indices[1]; j=indices[2]; print i\"\\t\"j\"\\t\"j+1\"\\t\"aplus[ij] > FILENAME \".plus.rd.unique.bedgraph\"};for (ij in aminus) {split(ij,indices,SUBSEP); i=indices[1]; j=indices[2]; print i\"\\t\"j-1\"\\t\"j\"\\t-\"aminus[ij] > FILENAME \".minus.rd.unique.bedgraph\"}}' $i"; done

#8c.convert bedgraphs to bigWig
for i in $(find . -maxdepth 1 -name "*.bedgraph"); do /project/umw_phil_zamore/common/pipelines/ms_short_64G "sort -k1,1 -k2,2n $i > $i.srt; /project/umw_phil_zamore/common/pipelines/piPipes/bin/bedGraphToBigWig $i.srt mm10.genome $i.bigWig; rm $i.srt"; done