#get site types for Kd gN_Kmer matches: parallel computing
for i in $(find . -name "*.trimmed"); do echo $i; split -a 4 -d -l 4000000 $i $i.xx. ; done 
for f in $(ls *.trimmed.xx.????); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR%4==2){if((length(\$1)==20)&&(index(\$1,\"N\")==0)){sites=0;site_loc=\"\";full=\"GATC\"\$1\"TGGA\";
for(k=6;k<=11;k++){for (g=2;g<=(22-k);g++){num=split( full, temp, substr(\"ACTATACAACCTACTACCTC\",(23-k)-g,k)); if(num>2){sites=2}; if((num==2)&&(substr(temp[1],length(temp[1]),1)!=substr(\"NACTATACAACCTACTACCTCN\",(23-k)-g,1))&&(substr(temp[2],1,1)!=substr(\"NACTATACAACCTACTACCTCN\",24-g,1)))
{site_loc=\"g\"g\"-\"(g+k-1);sites++}}};if(sites==1){counts[site_loc]++};if(sites==0){counts[\"nosite\"]++}}}END{for (site in counts){print site\"\\t\"counts[site]}}' $f > $f.gN_Kmer.types"; done
for f in $(ls *.trimmed); do echo $f; awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' $f.xx.????.gN_Kmer.types > $f.Kmer.types; done

#get  site types for Kd t1N_g2_Kmer matches: parallel computing 
for f in $(ls *.trimmed.xx.????); do /project/umw_phil_zamore/common/pipelines/ms_1G_1core "awk 'BEGIN{FS=OFS=\"\\t\"}(FNR%4==2){if((length(\$1)==20)&&(index(\$1,\"N\")==0)){sites=0;site_loc=\"\";full=\"GATC\"\$1\"TGGA\";
for(k=6;k<=10;k++){num=split( full, temp, substr(\"ACTATACAACCTACTACCTC\",(23-k)-2,k)); if(num>2){sites=2}; if((num==2)&&(substr(temp[1],length(temp[1]),1)!=substr(\"NACTATACAACCTACTACCTCN\",(23-k)-2,1)))
{site_loc=\"g2-\"(2+k-1)\":t1\"substr(temp[2],1,1);sites++}};if(sites==1){counts[site_loc]++};if(sites==0){counts[\"nosite\"]++}}}END{for (site in counts){print site\"\\t\"counts[site]}}' $f > $f.t1N_g2_Kmer.types"; done
for f in $(ls *.trimmed); do echo $f; awk 'BEGIN{FS=OFS="\t"}{sites[$1]+=$2}END{for (site in sites){print site"\t"sites[site]}}' $f.xx.????.t1N_g2_Kmer.types > $f.t1N_g2_Kmer.types; done

