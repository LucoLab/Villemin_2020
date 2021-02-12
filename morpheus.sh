#!/bin/bash

## NF ClaudinLow == 31
metadata="./data/metadata/TCGA_BRCA_metadata.txt"
for d in  ./resuts_* ; do
    echo $d
    r_type=$( echo $d | awk ' { split($0, arr, "/"); split(arr[2], ar, "_"); print ar[2] }' )
    head -n 1 ./data/MatriceExon${r_type}_Patients.tsv > ${d}/MatriceExon${r_type}_Patients_morpheus.tsv
    awk 'BEGIN { FS=OFS="\t"} NR > 1 && FNR==NR { clow[$1]=$31 } 
    /^name/ { line="ClaudinLow"; for (i=2; i<= NF; i++) {  line = line "\t" clow[$i] } ; print line } 
    ' $metadata ./data/MatriceExon${r_type}_Patients.tsv >> ${d}/MatriceExon${r_type}_Patients_morpheus.tsv
    awk 'BEGIN { FS=OFS="\t"} NR > 1 && FNR==NR { label[$1]=$2 } 
    /^name/ { line="Label"; for (i=2; i<= NF; i++) {  line = line "\t" label[$i] } ; print line } 
    ' ${d}/transfer/patients_labels.tsv ./data/MatriceExon${r_type}_Patients.tsv >> ${d}/MatriceExon${r_type}_Patients_morpheus.tsv
    awk 'BEGIN { FS=OFS="\t"} NR == FNR { elite[$1]=1 } NR != FNR && $1 in elite {  print } '  $d/boruta/elite_features.tsv ./data/MatriceExon${r_type}_Patients.tsv >> ${d}/MatriceExon${r_type}_Patients_morpheus.tsv
done        
