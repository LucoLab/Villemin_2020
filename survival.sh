#!/bin/bash

for d in  ./resuts_* ; do
    r_type=$( echo $d | awk ' { split($0, arr, "/"); split(arr[2], ar, "_"); print ar[2] }' )
    awk '$2=="BASALB" {print $1}' ${d}/transfer/patients_labels.tsv  > ${d}/basalb.ls
    awk '$2=="BASALA" {print $1}' ${d}/transfer/patients_labels.tsv  > ${d}/basala.ls
    ./R_scripts/Survival_AvsB.R -s ./data/metadata/TCGA_CDR.csv -b ./resuts_${r_type}/basala.ls -a ./resuts_${r_type}/basalb.ls -o ./${r_type}
    rm ${d}/basalb.ls ${d}/basala.ls
done

