#!/bin/bash

awk 'BEGIN { FS=OFS="\t" } 
    NR == 1 { split("", arr); for ( i =2; i<=NF; i++ ) { arr[i]=$i } } 
    NR == FNR { print } 
    NR != FNR && FNR == 1 { split("", ar); for (i =2; i<=NF; i++ ) { ar[$i]=i }  } 
    NR != FNR && FNR > 1 { line=$1; for (idx in arr) { n=ar[arr[idx]] ;line = line OFS $n } ; print line } 
 ' ./MatriceExonPSI_Patients.tsv ./MatriceExonTPM_Patients.tsv > ./MatriceExonMIX_Patients.tsv
 
awk 'BEGIN { FS=OFS="\t" } 
    NR == 1 { split("", arr); for ( i =2; i<=NF; i++ ) { arr[i]=$i } } 
    NR == FNR { print } 
    NR != FNR && FNR == 1 { split("", ar); split("", all); for (i =2; i<=NF; i++ ) { if ( all[$i] == 1 ) { $i=$i ".1" } ; all[$i]=1; ar[$i]=i }  } 
    NR != FNR && FNR > 2 { line=$1; for (idx in arr) { n=ar[arr[idx]] ;line = line OFS $n } ; print line } 
 ' ./MatriceExonTPM_CellLines.tsv ./MatriceExonPSI_CellLines.tsv  > MatriceExonMIX_CellLines.tsv

