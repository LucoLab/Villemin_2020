
<p align="center">
<img align="center"   src="/img/cell-to-patient.png" alt="cell2patient Logo">
</p>


A cell-to-patient machine learning transfer approach uncovers novel basal-like breast cancer prognostic markers amongst alternative splice variants
=============


---

<p align="center">1. Main</p>


## Quick overview



```shell

python  classification_semi_expression.py \
	 -c ${dir_analysis}MatriceExonTPM_CellLines.csv \
	 -p ${dir_analysis}MatriceExonTPM_Patients.csv \
	 -t 0.6 \
	 -n 1000 

```

- **t :** Threshold for class probabilities.
- **c :** Path to a matrice with Expression/Splicing values for Cell Lines.
- **p :** Path to a matrice with Expression/Splicing values for Patients.
- **n :** Number of tree in the forest.



__Very Important__  

---
