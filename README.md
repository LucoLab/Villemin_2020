
<p align="center">
<img align="center"   src="/img/cell2patient.png" alt="cell2patient Logo">
</p>


A cell-to-patient machine learning transfer approach uncovers novel basal-like breast cancer prognostic markers amongst alternative splice variants
=============

---

Data can be accessed here : [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4117806.svg)](https://doi.org/10.5281/zenodo.4117806)


## Quick overview

Two scripts are given separetely for splicing and expression.  


```shell

python  classification_cell2patient_splicing.py \
	 -c MatriceExonPSI_CellLines.csv \
	 -p MatriceExonPSI_Patients.csv \
	 -t 0.6 \
	 -n 1000 \ 
	 -a /path2test/


```shell

python  classification_cell2patient_expression.py \
	 -c MatriceExonTPM_CellLines.csv \
	 -p MatriceExonTPM_Patients.csv \
	 -t 0.6 \
	 -n 1000 

```

- **t :** Threshold for class probabilities.
- **c :** Path to a matrice with Expression/Splicing values for Cell Lines.
- **p :** Path to a matrice with Expression/Splicing values for Patients.
- **n :** Number of tree in the forest.
- **r :** Random state.
- **r :** Path to directory where a matrice of PSI for cell lines (different projects from training) is found. (this option is only available for splicing)


---
