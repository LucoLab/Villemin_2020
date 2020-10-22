
<p align="center">
<img align="center"   src="/img/cell2patient.png" alt="cell2patient Logo">
</p>


A cell-to-patient machine learning transfer approach uncovers novel basal-like breast cancer prognostic markers amongst alternative splice variants
=============

---

All data to reproduce figures can be accessed here : [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4117806.svg)](https://doi.org/10.5281/zenodo.4117806)


## How to use

Two python(3) scripts are given separately for splicing and expression.  

In directory data, you will find the input files.  

They are based on the following version of scikit-learn (0.21.2.)  

NB: Imputer warnings when script start is not an error.

They call one R script to plot survival over the rounds of  classification.


```shell

python  classification_cell2patient_splicing.py \
	 -c {absolutepath2}/MatriceExonPSI_CellLines.csv \
	 -p {absolutepath2}/MatriceExonPSI_Patients.csv \
	 -t 0.6 \
	 -n 1000 \ 
```

```shell

python  classification_cell2patient_expression.py \
	 -c {absolutepath2}/MatriceExonTPM_CellLines.csv \
	 -p {absolutepath2}/MatriceExonTPM_Patients.csv \
	 -t 0.6 \
	 -n 1000 

```

- **t :** Threshold for class probabilities.
- **c :** Path to a matrice with Expression/Splicing values for Cell Lines.
- **p :** Path to a matrice with Expression/Splicing values for Patients.
- **n :** Number of tree in the forest.



The final file annotated is splicing_TCGA_BASAL_HEADER_ADDED.tsv.  
You can visualize using https://software.broadinstitute.org/morpheus.

The best features of interest are in outputBorutaPy.txt/.bed.  

---
