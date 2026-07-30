PRADA-RNA_Structuromics
============

Understanding the spatial organization of biomolecules in living organisms is a central but challenging goal in biology. A nontoxic proximity labeling approach Peroxidase Reactions Activated by D-amino Acids (PRADA) is developed, which uses an engineered oxidase to convert non-proteinogenic D-amino acids into H2O2 for in situ activation of a genetically-fused peroxidase. PRADA is applied in C. elegans, Drosophila, zebrafish, and mice and has unique RNA labeling activity to probe RNA structures in vitro and in vivo. PRADA-RNA_Structuromics shows the analysis for RNA structures from PRADA datasets. 

Downstream analysis can be found in `analysis`. 

Software Pre-requisites
------------------------

  1. ShapeMapper2
  2. Python
  3. R

Installation Requirments
--------------------------

   1. snakemake == 7.15.2
   2. scipy >= 1.7.1
   3. scikit-learn >= 0.24.2
   4. numpy <= 1.20.3
   5. pandas >= 1.5.0
   6. matplotlib >= 3.5.1
   7. seaborn >= 0.11.2
   8. biopython >= 1.79
   9. cairosvg >= 2.7.0
   10. statannot >= 0.2.3

This will be automatically done during `pip install`

Usage
--------------------------

After completing the ShapeMapper2 pipeline, run the notebooks in the analysis/ directory to reproduce the downstream analyses and figures.


