# Protein composition of extracellular vesicles from malignant ascites defines the key components of tumor microenvironment and predicts prognosis of ovarian cancer patients

## Authors

Vyhlídalová Kotrbová A, Gömöryová K, Mikulová A, Kravec M, Plešingerová H, Potěšil D, Blériot C, Bied M, Dunsmore G, Kotouček J, Bednaříková M, Hausnerová J, Minář L, Crha I, Jandáková E, Zdráhal Z, Ginhoux F, Weinberger V, Bryja V, Hlaváčková Pospíchalová V*

* Corresponding author: pospich@sci.muni.cz 

## Reproducing the data analysis from this article

### Deposition of raw data to PRIDE
Raw proteomic data were deposited to the ProteomeXchange Consortium via PRIDE and can be accessed by identifier [PXD041751](https://www.ebi.ac.uk/pride/archive?keyword=PXD041751).

### Analysis of MS data

This repository contains three folders required to obtain the tables and figures from the manuscript:

`data`: contains input tables for data analysis

  - proteinGroups.txt - table containing mass spectrometry data (output of MaxQuant)
  - DSV_analysis_input - DLS measurement data
  - MISEV2018_protein_categories - list of EV markers based on MISEV 2018 ([Thery et al., 2018](https://pubmed.ncbi.nlm.nih.gov/30637094/))
  - preys-latest - protein localization assignment based on the [Human Cell Map](https://humancellmap.org/) resource
  - RNAseq_S2_table - supplementary table 2 from [Izar et al. (2020)](https://www.nature.com/articles/s41591-020-0926-0)
  - Izar_cell_markers_updated - list of cell type markers provided by Izar
  
`code`: contains individual scripts

  - 01_DLS-data.R (*Fig. S1B*)
  - 02_data-cleaning_genenames-update.R 
  - 03_data-processing.R 
  - 04_MISEV-markers-mapping.R (*Fig. 1C*, *Suppl. Table 3*)
  - 05_filter-out-B-samples.R (*Fig. 1B*)
  - 06_methods-intersection.R (*Fig. 2B, S2A, S2A'*, *Suppl. Table 4*)
  - 07_intersections_101.R (*Fig. 2C, 2D, S2B*, *Suppl. Table 5*)
  - 08_intersections_392.R (*Fig. 2C, 2D, S2B*, *Suppl. Table 6*)
  - 09_controls-addition.R (*Suppl. Table 7*)
  - 10_RNAseq_comparison_a.R (*Fig. 4B, 4C*)
  - 10_RNAseq_comparison_b.R (*Fig. 5A, 5B', 5C*)
  - 11_flow_cytometric_data_transformation.R 
  - 12_flow_cytometric_data_clustering_and_annotation.R (*Fig. 5C'*, *Suppl. Fig. 3A*)
  - 13_flow_cytometric_data_generating_of_final_plots.R (*Fig. 5C, 5C'*)
  
`outputs`: resulting Figures/Tables (some have been adjusted further in Affinity Designer software)