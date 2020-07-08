# NCLDV_Biogeography
Code associated with Endo et al. 2020, Nature Ecology &amp; Evolution

This git repository contains the files and scripts that are used to generate normalized abundance profile of NCLDVs for the following manuscript

```
Biogeography of marine giant viruses reveals their interplay with eukaryotes and ecological functions
by
Endo H, Blanc-Mathieu R, Li Y, Salazar G, Henry N, Labadie K, de Vargas C, Sullivan MB, Bowler C, Wincker P, Karp-Boss L, Sunagawa S & Ogata H
Nature Ecology & Evolution (in revision)
```

# Prerequisites
- MAFFT (v7.453 or higher)
- Python (v3.7.3 or higher)
- pplacer (v1.1.alpha19)
- R (v.3.6.2 or higher)

# Data
Input files for pplacer.
- ```pplacer_ref211.msa```: Multiple alignments of PolB reference sequences for pplacer analysis
- ```RAxML_result.pplacer.nwk```: Phylogenetic tree (newick format) constructed from the reference sequences
- ```RAxML_info.pplacer.txt```: Information on the reference phylogenetic tree
- The original files of sequences and the abundance profile of the OM-RGC.v2 can be downloaded from https://www.ocean-microbiome.org

Input files for data standardization
- ```PolB_NCLDV_table.txt```: Raw frequency table of the NCLDV PolB genes


# Scripts
- ```run_pplacer.sh```: Main script for pplacer
  - ```parse_jplace.py```: Parse jplacer file and output sequence IDs
  - ```extract_target_line.py```: Extract target lines from the abundance profiles
- ```ncldv_stdz.R```: R script to standardize the dabundance profile

# Contact
- Hisashi Endo - endo[@]scl.kyoto-u.ac.jp
