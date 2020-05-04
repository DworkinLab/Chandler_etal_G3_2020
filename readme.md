# Data and code for the analysis of Chandler et al. 2020 G3

This is repository contains all of the data and `R` scripts to reproduce the results from
Chandler, C.H., Mammel, A. & Dworkin, I. 2020. Sexual selection does not increase the rate of compensatory adaptation to a mutation influencing a secondary sexual trait in *Drosophila melanogaster*. G3: Genes, Genomes. (10):1541-1551/ https://doi.org/10.1534/g3.119.400934. Link to paper [here](https://www.g3journal.org/content/10/5/1541).

In the data folder there are 5 `.csv` files. Details about each below.

## Figure 2

`plot_dup_data.csv` is the data on which genomic duplications rescued the scd phenotype. The `gene_data.csv` simply lists genes and positions for the plot for Figure 2 summarizing this (Figure2_plot_mapping_viewports.R).


## Figure 3

`scd_RAL_data_final.csv` is the data for the F1 males from crosses between a panel of different wild type strains and the progenitor *scd[1]* strain. The script Figure3_plot_RAL.R uses this data to reproduce figure 3. The data columns for this can be interpreted as follows:

genotype - The strain used for the male parent (females always were scd strain so male F1 would be hemizygous for the allele).
5070 Is the progenitor strain (*y[1] scd[1] ras[1] v[1] f[1]*).
ORE and SAM are two standard lab wild types (Oregon-R and Samarkand). The ORE_WT and SAM_WT were pure lines (not crossed to 5070). All RAL* (where * is the number for the RAL line) represent sequenced DGRP strains used as a panel for genetic variation. FVW (FennVille Winery) is a large natural outbred population collected in Fennville Michigan.

slide_date - Date of dissections.

t1_teeth - count of sex comb teeth for normal sex comb on tarsal segment 1 (basi-tarsus).

t1_gaps - Number of gaps observed in primary sex comb.

t1_abn - Number of abnormal teeth in primary sex comb.

t1_extra_comb - Presence of a secondary sex comb on T1 (basi-tarsus).

t2_teeth - Number of sex comb teeth on the additional sex comb.
 that commonly occurs in scd[1] individuals on T2 (second tarsal segments).

## Figure 4

`compensation_4_30.csv` Data for allele purging experiment for the scd[1] allele, this data is analyzed in the script Figure4_plot_fitness.R for Figure 4. Columns for the data are as follows:

rep - replicate lineage for selection.

gen - generation.

scd - number of scd[1] males observed.

wt - number of wild type males observed.

unsure - individuals who may have been either scd[1] or wild type (scd[+]).

## Figure 5

`sex_comb_data_uptogen24.csv` contains the data for the compensatory evolution experiment. The script analyzing this experiment is Figure5_plot_main_EE.R reproducing figure 5. Please also note that the simulations for power analyses are also included in this script (supplemental figures).

Variables in this data set are:

rep - replicate lineage (NOT the same lineages as for the purging experiment).

treatment - experimental treatments for the evolution experiment (LV - low genetic variation, HSS - high sexual selection, LSS - low sexual selection, WTC - control lineages for lab domestication that are wild type for scd).

gen - generation.

id - individual.

leg - left or right leg.

primary_teeth - number of sex comb teeth in primary sex comb.

primary_gaps - presence of gaps in primary sex comb.

primary_partial	- ectopic sex comb on T1.

secondary_teeth	- Number of sex comb teeth on T2.

secondary_gaps - Gaps in secondary sex sex_comb_data_uptogen24	 
