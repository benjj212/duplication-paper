# duplication-paper
all script analysis of the paper "Extensive gene duplication in Arabidopsis revealed by pseudo-heterozygosity"

### Basic SNP analysis:
Original VCF matrix 
https://doi.org/10.5281/zenodo.6025134

### Running the GWAS:

Using pygwas https://github.com/timeu/PyGWAS 

##### example:

pygwas run phenotype.csv -a kw -g 1001genomes_genotype -o "file.csv"


## Circular plot: 
Getting the RData from 

https://doi.org/10.1101/2021.11.15.468652

script: Circular_plot.R

## GWAS analysis: 
script: Meta_GWAS_analysis.R

## Methylation analysis
The methylation status were generated using  and can be found in the additionnal file 1-8 of the publication

script: meythylation comparison_ref_vs_pacbio.R

