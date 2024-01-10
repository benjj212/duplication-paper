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

[https://doi.org/10.1101/2021.11.15.468652](https://doi.org/10.1101/2021.11.15.468652)

script: Circular_plot.R

## GWAS analysis: 
script: Meta_GWAS_analysis.R

## BLAST 

Blast software was download from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 

##### example:

/blast_path/bin/blastn -query /path_to_query/GENE.fasta -db /path_to_database/GENOME.fasta -task blastn -evalue 1 -out output.txt -perc_identity 70 -outfmt "7 qseqid qacc sacc evalue qstart qend sstart send sseqid" -penalty -2 -reward 1 -word_size 28 -gapopen 1 -gapextend 1


## Methylation analysis
The methylation status were generated as described in https://github.com/Gregor-Mendel-Institute/pyBsHap and can be found in the additionnal file 1-8 of the publication

script: meythylation comparison_ref_vs_pacbio.R

