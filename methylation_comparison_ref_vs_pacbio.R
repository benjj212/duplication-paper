### reading araport annotation
gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",  
                                "integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}

TAIR10 <- gffRead("/path_to_annotation/Araport11_GFF3_genes_transposons.201606.gff")

TAIR10_gene <- TAIR10[which(TAIR10$feature=="gene"),]
TAIR10_gene$attributes <- substr(TAIR10_gene$attributes,4,12)

#### table from paper "Extensive gene duplication in Arabidopsis revealed by pseudo-heterozygosity" additional file 1
CG_CHG <- read.table("/path/table.txt", header = T, row.names = 1)
### aplying the threshold for CG and CHG
CG <- as.matrix(CG_CHG[,c(1:7)])
CHG <- as.matrix(CG_CHG[,c(8:14)])
### CG
CG[which(CG>=0.05)] <- 1
CG[which(CG<0.05)] <- 0
## CHG
CHG[which(CHG>=0.03)] <- 1
CHG[which(CHG<0.03)] <- 0

CG_CHG <- cbind(CG,CHG)

### data mapped to mapped to pacbio
#### tables from paper "Extensive gene duplication in Arabidopsis revealed by pseudo-heterozygosity" additional file 2-8
setwd("/path_to_tables/")
methylation_calc <- dir()
### defining the accession based on the additional file order
pacbio_accessions <- c("1254", "5856", "6021","6024", "6909", "9412", "9470")

missing_per_pacbio <- c()
exactly_one <- c()
more_than_one <- c()
methylation_status_all <- c()
length_check_all <- c()
for(p in c(2:8)){
  genome_name <- methylation_calc[grep(p, methylation_calc)]
  methylation_calc_data <- read.table(genome_name, sep=",", header=T)
  ### removing obsolete genes
  length_check <- c()
  methylation_status <- c()
  for(h in c(1:length(both_CG_CHG[,1]))){
    length_check <- c(length_check,length(methylation_calc_data[grep(rownames(both_CG_CHG)[h], methylation_calc_data[,1]),1]))
    methylation_status <- c(methylation_status,list(methylation_calc_data[grep(rownames(both_CG_CHG)[h], methylation_calc_data[,1]),]))
  }
  methylation_status_all <- c(methylation_status_all, list(methylation_status))
  length_check_all <- c(length_check_all, list(length_check))
  missing_per_pacbio <- c(missing_per_pacbio,length(which(length_check==0)))
  exactly_one <- c(exactly_one, length(which(length_check==1)))
  more_than_one <- c(more_than_one, length(which(length_check>1)))
}

### single copy comparison
pacbio_CHG <- c()
pacbio_CG <- c()
missing_pacbio <- c()

not_in_list <- c()
for(p in c(1:length(pacbio_accessions))){
  gene_list <- methylation_status_all[[p]]
  length_number <- length_check_all[[p]]
  indices <- which(length_number==1)
  CG_status <- c()
  CHG_status <- c()
  missing_data <- c()
  chr_match <- c()
  for(g in indices){
    gene_X <- gene_list[[g]]
      if(is.na(gene_X[1,2])){
        missing_data <- c(missing_data, g)
      }else{
        ### applying threshold to data mapped to pacbio
        gene_X[which(gene_X[,2]>=0.05),2] <- 1
        gene_X[which(gene_X[,2]<0.05),2] <- 0
        gene_X[which(gene_X[,3]<0.03),3] <- 0
        gene_X[which(gene_X[,3]>=0.03),3] <- 1
        gene_X <- rbind(c("gene",both_CG_CHG[which(row.names(both_CG_CHG)==gene_X[1,1]),c(p, p+7)],0), gene_X)
        gene_X <- as.matrix(gene_X)
        gene_X[which(is.na(gene_X))] <- 0
        if(gene_X[1,2]==gene_X[2,2]){
          CG_status <- c(CG_status, "match")
        }else{
          CG_status <- c(CG_status, "no")
        }
        if(gene_X[1,3]==gene_X[2,3]){
          CHG_status <- c(CHG_status, "match")
        }else{
          CHG_status <- c(CHG_status, "no")
        }
      }
    }
  pacbio_CG <- c(pacbio_CG, length(which(CG_status=="no")))
  pacbio_CHG <- c(pacbio_CHG, length(which(CHG_status=="no")))
  missing_pacbio <- c(missing_pacbio, length(missing_data))
}

### multiple copies 
pacbio_CHG_m <- c()
pacbio_CG_m <- c()


for(p in c(1:length(pacbio_accessions))){
  gene_list <- methylation_status_all[[p]]
  length_number <- length_check_all[[p]]
  indices <- which(length_number>1)
  CG_status <- c()
  CHG_status <- c()
  for(g in indices){
    gene_X <- gene_list[[g]]
    gene_X[which(gene_X[,2]>0.05),2] <- 1
    gene_X[which(gene_X[,2]<0.05),2] <- 0
    gene_X[which(gene_X[,3]<0.03),3] <- 0
    gene_X[which(gene_X[,3]>0.03),3] <- 1
    epi_status <- both_CG_CHG[g,c(p, p+7)]
    if(length(which(gene_X[,2]==epi_status[1]))==length(gene_X[,1])){
      CG_status <- c(CG_status, "match")
    }else{
      CG_status <- c(CG_status, "no")
    }
    if(length(which(gene_X[,3]==epi_status[2]))==length(gene_X[,1])){
      CHG_status <- c(CHG_status, "match")
    }else{
      CHG_status <- c(CHG_status, "no")
    }
  }
  pacbio_CG_m <- c(pacbio_CG_m, length(which(CG_status=="no")))
  pacbio_CHG_m <- c(pacbio_CHG_m, length(which(CHG_status=="no")))
}

pacbio_CG_m/more_than_one
pacbio_CHG_m/more_than_one
