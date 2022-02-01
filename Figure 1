### installing package circlize
install.packages("circlize")
library("circlize")

### setting path 
Path_to_Rdata <- "/groups/nordborg/user/benjamin.jaegle/Documents/001_DATA/101_PAPER_DUPLICATION/GITHUB/FIGURE1/Figure1.RData"
Path_to_plot <- "/groups/nordborg/user/benjamin.jaegle/Documents/007_MANUSCRIPT/004_DUPLICATION/FIGURES/FIGURE1/circular_plot_gene_te_100mb_no_overlap.png"
Path_to_gff <- "/groups/nordborg/user/benjamin.jaegle/Documents/001_DATA/101_PAPER_DUPLICATION/GITHUB/FIGURE1/Araport11_GFF3_genes_transposons.201606.gff"
Path_for_temp <- "/groups/nordborg/user/benjamin.jaegle/Documents/007_MANUSCRIPT/004_DUPLICATION/FIGURES/FIGURE1/cytoband.txt"

### Loading 
load(Path_to_Rdata)

### circus plot 

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
araport <- gffRead(Path_to_gff)
araport_gene <- araport[which(araport$feature=="gene"),]
araport_TE <- araport[grep("transpos", araport$feature),]

## setting up chromosome coordinates
xlim=cbind(0, c(34964571,22037565,25499034,20862711,31270811))
cytoband.arabidopsis <- rbind(
  c("Chr1",0,32000000, "a", "gneg"),
  c("Chr2",0,20000000, "b", "gneg"),
  c("Chr3",0,24000000, "c", "gneg"),
  c("Chr4",0,20000000, "d", "gneg"),
  c("Chr5",0,28000000, "e", "gneg"))
write.table(cytoband.arabidopsis, Path_for_temp, sep="\t", row.names = F, col.names = F, quote = F)
cytoband.arabidopsis = read.table(Path_for_temp, colClasses = c("character", "numeric","numeric", "character", "character"), sep = "\t")
circos.par("clock.wise"=TRUE,start.degree=90)
circos.initializeWithIdeogram(cytoband.arabidopsis)
circos.clear()


window=100000 
by=100000 
list_threshold <- c(0,0.01,0.02,0.05,0.1,1)
centro_start <- c(14364752, 3602775, 12674550, 2919690, 11668616)
centro_end   <- c(15750321, 3735247, 13674767, 4011692, 12082583)
list_color <- c("blue","darkred","orange", "yellow", "white")
list_threshold <- c(0,1)

### getting pseudoSNPs density
out_all <- c()
#setting up the threshold of the pseudo-SNPs frequency: here between 0 and 1
  t=1
  threshold <- list_threshold[t]
  CHR_threshold <- CHR_all[which(CHR_all[,11]>threshold & CHR_all[,11]<list_threshold[t+1]),]
  win <- seq(0,35000000, by=by)
  # out is the overlapping sliding window data
  out <- as.data.frame(cbind(rep(win,5),rep(1:5,each=length(win)))) ; colnames(out) <- c('pos','chr') ; out$pi <- NA ; out$freq_hete <- NA 
  for (p in 1:nrow(out)){
    tmp <- CHR_threshold[CHR_threshold[,'CHR']==paste("Chr",out[p,'chr'],sep="") & CHR_threshold[,'POS']>(out[p,'pos'] - window/2) & CHR_threshold[,'POS']<(out[p,'pos'] + window/2),]
    out[p,'num_hete'] <-          length(tmp[,'CHR'])
  }
  out_all <- c(out_all, list(out))

#### getting gene density
out_gene <- as.data.frame(cbind(rep(win,5),rep(1:5,each=length(win)))) ; colnames(out_gene) <- c('pos','chr') ; out_gene$pi <- NA ; out_gene$freq_hete <- NA 
GENE_all <- cbind(araport_gene[,1],araport_gene[,4]+(araport_gene[,5]-araport_gene[,4]))
colnames(GENE_all) <- c("CHR", "POS")
for (p in 1:nrow(out_gene)){
  tmp <- GENE_all[GENE_all[,'CHR']==paste("Chr",out_gene[p,'chr'],sep="") & as.numeric(GENE_all[,'POS'])>(out_gene[p,'pos'] - window/2) & as.numeric(GENE_all[,'POS'])<(out_gene[p,'pos'] + window/2),]
  if(length(tmp)==2){
    out_gene[p,'num_gene'] <- 1
  }else{
    out_gene[p,'num_gene'] <- length(tmp[,'CHR'])
  }
}
out_all <- c(out_all, list(out_gene))


#### getting TE density
out_TE <- as.data.frame(cbind(rep(win,5),rep(1:5,each=length(win)))) ; colnames(out_TE) <- c('pos','chr') ; out_TE$pi <- NA ; out_TE$freq_hete <- NA
GENE_all <- cbind(araport_TE[,1],araport_TE[,4]+(araport_TE[,5]-araport_TE[,4]))
colnames(GENE_all) <- c("CHR", "POS")
for (p in 1:nrow(out_TE)){
  tmp <- GENE_all[GENE_all[,'CHR']==paste("Chr",out_TE[p,'chr'],sep="") & as.numeric(GENE_all[,'POS'])>(out_TE[p,'pos'] - window/2) & as.numeric(GENE_all[,'POS'])<(out_TE[p,'pos'] + window/2),]
  if(length(tmp)==2){
    out_gene[p,'num_gene'] <- 1
  }else{
    out_TE[p,'num_gene'] <-          length(tmp[,'CHR'])
  }
}
out_all <- c(out_all, list(out_TE))

### getting PI
out_PI <- as.data.frame(cbind(rep(win,5),rep(1:5,each=length(win)))) ; colnames(out) <- c('pos','chr') ; out$pi <- NA ; out$freq_hete <- NA

### plotting circulize plot
png(Path_to_plot, width=1000, height=1000)

circos.par("clock.wise"=TRUE,start.degree=90)
circos.genomicInitialize(cytoband.arabidopsis, sector.names = NULL, major.by = NULL,
                         plotType = c("axis", "labels"), tickLabelsStartFromZero = TRUE,
                         axis.labels.cex = 1, labels.cex = 1.5,
                         track.height = NULL)

for(t in c(1,2,3)){
  out_threshold <- out_all[[t]]
  circos.track(factors=NULL, ylim=c(1,max(out_threshold[,5])),track.height=0.12)
  for(i in c(1:5)){
    ### centromere
    circos.rect(xleft = centro_start[i], xright = centro_end[i], ybottom = 1, ytop = max(out_threshold[,5]), col="lightgrey", sector.index = paste("Chr",i,sep=""))
    ### data 
    circos.lines(out_threshold[which(out_threshold[,2]==i),1], out_threshold[which(out_threshold[,2]==i),5],
                 sector.index = paste("Chr",i,sep=""), type="l", area=T, track.index = t+1, col=list_color[t])
  }
}

dev.off()
