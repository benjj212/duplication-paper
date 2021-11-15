list_accession <- c(1254,5856,6021,6909,9412,9470)

mar_set <- c(2,2,4)
for(k in c(1:length(list_accession))){
png(paste("/home/GMI/benjamin.jaegle/Documents_GMI/001_DATA/101_PAPER_DUPLICATION/EXAMPLE_AT1G20400/",list_accession[k],".methylation.png",sep=""), width = 1000, height = 1000)
par(mfrow=c(3,1))
list_context <- c("CG", "CHG", "CHH")
for(j in c(1:3)){
setwd("/home/GMI/benjamin.jaegle/Documents_GMI/001_DATA/101_PAPER_DUPLICATION/METHYLATION_DATA/")
dir()
accession <- list_accession[k]
context <- list_context[j]
GENE <- "AT1G20400"
### methylation data
X <- read.table(paste(accession,"_10c.",context,".bg", sep=""))

### Blast data 
setwd(paste("/home/GMI/benjamin.jaegle/Documents_GMI/001_DATA/101_PAPER_DUPLICATION/BLAST/",accession,sep=""))
BLAST <- read.table(file = paste(GENE,".fasta.",accession,".70.txt", sep=""))

selected_BLAST <- BLAST[which(BLAST[,6]-BLAST[,5]>1000),]
par(cex=2, mar =c(mar_set[j],4,1,1))
plot(0, type = "n", ylim=c(0,1), xlim=c(-2000,6000), main="", xlab="position around the insertion in bp", ylab="Methylation level")
colors_lines <- c("darkblue", "darkgreen", "darkred", "orange")
for(i in c(1:length(selected_BLAST[,1]))){
  START <- min(selected_BLAST[i,c(7,8)])-2000
  END <- max(selected_BLAST[i,c(7,8)])+2000
  CHR <- selected_BLAST[i,9]
  X_selected <- X[which(X[,2]>START & X[,3]<END & X[,1]==as.character(CHR)),]
  points(X_selected[,2]-min(selected_BLAST[i,c(7,8)]), X_selected[,4], type = "l", col=colors_lines[i], lwd=2)
}
}
dev.off()
}

X[c(1:5),]
dev.off()
