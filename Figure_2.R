### packages to install
library("plotrix")
library("nucleR")
### Setting up the path
path_to_file <- "/groups/nordborg/user/benjamin.jaegle/Documents/001_DATA/101_PAPER_DUPLICATION/ANALYSIS/017_plotting_GWAS_vs_SNPs/GWAS_plot.pdf"
path_to_csv <- "/groups/nordborg/user/benjamin.jaegle/Documents/001_DATA/015_duplication/013_GWAS/HETE/CSV_2/"
### loading the RData
load("/groups/nordborg/user/benjamin.jaegle/Documents/001_DATA/101_PAPER_DUPLICATION/ANALYSIS/017_plotting_GWAS_vs_SNPs/001_script/plotting_GWAS.RData")
##### ploting SNP vs pvalue
pdf(path_to_file, width=10, height = 10)
par(fig=c(0,1,0,1), mar=c(0,0,0,0))
plot(0,xlim=c(0, 140000000), ylim=c(0, 140000000), type="n", axes=F)
setwd(path_to_csv)
all <- dir()
chr_length <- c(0,34964571, 34964571+22037565, 34964571+22037565+25499034, 34964571+22037565+25499034+20862711)
csv_2 <- all[grep("_2_", all)]
SNP <- tools::file_path_sans_ext(csv_2)
for(i in c(1:length(SNP))){
  chr <- strsplit(SNP[i],"_")[[1]][1]
  pos <- as.numeric(strsplit(SNP[i],"_")[[1]][2])
  file <- read.table(paste(SNP[i],".csv",sep=""), header = T)

  ### selection based on fourrier transform
   all_peaks <- c()
   for(k in c(1:5)){
   file_CHR <- file[which(file$chromosomes==k),]
   if(length(file_CHR[,1])>1){
   extended <- cbind(k,c(1:100), 0.1, 0.1, 0.1,0.1)
   colnames(extended) <- colnames(file_CHR)
   file_CHR <- rbind(file_CHR, extended)
   htseq_raw <- as.vector(-log10(file_CHR$pvals))
   htseq_fft <- filterFFT(htseq_raw, pcKeepComp=0.002)
   peaks <- peakDetection(htseq_fft, threshold="10%", score=FALSE)
   all_peaks <- c(all_peaks,length(peaks))
   }
   }
   if(sum(all_peaks)<20){

  if(length(which(file$positions>pos-500 & file$positions<pos+500 & -log10(file$pvals)>10))>1){
  if(length(file[,1])>0){
    ###selecting based on the mafs
  file <- file[which(file$mafs>0.2),]
  if(length(file[,1])>0){
    ###selecting based on the -log10(pvalue)
  file <- file[which(-log10(file$pvals)>10),]
  file_100 <- file
  chr_x <- file_100
  if(length(which(all_CHR[,1]==chr & all_CHR[,2]==pos))==1){
  if(length(levels(as.factor(chr_x$pvals)))>2){
  points(matrix(pos+chr_length[as.numeric(substr(chr,4,5))],nrow=length(chr_x[,1]),ncol=1)[,1],
         chr_x$positions+chr_length[chr_x$chromosomes],
         pch=16, cex=0.5,
         col=color.scale(chr_x$pvals, extremes=c("red", "red")))
  }else{
    points(matrix(pos+chr_length[as.numeric(substr(chr,4,5))],nrow=length(chr_x[,1]),ncol=1)[,1],
           chr_x$positions+chr_length[chr_x$chromosomes],
           pch=16, cex=0.5,
           col="red")
  }
  print(c(i, length(chr_x[,1])))
  }else{
    if(length(levels(as.factor(chr_x$pvals)))>2){
      points(matrix(pos+chr_length[as.numeric(substr(chr,4,5))],nrow=length(chr_x[,1]),ncol=1)[,1],
             chr_x$positions+chr_length[chr_x$chromosomes],
             pch=16, cex=0.5,
             col=color.scale(chr_x$pvals, extremes=c("lightgrey", "grey")))
    }else{
      points(matrix(pos+chr_length[as.numeric(substr(chr,4,5))],nrow=length(chr_x[,1]),ncol=1)[,1],
             chr_x$positions+chr_length[chr_x$chromosomes],
             pch=16, cex=0.5,
             col="darkgrey")
    }
    print(c(i, length(chr_x[,1])))

  }
  }
}
}
}
}
abline(h=chr_length)
abline(v=11400000)
centro_start <- c(14364752, 3602775, 12674550, 2919690, 11668616)
centro_end   <- c(15750321, 3735247, 13674767, 4011692, 12082583)
rect(xleft = 0, xright = 140000000, ybottom = chr_length+centro_start, ytop = chr_length+centro_end, col=rgb(150,150,150,alpha=100, maxColorValue = 255),border = F)
dev.off()
