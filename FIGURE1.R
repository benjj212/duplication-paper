load("/net/gmi.oeaw.ac.at/nordborg/user/benjamin.jaegle/Documents/002_SCRIPT/001_R/ashley script/Pi_Div_exp_500k_window_100k_step.RData")
load("/net/gmi.oeaw.ac.at/nordborg/user/benjamin.jaegle/Desktop/001_DATA/015_duplication/all_duplication_X_ALL_10.RData")
setwd("/net/gmi.oeaw.ac.at/nordborg/user/benjamin.jaegle/Documents/001_DATA/015_duplication/012_HOMO_HETE/")
HWE <- read.table("out_new.hwe", header=T, sep="\t")
threshold <- c(0,10,20,50,100)

for(p in c(1:5)){

HWE[c(1:10),]
HWE_NEW <- HWE[which(HWE[,4]>0),]
HWE_NEW <- cbind(HWE_NEW, HWE_NEW[,4]/(HWE_NEW[,5]+HWE_NEW[,3]+ HWE_NEW[,4]))
CHR_all <- c() 
for(i in c(1:5)){
  CHR <- HWE_NEW[which(HWE_NEW[,1]==paste("Chr",i,sep="")),]
  CHR_all <- rbind(CHR_all, CHR)
  print(i)
}
}
#### start here ###

load("/groups/nordborg/user/benjamin.jaegle/Documents/007_MANUSCRIPT/004_DUPLICATION/FIGURE1_raw.RData")

length()
length(which(HWE_NEW[,4]==1))
### 1258451 singleton
length(which(HWE_NEW[,11]>0.02 & HWE_NEW[,11]<0.05))



#####    A   #### histogram with the frequency of heterozygous SNPs
X11()
png("/net/gmi.oeaw.ac.at/nordborg/user/benjamin.jaegle/Documents/007_MANUSCRIPT/004_DUPLICATION/FIGURES/FIGURE1/FIGURE_1_A.png", width=1000, height = 1000)
par(fig=c(0,1,0,1), cex=1.5)
hist(as.numeric(CHR_all[,11]), breaks = 1000, xlab="frequency of heterozygots in the population", main="", ylab="number of SNPs")
dev.off()
par(fig=c(0.5,0.95,0.5,0.95), new=T, mai=c(0,0,0.3,0.3))
X11()
png("/net/gmi.oeaw.ac.at/nordborg/user/benjamin.jaegle/Documents/007_MANUSCRIPT/004_DUPLICATION/FIGURES/FIGURE1/FIGURE_1_A_freq_zoom_no_col.png", width=1000, height = 1000)
par(cex=1.5)
hist(as.numeric(CHR_all[,11]), breaks = 1000, xlim=c(0, 0.1),ylab="number of SNPs", xlab="frequency of heterozygots in the population", axes = T, main="")
#rect(xleft = list_threshold[c(1:5)], xright = list_threshold[c(2:6)], ybottom = 1170000, ytop = 1200000, col=c("darkred", "darkgreen","orange", "yellow", "white"))
dev.off()
#####    B  #### genome wide distribution of heterozygous SNPs (windows)
X11()
#png("/net/gmi.oeaw.ac.at/nordborg/user/benjamin.jaegle/Documents/007_MANUSCRIPT/004_DUPLICATION/FIGURE_1_B.png", width=500, height = 500)
#layout(matrix(c(1,1,2,1,1,3), 2, 3, byrow = TRUE))

cum_chr_length <- c(0,34964571, 34964571+22037565, 34964571+22037565+25499034, 34964571+22037565+25499034+20862711, 140000000)
Chr_length <- c(34964571, 22037565, 25499034, 20862711, 31270811)
centro_start <- c(14364752, 3602775, 12674550, 2919690, 11668616)
centro_end   <- c(15750321, 3735247, 13674767, 4011692, 12082583)


cen <- c(15088487,3607927,13591500,3956019,11742255)

### normalization of the coordinates
HWE_NEW$relpos <- NA ; for (i in 1:5){HWE_NEW[HWE_NEW$CHR == paste("Chr",i,sep=""),'relpos'] <- HWE_NEW[HWE_NEW$CHR == paste("Chr",i,sep=""),'POS']-cen[i]}
HWE$relpos <- NA ; for (i in 1:5){HWE[HWE$CHR == paste("Chr",i,sep=""),'relpos'] <- HWE[HWE$CHR == paste("Chr",i,sep=""),'POS']-cen[i]}
HWE$percent_hete <- HWE[,4]/(HWE[,3]+HWE[,5])
offset=c(.985,.79,.61,.42,.23) ; offset2=offset-.04 ; W=0.09; WW=0.15 ; x1=-15300000 ; x2=16400000
par(xpd=F, mar=c(2,3,.5,3) , fig=c(0,1,0,1), new=F,mgp=c(4,0.5,0))
plot(NA, ylim=c(0,0.04),xlim=c(x1,x2),xaxs='i',xaxt='n',yaxt='n') ; axis(1,at=c(-15,-10,-5,0,5,10,15)*10^6,c('-15','-10','-5','CEN','5','10','15')) ; abline(v=1,col='#000080',lwd=2)
mtext('Heterozygotes SNPs density',2,2)
threshold <- c(0,0.01,0.05,0.1,0.25)
colors_list <- c("darkred", "darkgreen", "orange", "yellow", "white")

window=500000 # kb
by=100000 # kb
win <- seq(-16000000,16100000, by=by)
# with all SNPs
out_all <- as.data.frame(cbind(rep(win,5),rep(1:5,each=length(win)))) ; colnames(out_all) <- c('pos','chr')
for (p in 1:nrow(out_all)){
  out[100,]
  tmp <- HWE[HWE[,'CHR']==paste("Chr",out_all[p,'chr'],sep="") & HWE[,'relpos']>(out_all[p,'pos'] - window/2) & HWE[,'relpos']<(out_all[p,'pos'] + window/2),]
  out_all[p,'um_hete'] <-          length(tmp[,'CHR'])
}
out_new_a <- out_all[-which(out_all[,3]<=0),]
for (i in 1:5){
  par(xpd=F, mar=c(0,3,0,3) , fig=c(0,1,offset[i]-W,offset[i]), new=TRUE)
  out_plot <- out_new_a[which(out_new_a[,2]==i),]
  one <- c(out_plot[1,1]-1, as.character(out_plot[1,2]), 0)
  two <- c(out_plot[length(out_plot[,1]),1]+1, as.character(out_plot[1,2]), 0)
  out_plot <- rbind(one, out_plot, two)
  plot(xlim=c(x1,x2),xaxs='i',xaxt='n',yaxt='n',bty='n', out_plot[out_plot$chr==i,'um_hete']  ~out_plot[out_plot$chr==i,'pos'] ,col='#000000', type='l', ylim=c(0,5000))
  polygon( out_plot[out_plot$chr==i,'pos'],out_plot[out_plot$chr==i,'um_hete'], col="red")
}

### with threshold
for(q in c(1:5)){
  HWE_NEW <- HWE[which(HWE[,12]>threshold[q]),]
  HWE_NEW <- cbind(HWE_NEW, HWE_NEW[,4]/(HWE_NEW[,5]+HWE_NEW[,3]+ HWE_NEW[,4]))
  CHR_all <- c() 
  for(i in c(1:5)){
    CHR <- HWE_NEW[which(HWE_NEW[,1]==paste("Chr",i,sep="")),]
    CHR_all <- rbind(CHR_all, CHR)
    print(i)
  }


# out is the overlapping sliding window data
out <- as.data.frame(cbind(rep(win,5),rep(1:5,each=length(win)))) ; colnames(out) <- c('pos','chr') ; out$pi <- NA ; out$mean_div <- NA 
for (p in 1:nrow(out)){
  out[100,]
  tmp <- HWE_NEW[HWE_NEW[,'CHR']==paste("Chr",out[p,'chr'],sep="") & HWE_NEW[,'relpos']>(out[p,'pos'] - window/2) & HWE_NEW[,'relpos']<(out[p,'pos'] + window/2),]
  out[p,'um_hete'] <-          length(tmp[,'CHR'])
}

out_new <- out[-which(out[,5]<=0),]
for (i in 1:5){
  par(xpd=F, mar=c(0,3,0,3) , fig=c(0,1,offset[i]-W,offset[i]), new=TRUE)
  out_plot <- out_new[which(out_new[,2]==i),]
  one <- c(out_plot[1,1]-1, as.character(out_plot[1,c(2,3,4)]), 0)
  two <- c(out_plot[length(out_plot[,1]),1]+1, as.character(out_plot[1,c(2,3,4)]), 0)
  out_plot <- rbind(one, out_plot, two)
  plot(xlim=c(x1,x2),xaxs='i',xaxt='n',yaxt='n',bty='n', out_plot[out_plot$chr==i,'um_hete']  ~out_plot[out_plot$chr==i,'pos'] ,col='#000000', type='l',ylim=c(0,5000))
  polygon( out_plot[out_plot$chr==i,'pos'],log(as.numeric(out_plot[out_plot$chr==i,'um_hete'])), col=colors_list[q])
  print(max(log(as.numeric(out_plot[,'um_hete']))))
}
}
dev.off()
#### selected 26000 snp

window=500000 # kb
by=100000 # kb
win <- seq(-16000000,16100000, by=by)
# out is the overlapping sliding window data
out <- as.data.frame(cbind(rep(win,5),rep(1:5,each=length(win)))) ; colnames(out) <- c('pos','chr') ; out$pi <- NA ; out$mean_div <- NA 
for (p in 1:nrow(out)){
  out[100,]
  tmp <- CHR_all[CHR_all[,'CHR']==paste("Chr",out[p,'chr'],sep="") & CHR_all[,'relpos']>(out[p,'pos'] - window/2) & CHR_all[,'relpos']<(out[p,'pos'] + window/2),]
  out[p,'um_hete'] <-          length(tmp[,'CHR'])
}

out_new <- out[-which(out[,5]<=0),]
for (i in 1:5){
  par(xpd=F, mar=c(0,3,0,3) , fig=c(0,1,offset[i]-W,offset[i]), new=TRUE)
  out_plot <- out_new[which(out_new[,2]==i),]
  one <- c(out_plot[1,1]-1, as.character(out_plot[1,c(2,3,4)]), 0)
  two <- c(out_plot[length(out_plot[,1]),1]+1, as.character(out_plot[1,c(2,3,4)]), 0)
  out_plot <- rbind(one, out_plot, two)
  plot(xlim=c(x1,x2),xaxs='i',xaxt='n',yaxt='n',bty='n', as.numeric(out_plot[out_plot$chr==i,'um_hete'])  ~out_plot[out_plot$chr==i,'pos'] ,col='#000000', type='l')
  polygon( out_plot[out_plot$chr==i,'pos'],as.numeric(out_plot[out_plot$chr==i,'um_hete']), col=rgb(1,0,0, alpha = 0.3, maxColorValue = 1))
}



length(levels(as.factor(CHR_all[,2])))




CHR_all[c(1:10),]
HWE[c(1:10),]

ALL_match <- c()
for(chr in c(1:5)){
  chr=1
  CHR_vcf_1 <- CHR_all[which(CHR_all[,1]==paste("Chr",chr,sep="")),]
  CHR_vcf_2 <- HWE[which(HWE[,1]==paste("Chr",chr,sep="")),]
  TEMP <- cbind(CHR_vcf_1, CHR_vcf_2[match(CHR_vcf_1[,2], as.numeric(CHR_vcf_2[,2])),])
  TEMP[c(1:10),]
  ALL_match <- rbind(ALL_match, TEMP[-which(is.na(TEMP[,95])),])
}




# with all SNPs

out_all <- as.data.frame(cbind(rep(win,5),rep(1:5,each=length(win)))) ; colnames(out_all) <- c('pos','chr')
for (p in 1:nrow(out_all)){
  out[100,]
  tmp <- HWE[HWE[,'CHR']==paste("Chr",out_all[p,'chr'],sep="") & HWE[,'relpos']>(out_all[p,'pos'] - window/2) & HWE[,'relpos']<(out_all[p,'pos'] + window/2),]
  out_all[p,'um_hete'] <-          length(tmp[,'CHR'])
}

out_all[c(1:10),]

out_new_a <- out_all[-which(out_all[,3]<=0),]
for (i in 1:5){
  par(xpd=F, mar=c(0,3,0,3) , fig=c(0,1,offset[i]-W,offset[i]), new=TRUE)
  out_plot <- out_new_a[which(out_new_a[,2]==i),]
  one <- c(out_plot[1,1]-1, as.character(out_plot[1,2]), 0)
  two <- c(out_plot[length(out_plot[,1]),1]+1, as.character(out_plot[1,2]), 0)
  out_plot <- rbind(one, out_plot, two)
  plot(xlim=c(x1,x2),xaxs='i',xaxt='n',yaxt='n',bty='n', out_plot[out_plot$chr==i,'um_hete']  ~out_plot[out_plot$chr==i,'pos'] ,col='#000000', type='l')
  polygon( out_plot[out_plot$chr==i,'pos'],out_plot[out_plot$chr==i,'um_hete'], col="darkblue")
}

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

araport <- gffRead("/groups/nordborg/user/benjamin.jaegle/Documents/001_DATA/003_annotations/001_TAIR/Araport11_GFF3_genes_transposons.201606.gff")
araport_gene <- araport[which(araport$feature=="gene"),]
araport_TE <- araport[grep("transpos", araport$feature),]

levels(as.factor(araport$feature))




library("circlize", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
xlim=cbind(0, c(34964571,22037565,25499034,20862711,31270811))
cytoband.arabidopsis <- rbind(
  c("Chr1",0,32000000, "a", "gneg"),
  c("Chr2",0,20000000, "b", "gneg"),
  c("Chr3",0,24000000, "c", "gneg"),
  c("Chr4",0,20000000, "d", "gneg"),
  c("Chr5",0,28000000, "e", "gneg"))
write.table(cytoband.arabidopsis, "/groups/nordborg/user/benjamin.jaegle/Documents/007_MANUSCRIPT/004_DUPLICATION/FIGURES/FIGURE1/cytoband.txt", sep="\t", row.names = F, col.names = F, quote = F)
cytoband.arabidopsis = read.table("/groups/nordborg/user/benjamin.jaegle/Documents/007_MANUSCRIPT/004_DUPLICATION/FIGURES/FIGURE1/cytoband.txt", 
                                  colClasses = c("character", "numeric","numeric", "character", "character"), sep = "\t")

circos.par("clock.wise"=TRUE,start.degree=90)
circos.initializeWithIdeogram(cytoband.arabidopsis)

circos.clear()

#circos.par("clock.wise"=TRUE,start.degree=90)
#circos.initialize(factors=as.factor(c("Chr1","Chr2","Chr3","Chr4","Chr5")),xlim=xlim)

CHR_all[c(1:10),]

window=100000 # kb
by=100000 # kb
list_threshold <- c(0,0.01,0.02,0.05,0.1,1)
centro_start <- c(14364752, 3602775, 12674550, 2919690, 11668616)
centro_end   <- c(15750321, 3735247, 13674767, 4011692, 12082583)
list_color <- c("darkred", "darkgreen","orange", "yellow", "white")

list_threshold <- c(0,1)

out_all <- c()
for(t in c(1:2)){
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
}
#### getting gene density
out_gene <- as.data.frame(cbind(rep(win,5),rep(1:5,each=length(win)))) ; colnames(out_gene) <- c('pos','chr') ; out_gene$pi <- NA ; out_gene$freq_hete <- NA 
TAIR10_GENE[c(1:10),]

GENE_all <- cbind(araport_gene[,1],araport_gene[,4]+(araport_gene[,5]-araport_gene[,4]))
colnames(GENE_all) <- c("CHR", "POS")
for (p in 1:nrow(out_gene)){
  tmp <- GENE_all[GENE_all[,'CHR']==paste("Chr",out_gene[p,'chr'],sep="") & as.numeric(GENE_all[,'POS'])>(out_gene[p,'pos'] - window/2) & as.numeric(GENE_all[,'POS'])<(out_gene[p,'pos'] + window/2),]
  if(length(tmp)==2){
    out_gene[p,'num_gene'] <- 1
  }else{
  out_gene[p,'num_gene'] <-          length(tmp[,'CHR'])
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




png("/groups/nordborg/user/benjamin.jaegle/Documents/007_MANUSCRIPT/004_DUPLICATION/FIGURES/FIGURE1/circular_plot_gene_te_100mb_no_overlap.png", width=1000, height=1000)

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

circos.genomicInitialize(cytoband.arabidopsis, sector.names = NULL, major.by = NULL,
                         plotType = c("axis", "labels"), tickLabelsStartFromZero = TRUE,
                         axis.labels.cex = 1.2, labels.cex = 1.2,
                         track.height = NULL)




circos.initializeWithIdeogram(plotType = c("axis", "labels"))

circos.initializeWithIdeogram(cytoband.arabidopsis)

circos.par("clock.wise"=TRUE,start.degree=90)
circos.initialize(factors=as.factor(c("Chr1","Chr2","Chr3","Chr4","Chr5")),xlim=xlim,plotType = c("axis", "labels") )
text(0, 0, "plotType = c('axis', 'labels')", cex = 1)


circos.initializeWithIdeogram()
circos.info()
circos.clear()

### polymorphism data
syn[c(1:10),]
CHR_all[c(1:10),]
### check snp between hete and binary 

CHR_all_5 <- CHR_all[-which(CHR_all$HET<5),]
CHR_SNP_hete <- CHR_all_5[which(CHR_all_5[,1]=="Chr1"),]
CHR_SNP_bina <- syn[which(syn[,1]==1),]

61426/length(CHR_SNP_bina[,1])

### matching SNPs between hete and binary and vis versa
matching <- match(CHR_SNP_hete$POS, CHR_SNP_bina$pos)
matching_no_na <- matching[-which(is.na(matching))]

matching_3 <- match(CHR_SNP_bina$pos,CHR_SNP_hete$POS)
matching_3_no_na <- matching_3[-which(is.na(matching_3))]

### selecting matching SNPs
CHR_all_bina_hete <- CHR_all[matching_3_no_na,]
### checking / plot frequency
hist(CHR_all_bina_hete$HET, breaks=300, ylim=c(0,100))
hist(as.numeric(CHR_all[,4]), breaks = 300, xlab="frequency of heterozygots in the population", ylim=c(0,100))

### selecting match snp in Binary set 
SNP_bina_hete <- CHR_SNP_bina[matching_no_na,]
syn_2 <- syn[-match(SNP_bina_hete$pos, syn$pos),]

### removing genes detected to be duplicated
duplicated_genes <- read.table("/net/gmi.oeaw.ac.at/nordborg/user/benjamin.jaegle/Documents/001_DATA/015_duplication/ALL_DUPLICATED.txt", header=T, as.is=T)

for(j in c(1:length(duplicated_genes[,1]))){
  
  syn_3 <- syn[-which(syn$gene==duplicated_genes[j,1]),]
  
}



### plotting

window=500000 # kb
by=100000 # kb
win <- seq(0,35000000, by=by)
# out is the overlapping sliding window data
hete <- as.data.frame(cbind(rep(win,5),rep(1:5,each=length(win)))) ; colnames(hete) <- c('pos','chr') ; out$all <- NA ; out$ab_5 <- NA 
for (p in 1:nrow(hete)){
  tmp <- CHR_all[CHR_all[,'CHR']==paste("Chr", hete[p,'chr'], sep="") & CHR_all[,'POS']>(hete[p,'pos'] - window/2) & CHR_all[,'POS']<(hete[p,'pos'] + window/2),]
  tmp_1 <- above_5[above_5[,'CHR']==paste("Chr", hete[p,'chr'], sep="") & above_5[,'POS']>(hete[p,'pos'] - window/2) & above_5[,'POS']<(hete[p,'pos'] + window/2),]
  tmp_2 <- above_10[above_10[,'CHR']==paste("Chr", hete[p,'chr'], sep="") & above_10[,'POS']>(hete[p,'pos'] - window/2) & above_10[,'POS']<(hete[p,'pos'] + window/2),]
  tmp_3 <- syn_3[syn_3[,'chr']==hete[p,'chr'] & syn_3[,'pos']>(hete[p,'pos'] - window/2) & syn_3[,'pos']<(hete[p,'pos'] + window/2),]
  hete[p,'all'] <-          length(tmp[,1])
  hete[p,'ab_5'] <-          length(tmp_1[,1]) 
  hete[p,'ab_10'] <-          length(tmp_2[,1])
  hete[p,'pi'] <- mean(tmp_3[,'pi.syn.th'],na.rm=T)
  hete[p,'snp_dens'] <- length(tmp_3[,'pi.syn.th'])
  print(p)
}


hete_chr <- hete[which(hete$chr==1),]
plot(hete_chr$pos, hete_chr$all, type="l")
points(hete_chr$pos, hete_chr$ab_5, type="l")
points(hete_chr$pos, hete_chr$ab_10, type="l")
points(hete_chr$pos, hete_chr$pi*200000, type="l", col="red", cex=2)
points(hete_chr$pos, hete_chr$snp_dens, type="l", col="blue", cex=2)





  out[p,'pi_med'] <-    median(tmp[,'pi.syn.th']*1000)
  out[p,'mean_div'] <-    mean(tmp[,'diff.fix'])
  out[p,'med_div'] <-   median(tmp[,'diff.fix'])
  out[p,'sites_count']<-length(tmp[,'pi.syn.th'])
  out[p,'genes_count']<-length(unique(tmp[,'gene']))
  out[p,'mean_div'] <-    mean(tmp[,'diff.fix'])
  out[p,'mean_na'] <-     mean(tmp[,'no.calls.th'])
  
  tmp <- exp[exp[,'chr']==out[p,'chr'] & exp[,'relpos']>(out[p,'pos'] - window/2) & exp[,'relpos']<(out[p,'pos'] + window/2),]
  out[p,'col'] <- median(tmp[,'col'])
  out[p,'lyr'] <- median(tmp[,'lyr'])
  out[p,'cap'] <- median(tmp[,'cap'])
  out[p,'med'] <- median(tmp[,'med'])
  out[p,'mean'] <- mean(tmp[,'mean'])
  out[p,'off'] <- mean(tmp[,'off'])
  out[p,'var'] <- median(tmp[,'var'])
  out[p,'sd'] <- mean(tmp[,'sd'])
  
  print(p)
}
