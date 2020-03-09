library(ggplot2)

#======================
# Calculating GC ratio 
#======================
fasta <- read.csv("570.PPYXA570.1_contig.fasta",header = F)
WGS <- paste(fasta[-grep("^>", fasta$V1),], collapse = "")
rm(fasta)

# Divide by 10000bp 
WGSdivide <- function(data){
  a <- as.integer(nchar(data)/10000)+1
  splited <- c()
  for(i in 1:a){
    b <- substr(data,1+(i-1)*10000,i*10000)
    splited <-c(splited,b)
  }
  print(as.data.frame(splited))
}

WGSdivide_result <- WGSdivide(WGS)

# Calcuating GC ratio per 10000bp
require(dplyr)
require(stringr)

GCratio <- function(contigseq)
{
  ratio <- c()
  for(i in 1:nrow(contigseq))
  {
    a <- length(gregexpr("[C]|[G]|[c]|[g]",contigseq[i,1])[[1]])
    b <- str_length(contigseq[i,1])
    c <- a/b
    ratio <- c(ratio,c)
  }
  print(ratio)
}

wgsGCratio <- GCratio(WGSdivide_result)
wgsGCratio <- data.frame(wgsGCratio,1:length(wgsGCratio)*10000)
colnames(wgsGCratio) <- c("GCratio", "length")
max_GC<-round(max(wgsGCratio$GCratio)*100,1)
min_GC<-round(min(wgsGCratio$GCratio)*100,1)

gc_ratio_whole <- round(100*GCratio(as.data.frame(WGS)),1)

wgsGCratio$range <- ifelse(wgsGCratio$GCratio>=gc_ratio_whole/100, paste0(max_GC,"~",gc_ratio_whole,"%"), paste0(gc_ratio_whole,"~",min_GC,"%"))

#==========================
# contig length information
#==========================
contig

contig <- read.csv("contig info.csv")
contig$seg.cumulative.End <- cumsum(contig$seg.End)
contig$seg.cumulative.Start <- c(1, contig$seg.cumulative.End+1)[1:length(contig$seg.cumulative.End)]

contig$seg.cumulative.End


#==========================
# CDS, RNA etc information
#==========================
gene <- read.csv("gene information.csv")

gene$ContigNo <- gene$Contig.accession
gene$ContigNo <- as.integer(gsub("PPYXA570[.][0]{1,}","",gene$ContigNo))
gene$ContigNo

gene$wholeBegin <- NA
gene$wholeEnd <- NA

for(i in 1:nrow(contig)){
gene[which(gene$ContigNo==i),]$wholeBegin <- gene[which(gene$ContigNo==i),]$Begin + (contig$seg.cumulative.Start[i]-1)
gene[which(gene$ContigNo==i),]$wholeEnd <- gene[which(gene$ContigNo==i),]$End + (contig$seg.cumulative.Start[i]-1)
}

#=================
#visulaization
#================

ylabel <- c("Contigs", "+Strand", "-Strand", "GC ratio")
yloc <- c(.55, 0.425,0.325,0.05)
xloc <- c(6139999.5, 6139999.5, 6139999.5, 6139999.5)
yaxislabel <- data.frame(ylabel,yloc,xloc)

DNAsizebreaks <- data.frame(x_axis=0:12*500000)

ggplot()+
  geom_bar(data = wgsGCratio, aes(length, y=(GCratio-gc_ratio_whole/100)*2, fill=range), size=2, stat="identity",inherit.aes = FALSE)+
  geom_rect(data=contig, aes(xmin=seg.cumulative.Start, xmax= seg.cumulative.End, ymin=0.5, ymax=0.6), color="white", size=0.1)+
  geom_rect(data=subset(gene, Strand=="+"&Feature=="CDS"), aes(xmin=wholeBegin, xmax=wholeEnd, ymin=0.4, ymax=0.45), fill="darksalmon", alpha=0.8)+
  geom_rect(data=subset(gene, Strand=="-"&Feature=="CDS"), aes(xmin=wholeBegin, xmax=wholeEnd, ymin=0.3, ymax=0.35), fill="darksalmon", alpha=0.8)+
  geom_point(data=subset(gene, Strand=="+"&Feature!="CDS"), aes(x=(wholeBegin+wholeEnd)/2, y=0.425, color=Feature), size=0.8)+
  geom_point(data=subset(gene, Strand=="-"&Feature!="CDS"), aes(x=(wholeBegin+wholeEnd)/2, y=0.325, color=Feature), size=0.8)+
  geom_segment(data=DNAsizebreaks, aes(x=x_axis, xend=x_axis, y=0.75,yend=0.8))+
  geom_text(data = yaxislabel, aes(xloc, yloc, label = ylabel), inherit.aes = F, size = 2.8)+
  scale_x_continuous(limits = c(0,6300000), expand = c(0, 0),
                     breaks =c(0:12*500000),
                     labels = c("0.0M","0.5M","1.0M","1.5M","2.0M","2.5M","3.0M","3.5M","4.0M",
                                "4.5M","5.00M","5.5M","6.0M"))+
  scale_y_continuous(limits = c(-1,0.8), expand = c(0, 0))+
  scale_color_manual(values=c("purple1","limegreen"))+
  xlab("")+
  ylab("")+
  labs(fill="GC ratio", color="Feature")+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_line(size = 0),
        axis.line = element_line(size = 0),
        axis.ticks.x = element_line(color = "black",size=0.1))+
  coord_polar(start = pi / 23.5)

