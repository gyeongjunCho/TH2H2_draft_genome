library(ggplot2)
library(dplyr)
library(stringr)

###############################
# Figure 1 A
###############################
#==================================
#  GC ratio and Gc skew Calculating
#==================================
# download TH2H2 fasta file from NCBI
# If you are unable to download due to changes in the NCBI system, search directly from NCBI to get TH2H2 fasta file.
download.file(url = "https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/RP/DG/RPDG01/RPDG01.1.fsa_nt.gz",
              destfile= "RPDG01.1.fsa_nt.gz")
fasta <- read.csv(gzfile(description="RPDG01.1.fsa_nt.gz", open = "RPDG01.1.fsa_nt"), header = F)[,1]
fasta <- gsub("^>.{1,}",">", fasta)
WGS <- paste(fasta, collapse = "")

# check WGS bp size
nchar(WGS)

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
max_GCratio<-round(max(wgsGCratio$GCratio)*100,1)
min_GCratio<-round(min(wgsGCratio$GCratio)*100,1)
whole_gc_ratio <- round(100*GCratio(as.data.frame(WGS)),1)

wgsGCratio$range <- ifelse(wgsGCratio$GCratio>=whole_gc_ratio/100, paste0(max_GCratio," ~ ",whole_gc_ratio," %"), paste0(whole_gc_ratio," ~ ",min_GCratio," %"))



# Calcuating GC skew per 10000bp
GCskew <- function(contigseq)
{
  gcskew<- c()
  for(i in 1:nrow(contigseq))
  {
    g <- length(gregexpr("[G]|[g]",contigseq[i,1])[[1]])
    c <- length(gregexpr("[C]|[c]",contigseq[i,1])[[1]])
    a <- (g-c)/(g+c)
    gcskew <- c(gcskew, a)
  }
  print(gcskew)
}

wgsGCskew <- GCskew(WGSdivide_result)
wgsGCskew <- data.frame(wgsGCskew,1:length(wgsGCskew)*10000)
colnames(wgsGCskew) <- c("GCskew", "length")
max_GCskew<-round(max(wgsGCskew$GCskew),1)
min_GCskew<-round(min(wgsGCskew$GCskew),1)
max_GCskew
whole_gc_skew <- GCskew(as.data.frame(WGS))
whole_gc_skew <- round(whole_gc_skew,3)
wgsGCskew$range <- ifelse(wgsGCskew$GCskew>=whole_gc_skew, paste0(max_GCskew," ~ ",whole_gc_skew," bp/bp"), paste0(whole_gc_skew," ~ ",min_GCskew," bp/bp"))

d#==========================
# contig length information
#==========================
contig <- read.csv("row data/TH2H2 contig info.csv")
contig$seg.cumulative.End <- cumsum(contig$seg.End)
contig$seg.cumulative.Start <- c(1, contig$seg.cumulative.End+1)[1:length(contig$seg.cumulative.End)]

contig$seg.cumulative.End


#==========================
# CDS, RNA etc information
#==========================
gene <- read.csv("row data/gene information.csv") # CLgenomics 1.55 version results
gene$ContigNo <- gene$Contig.accession
gene$ContigNo <- as.integer(gsub("Contig_[0]{0,}","",gene$ContigNo))

gene$wholeBegin <- NA
gene$wholeEnd <- NA
length(contig$seg.cumulative.Start)
for(i in 1:nrow(contig)){
  gene[which(gene$ContigNo==i),]$wholeBegin <- gene[which(gene$ContigNo==i),]$Begin + (contig$seg.cumulative.Start[i]-1)
  gene[which(gene$ContigNo==i),]$wholeEnd <- gene[which(gene$ContigNo==i),]$End + (contig$seg.cumulative.Start[i]-1)
}

#=================
# Visulaization
#================

ylabel <- c("Contigs", "+Strand", "-Strand", "GC ratio", "GC skew")
yloc <- c(.55, 0.425,0.325,whole_gc_ratio/100-0.28,-0.08)
xloc <- rep(6139999.5, 5)
yaxislabel <- data.frame(ylabel,yloc,xloc)

plotDNAsizebreaks <- data.frame(x_axis=0:12*500000,
                                labels = c("0.0Mbp","0.5Mbp","1.0Mbp","1.5Mbp","2.0Mbp","2.5Mbp","3.0Mbp","3.5Mbp","4.0Mbp",
                                           "4.5Mbp","5.00Mbp","5.5Mbp","6.0Mbp"))
whole_gc_ratio/100
whole_gc_skew

wgsGCratio$GCratio+whole_gc_ratio/100

pesudogenome_map<- ggplot()+
  geom_segment(data = wgsGCratio, aes(x=as.integer(length), xend=as.integer(length), y=whole_gc_ratio/100-0.3, yend=GCratio-0.3, color=range), size=0.5, stat="identity")+
  geom_segment(data = wgsGCskew, aes(x=as.integer(length), xend=as.integer(length), y=(whole_gc_skew/2)-0.08,yend=GCskew/2-0.08, color=range), size=0.5, stat="identity")+
  geom_rect(data=contig, aes(xmin=seg.cumulative.Start, xmax= seg.cumulative.End, ymin=0.5, ymax=0.6, fill=ifelse(NO%%2==0,"Even","Odd"), linetype=ifelse(NO%%2==0,"Even","Odd")), color="white", size=0.1)+
  geom_rect(data=subset(gene, Strand=="+"&Feature=="CDS"), aes(xmin=as.integer(wholeBegin), xmax=as.integer(wholeEnd), ymin=0.4, ymax=0.45), fill="grey", color="grey",size=0.2, alpha=0.8)+
  geom_rect(data=subset(gene, Strand=="-"&Feature=="CDS"), aes(xmin=as.integer(wholeBegin), xmax=as.integer(wholeEnd), ymin=0.3, ymax=0.35), fill="grey",  color="grey",size=0.2, alpha=0.8)+
  geom_point(data=subset(gene, Strand=="+"&Feature!="CDS"), aes(x=(wholeBegin+wholeEnd)/2, y=0.425, color=Feature, alpha=Feature), size=2, shape=1)+
  geom_point(data=subset(gene, Strand=="-"&Feature!="CDS"), aes(x=(wholeBegin+wholeEnd)/2, y=0.325, color=Feature), size=2, shape=1)+
  geom_segment(data=plotDNAsizebreaks, aes(x=x_axis, xend=x_axis, y=0.62,yend=0.67))+
  geom_text(data=plotDNAsizebreaks, aes(x=x_axis, y=0.8, label=labels), size=3)+
  geom_text(data=yaxislabel, aes(xloc, yloc, label = ylabel), inherit.aes = F, size = 2.8)+
  annotate(geom="text", x=6139999.5, y=-0.6, label="  Paenibacillus polymyxa", fontface="italic", size=6.5)+
  annotate(geom="text", x=6139999.5, y=-0.75, label=" TH2H2", size=6.5)+
  scale_y_continuous(limits = c(-1,0.85), expand = c(0, 0))+
  scale_x_continuous(limits = c(0,6300000), expand = c(0, 0))+
  scale_color_manual(values=c("tomato","steelblue","yellowgreen","orange","deeppink","purple1"), guide=FALSE)+
  scale_alpha_manual(values = c(1,1), 
                     guide=guide_legend(order=1,
                                        override.aes = list(color=c("deeppink","purple1"),size=2, stroke=0.8)))+
  scale_linetype_manual(values = c(0, 0),
                        name = "GC ratio",
                        labels=c(paste0(whole_gc_ratio," ~ ",min_GCratio," %"),
                                 paste0(max_GCratio," ~ ",whole_gc_ratio," %")),
                        guide=guide_legend(order=2,
                                           override.aes = list(fill = c("yellowgreen","orange"))))+
  scale_fill_manual(values=c("gray25", "gray50"),
                    name="GC skew",
                    labels=c(paste0(whole_gc_skew," ~ ",min_GCskew," bp/bp"),
                             paste0(max_GCskew," ~ ",whole_gc_skew," bp/bp")),
                    guide=guide_legend(order=3,
                                       override.aes = list(fill = c("tomato", "steelblue"))))+
  theme_void()+
  ggtitle("A")+
  theme(plot.title.position = "plot",
        plot.title = element_text(margin = margin(t = 10, b = -100), size=40, hjust=0.1, face = "bold"),
        legend.position = c(0.5, 0.42),
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.text = element_text(size=9),
        legend.title.align = 0.5,
        legend.key.width = unit(3,units = "pt"),
        plot.margin = unit(c(0,-100,-80,10), "pt"))+
  coord_polar(start = pi / 23.5)
pesudogenome_map

###############################
# Figure 1 B
###############################
eggNOG <- read.csv("TH2H2 EggNOG.csv")
eggNOG_plot <- ggplot(eggNOG, aes(x=EggNOG_category, y=counts, label=counts, fill=description))+
  geom_bar(stat="identity")+
  geom_text(vjust=-1)+
  scale_y_continuous(limits=c(0, 1750), expand=c(0,0))+
  scale_fill_manual(values=rep("grey20", nrow(eggNOG)), name="Category discription",
                    guide=guide_legend(ncol=2,title.position = "top", override.aes = list(fill="white")))+
  xlab("EggNOG category")+
  ylab("")+
  theme_classic()+
  ggtitle("B")+
  theme(legend.position = "bottom",
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        legend.key.width = unit(0,"pt"),
        plot.title.position = "plot",
        plot.title = element_text(size=40, face = "bold"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15),
        plot.margin = unit(c(10,10,10,10), "pt"))

eggNOG_plot

###############################
# Merge to figure 1
##########################
library(gridExtra)
Figure1 <- grid.arrange(pesudogenome_map,eggNOG_plot, nrow=1)

ggsave(file="Figure1.svg", plot=Figure1, width=21, height=9)

