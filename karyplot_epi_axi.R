load("gene_density_windows.Rdata")
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

#gene density plot with methylation
library(karyoploteR)

#load fasta and gff
almond.fasta<- readDNAStringSet("/fs/scratch/PAS1755/Will/EpiAxi/Genomes/Almond/Almond.fasta")
almond.gff<- import.gff("/fs/scratch/PAS1755/Will/EpiAxi/Visualization/In_DMRcaller/context_report/genomic.gff")


#adjust the RLE names in the GFF so they match with the Granges methylation object
seqlevels(almond.gff)<- sub("CM037988.1", "chr001", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037989.1", "chr002", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037990.1", "chr003", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037991.1", "chr004", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037992.1", "chr005", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037993.1", "chr006", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037994.1", "chr007", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037995.1", "chr008", seqlevels(almond.gff))

chr_sizes<-seqlengths(almond.fasta)
chr_sizes<-chr_sizes[1:8]

almond.genome.size<-toGRanges(data.frame(chr=c("chr001", "chr002","chr003","chr004","chr005","chr006","chr007","chr008"), start= rep(1,8), end=chr_sizes))

#fix the window column names

colnames(mcols(windows.100k)) <- make.names(colnames(mcols(windows.100k)))

colnames(mcols(windows.100k)) <- gsub("^([0-9])", "_\\1", colnames(mcols(windows.100k)))

#load  DMRs from the annotation calls

#CpG 8-160
load("/users/PAS1755/whazzard/DMRCGs160_annotation.RData")
DMRs_160_cpg<-as.data.frame( c(DMRsGenesCG1, DMRsGenesCG2, DMRsGenesCG3, DMRsGenesCG4, DMRsGenesCG5, DMRsGenesCG6, DMRsGenesCG7, DMRsGenesCG, DMRsGenesCGmit)) %>%
  mutate(context="CpG")
DMRs.160<- DMRs_160_cpg %>%
  filter(type == c("CDS"))%>%
  select(seqnames, start, end, strand)

DMRs.160.granges<-makeGRangesFromDataFrame(DMRs.160, keep.extra.columns = T)

#CpG 8-201
load("/users/PAS1755/whazzard/DMRCGs201_annotation.RData")
DMRs_201_cpg<-as.data.frame( c(DMRsGenesCG1, DMRsGenesCG2, DMRsGenesCG3, DMRsGenesCG4, DMRsGenesCG5, DMRsGenesCG6, DMRsGenesCG7, DMRsGenesCG, DMRsGenesCGmit)) %>%
  mutate(context="CpG")
DMRs.201<- DMRs_201_cpg %>%
  filter(type == c("CDS"))%>%
  select(seqnames, start, end, strand)


#CHG 8-160
load("/users/PAS1755/whazzard/DMRCHG_annotation.RData")
DMRs_160_CHG<-as.data.frame( c(DMRsGenesCG1, DMRsGenesCG2, DMRsGenesCG3, DMRsGenesCG4, DMRsGenesCG5, DMRsGenesCG6, DMRsGenesCG7, DMRsGenesCG)) %>%
  mutate(context="CHG")
DMRs.160.CHG<- DMRs_160_CHG %>%
  filter(type == c("CDS"))%>%
  select(seqnames, start, end, strand)

DMRs.160.CHG.granges<-makeGRangesFromDataFrame(DMRs.160.CHG, keep.extra.columns = T)

#CHG 8-201
load("/users/PAS1755/whazzard/DMRCHG201_annotation.RData")
DMRs_201_CHG<-as.data.frame( c(DMRsGenesCG1, DMRsGenesCG2, DMRsGenesCG3, DMRsGenesCG4, DMRsGenesCG5, DMRsGenesCG6, DMRsGenesCG7, DMRsGenesCG, DMRsGenesCGmit)) %>%
  mutate(context="CHG")
DMRs.201.CHG<- DMRs_201_CHG %>%
  filter(type == c("CDS"))%>%
  select(seqnames, start, end, strand)

DMRs.201.CHG.granges<-makeGRangesFromDataFrame(DMRs.201.CHG, keep.extra.columns = T)


#CHH 8-160
load("/users/PAS1755/whazzard/DMRCHH160_annotation.RData")
DMRs_160_CHH<-as.data.frame( c(DMRsGenesCG1, DMRsGenesCG2, DMRsGenesCG3, DMRsGenesCG4, DMRsGenesCG5, DMRsGenesCG6, DMRsGenesCG7, DMRsGenesCG, DMRsGenesCGmit)) %>%
  mutate(context="CHH")
DMRs.160.CHH<- DMRs_160_CHH %>%
  filter(type == c("CDS"))%>%
  select(seqnames, start, end, strand)

DMRs.160.CHH.granges<-makeGRangesFromDataFrame(DMRs.160.CHH, keep.extra.columns = T)

#CHH 8-201
load("/users/PAS1755/whazzard/DMRCHH201_annotation.RData")
DMRs_201_CHH<-as.data.frame( c(DMRsGenesCG1, DMRsGenesCG2, DMRsGenesCG3, DMRsGenesCG4, DMRsGenesCG5, DMRsGenesCG6, DMRsGenesCG7, DMRsGenesCG, DMRsGenesCGmit)) %>%
  mutate(context="CHH")
DMRs.201.CHH<- DMRs_201_CHH %>%
  filter(type == c("CDS"))%>%
  select(seqnames, start, end, strand)

DMRs.201.CHH.granges<-makeGRangesFromDataFrame(DMRs.201.CHH, keep.extra.columns = T)

#superplot 


#overall plot
kp <- plotKaryotype(plot.type=2, genome = almond.genome.size, cex = 0.6)  

kpDataBackground(kp, data.panel = 1, col="grey80")
kpDataBackground(kp, data.panel = 2, col="grey88")


kpAddLabels(kp, labels="UCD 8-160", cex=0.5, data.panel = 1)
kpAddLabels(kp, labels="UCD8-201", cex=0.5, data.panel = 2)


# Gene density as cytobands
kp <- kpPlotDensity(kp, data = almond.gff, window.size = 100000, r0 = 0, r1 = 1, col = "black", border = NA, data.panel = "ideogram")

# Add 8-160 epiaxi methylation for CpG on the top  
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X160_CpG_axi_Proportion, col="chocolate", pch=16, data.panel = 1)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X160_CpG_epi_Proportion, col="blue", pch=16, data.panel = 1, lty=2)

#Add 8-160 CHG
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X160_CHG_axi_Proportion, col="gold2", pch=16, data.panel = 1)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X160_CHG_epi_Proportion, col="dodgerblue", pch=16, data.panel = 1, lty=2)

#Add 8-160 CHH
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X160_CHH_axi_Proportion, col="goldenrod2", pch=16, data.panel = 1)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X160_CHH_epi_Proportion, col="cornflowerblue", pch=16, data.panel = 1, lty=2)

# Add 8-201 epiaxi methylation for CpG on the bottom  
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X201_CpG_axi_Proportion, col="chocolate", pch=16, data.panel = 2, r0=1, r1=0)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X201_CpG_epi_Proportion, col="blue", pch=16, data.panel = 2, r0=1, r1=0, lty=2)

#Add 201 CHG
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X201_CHG_axi_Proportion, col="gold2", pch=16, data.panel = 2, r0=1, r1=0)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X201_CHG_epi_Proportion, col="dodgerblue", pch=16, data.panel = 2, r0=1, r1=0, lty=2)

#Add 201 CHH
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X201_CHH_axi_Proportion, col="goldenrod2", pch=16, data.panel = 2, r0=1, r1=0)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X201_CHH_epi_Proportion, col="cornflowerblue", pch=16, data.panel = 2, r0=1, r1=0, lty=2)

#add the DMRs for gene areas

kpRect(kp, chr = as.character(seqnames(DMRs.160.granges)), x0 = start(DMRs.160.granges), x1 = end(DMRs.160.granges),
       y0 = 0, y1 = 1, col = "#FFA50080", border = "orange", lwd = 2, data.panel = 1)

kpRect(kp, chr = as.character(seqnames(DMRs.160.CHG.granges)), x0 = start(DMRs.160.CHG.granges), x1 = end(DMRs.160.CHG.granges),
       y0 = 0, y1 = 1, col = "#FFA50080", border = "sienna", lwd = 2, data.panel = 1)

kpRect(kp, chr = as.character(seqnames(DMRs.160.CHH.granges)), x0 = start(DMRs.160.CHH.granges), x1 = end(DMRs.160.CHH.granges),
       y0 = 0, y1 = 1, col = "#FFA50080", border = "orangered", lwd = 2, data.panel = 1)

kpRect(kp, chr = as.character(seqnames(DMRs.201.granges)), x0 = start(DMRs.201.granges), x1 = end(DMRs.201.granges),
       y0 = 0, y1 = 1, col = "#FFA50080", border = "orange", lwd = 2, data.panel = 2)

kpRect(kp, chr = as.character(seqnames(DMRs.201.CHG.granges)), x0 = start(DMRs.201.CHG.granges), x1 = end(DMRs.201.CHG.granges),
       y0 = 0, y1 = 1, col = "#FFA50080", border = "sienna", lwd = 2, data.panel = 2)

kpRect(kp, chr = as.character(seqnames(DMRs.201.CHH.granges)), x0 = start(DMRs.201.CHH.granges), x1 = end(DMRs.201.CHH.granges),
       y0 = 0, y1 = 1, col = "#FFA50080", border = "orangered", lwd = 2, data.panel = 2)





#legend
legend(x = "bottomright",
       lty = c(1, 2, 1, 2, 1, 2), 
       col = c("chocolate", "blue", "gold2", "dodgerblue", "goldenrod2", "cornflowerblue"),
       legend = c("CpG Canopy", "CpG Epicormic", "CHG Canopy", "CHG Epicormic", "CHH Canopy", "CHH Epicormic"))

legend(x = "bottomright",
       fill = c("orange", "sienna", "orangered"),
       legend = c("CpG DMR", "CHG DMR", "CHH DMR"))




#NOTE HERE: right now everything is a little too squished together, or when the panels are seperate it comprises the ability to discern between them, 
#even more sor than the panels right on top of each other. lets keep the first plot as is and try to reform a little. 
#let's try doing it as multiple panels

kp <- plotKaryotype(plot.type=2, genome = almond.genome.size, cex = 0.6)  
kpDataBackground(kp, data.panel = 1, col="lightblue1")
kpDataBackground(kp, data.panel = 2, col="grey84")

# Gene density as cytobands
kp <- kpPlotDensity(kp, data = almond.gff, window.size = 100000, r0 = 0, r1 = 1, col = "black", border = NA, data.panel = "ideogram")

# Add 8-160 epiaxi methylation for CpG on the top 

at<-autotrack(current.track=3, total.tracks = 3)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X160_CpG_axi_Proportion, r0=at$r0, r1=at$r1, col="chocolate", pch=16, data.panel = 1)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X160_CpG_epi_Proportion, r0=at$r0, r1=at$r1, col="blue", pch=16, data.panel = 1, lty=2)
#Add 8-160 CHG

at<-autotrack(current.track=2, total.tracks = 3)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X160_CHG_axi_Proportion, r0=at$r0, r1=at$r1, col="gold2", pch=16, data.panel = 1)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X160_CHG_epi_Proportion, r0=at$r0, r1=at$r1, col="dodgerblue", pch=16, data.panel = 1, lty=2)

#Add 8-160 CHH
at<-autotrack(current.track=1, total.tracks = 3)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X160_CHH_axi_Proportion,r0=at$r0, r1=at$r1,  col="goldenrod2", pch=16, data.panel = 1)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X160_CHH_epi_Proportion,r0=at$r0, r1=at$r1, col="cornflowerblue", pch=16, data.panel = 1, lty=2)

#Add 201 CpG
at<-autotrack(current.track=3, total.tracks = 3)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X201_CpG_axi_Proportion, col="chocolate", pch=16, data.panel = 2, r0=at$r0, r1=at$r1)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X201_CpG_epi_Proportion, col="blue", pch=16, data.panel = 2, r0=at$r0, r1=at$r1, lty=2)

#Add 201 CHG
at<-autotrack(current.track=2, total.tracks = 3)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X201_CHG_axi_Proportion, col="gold2", pch=16, data.panel = 2, r0=at$r0, r1=at$r1)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X201_CHG_epi_Proportion, col="dodgerblue", pch=16, data.panel = 2, r0=at$r0, r1=at$r1, lty=2)

#Add 201 CHH
at<-autotrack(current.track=1, total.tracks = 3)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X201_CHH_axi_Proportion, col="goldenrod2", pch=16, data.panel = 2, r0=at$r0, r1=at$r1)
kpLines(kp, data=windows.100k, y=mcols(windows.100k)$X201_CHH_epi_Proportion, col="cornflowerblue", pch=16, data.panel = 2, r0=at$r0, r1=at$r1, lty=2)
