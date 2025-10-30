library(DMRcaller)
library(GenomicRanges)
library(genomation)

load("/users/PAS1755/whazzard/Epicormic_Axillary/8_201_pooled_samples_fixed.Rdata")

#create epiaxi list for 8-201
epiaxilist201<-GRangesList("Epicormic" = epicormic.8201.billings,
                           epicormic.8201.salida,
                           epicormic.8201.chico,
                           epicormic.8201.chowchilla,
                          "Axillary" = axillary.8201.billings,
                          axillary.8201.chico,
                          axillary.8201.chowchilla,
                          axillary.8201.salida)

#make the granges objects for each chromosome
chr1 <- GRanges(seqnames = Rle("chr001"), ranges = IRanges(1,53037566))
chr2 <- GRanges(seqnames = Rle("chr002"), ranges = IRanges(1,33982222))
chr3 <- GRanges(seqnames = Rle("chr003"), ranges = IRanges(1,30537282))
chr4 <- GRanges(seqnames = Rle("chr004"), ranges = IRanges(1,30643767))
chr5 <- GRanges(seqnames = Rle("chr005"), ranges = IRanges(1,21032966))
chr6 <- GRanges(seqnames = Rle("chr006"), ranges = IRanges(1,32793846))
chr7 <- GRanges(seqnames = Rle("chr007"), ranges = IRanges(1,26338084))
chr8 <- GRanges(seqnames = Rle("chr008"), ranges = IRanges(1,26044989))

#create a granges list using those object
chr_list<- list(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8)

#create a seperate vector of strings containing the names of the chromosomes
chr_names<- list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")

#load gff of NP geneome fro NCBI
almond.gff<- gffToGRanges("/fs/scratch/PAS1755/Will/EpiAxi/Visualization/In_DMRcaller/context_report/genomic.gff")

#adjust the RLE names in the GFF so they match with the Granges methylation object
seqlevels(almond.gff)<- sub("CM037988.1", "chr001", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037989.1", "chr002", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037990.1", "chr003", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037991.1", "chr004", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037992.1", "chr005", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037993.1", "chr006", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037994.1", "chr007", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037995.1", "chr008", seqlevels(almond.gff))


#subset by chromosome
anno_chr_1<- almond.gff[seqnames(almond.gff)=="chr001"]
anno_chr_2<- almond.gff[seqnames(almond.gff)=="chr002"]
anno_chr_3<- almond.gff[seqnames(almond.gff)=="chr003"]
anno_chr_4<- almond.gff[seqnames(almond.gff)=="chr004"]
anno_chr_5<- almond.gff[seqnames(almond.gff)=="chr005"]
anno_chr_6<- almond.gff[seqnames(almond.gff)=="chr006"]
anno_chr_7<- almond.gff[seqnames(almond.gff)=="chr007"]
anno_chr_8<- almond.gff[seqnames(almond.gff)=="chr008"]

DMRsGenesCG1<-filterDMRs(epiaxilist201[["Axillary"]],
                   epiaxilist201[["Epicormic"]],
potentialDMRs = anno_chr_1,
context = "CHH",
test = "score",
pValueThreshold = 0.01,
minCytosinesCount = 4,
minProportionDifference = 0.4,
minReadsPerCytosine = 3,
cores = 28)

DMRsGenesCG2<-filterDMRs(epiaxilist201[["Axillary"]],
                   epiaxilist201[["Epicormic"]],
potentialDMRs = anno_chr_2[overlapsAny(anno_chr_2, chr2)],
context = "CHH",
test = "score",
pValueThreshold = 0.01,
minCytosinesCount = 4,
minProportionDifference = 0.4,
minReadsPerCytosine = 3,
cores = 28)

DMRsGenesCG3<-filterDMRs(epiaxilist201[["Axillary"]],
                   epiaxilist201[["Epicormic"]],
potentialDMRs = anno_chr_3[overlapsAny(anno_chr_3, chr3)],
context = "CHH",
test = "score",
pValueThreshold = 0.01,
minCytosinesCount = 4,
minProportionDifference = 0.4,
minReadsPerCytosine = 3,
cores = 28)

DMRsGenesCG4<-filterDMRs(epiaxilist201[["Axillary"]],
                   epiaxilist201[["Epicormic"]],
potentialDMRs = anno_chr_4[overlapsAny(anno_chr_4, chr4)],
context = "CHH",
test = "score",
pValueThreshold = 0.01,
minCytosinesCount = 4,
minProportionDifference = 0.4,
minReadsPerCytosine = 3,
cores = 28)

DMRsGenesCG5<-filterDMRs(epiaxilist201[["Axillary"]],
                   epiaxilist201[["Epicormic"]],
potentialDMRs = anno_chr_5[overlapsAny(anno_chr_5, chr5)],
context = "CHH",
test = "score",
pValueThreshold = 0.01,
minCytosinesCount = 4,
minProportionDifference = 0.4,
minReadsPerCytosine = 3,
cores = 28)

DMRsGenesCG6<-filterDMRs(epiaxilist201[["Axillary"]],
                   epiaxilist201[["Epicormic"]],
potentialDMRs = anno_chr_6[overlapsAny(anno_chr_6, chr6)],
context = "CHH",
test = "score",
pValueThreshold = 0.01,
minCytosinesCount = 4,
minProportionDifference = 0.4,
minReadsPerCytosine = 3,
cores = 28)

DMRsGenesCG7<-filterDMRs(epiaxilist201[["Axillary"]],
                   epiaxilist201[["Epicormic"]],
potentialDMRs = anno_chr_7[overlapsAny(anno_chr_7, chr7)],
context = "CHH",
test = "score",
pValueThreshold = 0.01,
minCytosinesCount = 4,
minProportionDifference = 0.4,
minReadsPerCytosine = 3,
cores = 28)

DMRsGenesCG<-filterDMRs(epiaxilist201[["Axillary"]],
                   epiaxilist201[["Epicormic"]],
potentialDMRs = anno_chr_8[overlapsAny(anno_chr_8, chr8)],
context = "CHH",
test = "score",
pValueThreshold = 0.01,
minCytosinesCount = 4,
minProportionDifference = 0.4,
minReadsPerCytosine = 3,
cores = 28)


DMRs8201CHH_all<- c(DMRsGenesCG1, DMRsGenesCG2, DMRsGenesCG3, DMRsGenesCG4, DMRsGenesCG5, DMRsGenesCG6, DMRsGenesCG7, DMRsGenesCG)


DMRs8201CHH_all_df<- as.data.frame(DMRs8201CHH_all)

write.csv(unlist(DMRs8201CHH_all_df), "/users/PAS1755/whazzard//DMRCHH201_annotation.csv")

save(DMRsGenesCG, DMRsGenesCG1, DMRsGenesCG2, DMRsGenesCG3, DMRsGenesCG4, DMRsGenesCG5, DMRsGenesCG6, DMRsGenesCG7, file = "/users/PAS1755/whazzard/DMRCHH201_annotation.RData")
