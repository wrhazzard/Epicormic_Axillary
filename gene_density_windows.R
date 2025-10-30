library(GenomicRanges)
library(genomation)
library(rtracklayer)
library(tidyverse)
library(Biostrings)
library(DMRcaller)


#load in the GFF and FASTA
almond.gff<- gffToGRanges("/fs/scratch/PAS1755/Will/EpiAxi/Visualization/In_DMRcaller/context_report/genomic.gff")
almond.fasta<- readDNAStringSet("/fs/scratch/PAS1755/Will/EpiAxi/Genomes/Almond/Almond.fasta")

#adjust the RLE names in the GFF so they match with the Granges methylation object
seqlevels(almond.gff)<- sub("CM037988.1", "chr001", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037989.1", "chr002", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037990.1", "chr003", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037991.1", "chr004", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037992.1", "chr005", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037993.1", "chr006", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037994.1", "chr007", seqlevels(almond.gff))
seqlevels(almond.gff)<- sub("CM037995.1", "chr008", seqlevels(almond.gff))

#create genome sizes from FASTA
chr_sizes<-seqlengths(almond.fasta)
chr_sizes<-chr_sizes[1:8]
  

#subset the GFF to only contain the genes, ranges and related ID
genes<- almond.gff[which(almond.gff$type=="gene")]
genes<-genes[seqnames(genes)==c("chr001", "chr002","chr003","chr004","chr005","chr006","chr007","chr008")]
genes.df<-as.data.frame(genes)
genes.df<-genes.df %>% select(seqnames, start, end, strand, ID)

genes<-makeGRangesFromDataFrame(genes.df, keep.extra.columns = T)

#create 1k, 10k, and 100k windows
windows.100k<-tileGenome(seqlengths = chr_sizes, tilewidth = 100000, cut.last.tile.in.chrom = T)
windows.10k<-tileGenome(seqlengths = chr_sizes, tilewidth = 10000, cut.last.tile.in.chrom = T)
windows.1k<-tileGenome(seqlengths = chr_sizes, tilewidth = 1000, cut.last.tile.in.chrom = T)


#count the number of genes in a given window
#100k
windows.100k$totgenes<- countOverlaps(windows.100k, genes)
overlaps <- findOverlaps(windows.100k, genes)
gene_size <- numeric(length(windows.100k))
gene_size[queryHits(overlaps)] <- tapply(width(genes)[subjectHits(overlaps)], 
                                             queryHits(overlaps), sum, default = 0)
mcols(windows.100k)$gene_size <-gene_size


#10k
windows.10k$totgenes<- countOverlaps(windows.10k, genes)
overlaps <- findOverlaps(windows.10k, genes)
gene_size <- numeric(length(windows.10k))
gene_size[queryHits(overlaps)] <- tapply(width(genes)[subjectHits(overlaps)], 
                                         queryHits(overlaps), sum, default = 0)
mcols(windows.10k)$gene_size <-gene_size

#1k
windows.1k$totgenes<-countOverlaps(windows.1k, genes)
overlaps <- findOverlaps(windows.1k, genes)
gene_size <- numeric(length(windows.1k))
gene_size[queryHits(overlaps)] <- tapply(width(genes)[subjectHits(overlaps)], 
                                         queryHits(overlaps), sum, default = 0)
mcols(windows.1k)$gene_size <-gene_size


#make the mean methylation profiles for epicormic and axillary for 100k windows

load("/users/PAS1755/whazzard/Epicormic_Axillary/8_160_pooled_samples_fixed.Rdata")

epiaxilist160<- GRangesList("Epicormic" = epicormic.8160.billings,
                            epicormic.8160.salida,
                            epicormic.8160.chico, 
                            epicormic.8160.chowchilla,
                            "Axillary" = axillary.8160.billings,
                            axillary.8160.chico, 
                            axillary.8160.salida)
                            
load("/users/PAS1755/whazzard/Epicormic_Axillary/8_201_pooled_samples_fixed.Rdata")

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


#100k loops

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist160[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 100000,
                                       context = "CG")
  
  object_name <- paste0("CpG.profile.axi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist160[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 100000,
                                       context = "CG")
  
  object_name <- paste0("CpG.profile.epi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CpG.profile.axi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CpG.profile.axi.160.", "", profile_names.axi))]

CpG.160.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CpG.profile.epi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CpG.profile.epi.160.", "", profile_names.epi))]

CpG.160.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CpG.160.axi.profiles)) {
  x <- CpG.160.axi.profiles[[i]]
  y <- CpG.160.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CpG.160.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CpG.160.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CpG.160.profile.final.", "", final_profile_names))]

CpG.160.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary
CpG.160.axi.profiles.tot <- Reduce(c, CpG.160.axi.profiles)
CpG.160.epi.profiles.tot <- unlist(GRangesList(CpG.160.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.100k, CpG.160.axi.profiles.tot)

mcols(windows.100k)$Proportion<-NA

mcols(windows.100k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CpG.160.axi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.100k))[colnames(mcols(windows.100k)) == "Proportion"] <- "160_CpG_axi_Proportion"


#epicormic


CpG.160.epi.profiles.tot <- Reduce(c, CpG.160.epi.profiles)
CpG.160.epi.profiles.tot <- unlist(GRangesList(CpG.160.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.100k, CpG.160.epi.profiles.tot)

mcols(windows.100k)$Proportion<-NA

mcols(windows.100k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CpG.160.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.100k))[colnames(mcols(windows.100k)) == "Proportion"] <- "160_CpG_epi_Proportion"


#10k loop

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist160[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 10000,
                                       context = "CG")
  
  object_name <- paste0("CpG.profile.axi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist160[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 10000,
                                       context = "CG")
  
  object_name <- paste0("CpG.profile.epi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CpG.profile.axi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CpG.profile.axi.160.", "", profile_names.axi))]

CpG.160.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CpG.profile.epi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CpG.profile.epi.160.", "", profile_names.epi))]

CpG.160.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CpG.160.axi.profiles)) {
  x <- CpG.160.axi.profiles[[i]]
  y <- CpG.160.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CpG.160.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CpG.160.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CpG.160.profile.final.", "", final_profile_names))]

CpG.160.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary
CpG.160.axi.profiles.tot <- Reduce(c, CpG.160.axi.profiles)
CpG.160.epi.profiles.tot <- unlist(GRangesList(CpG.160.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.10k, CpG.160.axi.profiles.tot)

mcols(windows.10k)$Proportion<-NA

mcols(windows.10k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CpG.160.axi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.10k))[colnames(mcols(windows.10k)) == "Proportion"] <- "160_CpG_axi_Proportion"


#epicormic


CpG.160.epi.profiles.tot <- Reduce(c, CpG.160.epi.profiles)
CpG.160.epi.profiles.tot <- unlist(GRangesList(CpG.160.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.10k, CpG.160.epi.profiles.tot)

mcols(windows.10k)$Proportion<-NA

mcols(windows.10k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CpG.160.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.10k))[colnames(mcols(windows.10k)) == "Proportion"] <- "160_CpG_epi_Proportion"

#1k loop

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist160[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 1000,
                                       context = "CG")
  
  object_name <- paste0("CpG.profile.axi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist160[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 1000,
                                       context = "CG")
  
  object_name <- paste0("CpG.profile.epi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CpG.profile.axi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CpG.profile.axi.160.", "", profile_names.axi))]

CpG.160.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CpG.profile.epi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CpG.profile.epi.160.", "", profile_names.epi))]

CpG.160.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CpG.160.axi.profiles)) {
  x <- CpG.160.axi.profiles[[i]]
  y <- CpG.160.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CpG.160.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CpG.160.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CpG.160.profile.final.", "", final_profile_names))]

CpG.160.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary

CpG.160.axi.profiles.tot <- Reduce(c, CpG.160.axi.profiles)
CpG.160.epi.profiles.tot <- unlist(GRangesList(CpG.160.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.1k, CpG.160.axi.profiles.tot)

mcols(windows.1k)$Proportion<-NA

mcols(windows.1k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CpG.160.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.1k))[colnames(mcols(windows.1k)) == "Proportion"] <- "160_CpG_axi_Proportion"

#epicormic


CpG.160.epi.profiles.tot <- Reduce(c, CpG.160.epi.profiles)
CpG.160.epi.profiles.tot <- unlist(GRangesList(CpG.160.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.1k, CpG.160.epi.profiles.tot)

mcols(windows.1k)$Proportion<-NA

mcols(windows.1k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CpG.160.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.1k))[colnames(mcols(windows.1k)) == "Proportion"] <- "160_CpG_epi_Proportion"



#CHG

#100k loops

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist160[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 100000,
                                       context = "CHG")
  
  object_name <- paste0("CHG.profile.axi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist160[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 100000,
                                       context = "CHG")
  
  object_name <- paste0("CHG.profile.epi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CHG.profile.axi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CHG.profile.axi.160.", "", profile_names.axi))]

CHG.160.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CHG.profile.epi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CHG.profile.epi.160.", "", profile_names.epi))]

CHG.160.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CHG.160.axi.profiles)) {
  x <- CHG.160.axi.profiles[[i]]
  y <- CHG.160.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CHG.160.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CHG.160.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CHG.160.profile.final.", "", final_profile_names))]

CHG.160.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary
CHG.160.axi.profiles.tot <- Reduce(c, CHG.160.axi.profiles)
CHG.160.axi.profiles.tot <- unlist(GRangesList(CHG.160.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.100k, CHG.160.axi.profiles.tot)

mcols(windows.100k)$Proportion<-NA

mcols(windows.100k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHG.160.axi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.100k))[colnames(mcols(windows.100k)) == "Proportion"] <- "160_CHG_axi_Proportion"


#epicormic


CHG.160.epi.profiles.tot <- Reduce(c, CHG.160.epi.profiles)
CHG.160.epi.profiles.tot <- unlist(GRangesList(CHG.160.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.100k, CHG.160.epi.profiles.tot)

mcols(windows.100k)$Proportion<-NA

mcols(windows.100k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHG.160.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.100k))[colnames(mcols(windows.100k)) == "Proportion"] <- "160_CHG_epi_Proportion"


#10k loop

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist160[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 10000,
                                       context = "CHG")
  
  object_name <- paste0("CHG.profile.axi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist160[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 10000,
                                       context = "CHG")
  
  object_name <- paste0("CHG.profile.epi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CHG.profile.axi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CHG.profile.axi.160.", "", profile_names.axi))]

CHG.160.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CHG.profile.epi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CHG.profile.epi.160.", "", profile_names.epi))]

CHG.160.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CHG.160.axi.profiles)) {
  x <- CHG.160.axi.profiles[[i]]
  y <- CHG.160.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CHG.160.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CHG.160.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CHG.160.profile.final.", "", final_profile_names))]

CHG.160.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary
CHG.160.axi.profiles.tot <- Reduce(c, CHG.160.axi.profiles)
CHG.160.epi.profiles.tot <- unlist(GRangesList(CHG.160.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.10k, CHG.160.axi.profiles.tot)

mcols(windows.10k)$Proportion<-NA

mcols(windows.10k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHG.160.axi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.10k))[colnames(mcols(windows.10k)) == "Proportion"] <- "160_CHG_axi_Proportion"


#epicormic


CHG.160.epi.profiles.tot <- Reduce(c, CHG.160.epi.profiles)
CHG.160.epi.profiles.tot <- unlist(GRangesList(CHG.160.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.10k, CHG.160.epi.profiles.tot)

mcols(windows.10k)$Proportion<-NA

mcols(windows.10k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHG.160.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.10k))[colnames(mcols(windows.10k)) == "Proportion"] <- "160_CHG_epi_Proportion"

#1k loop

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist160[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 1000,
                                       context = "CHG")
  
  object_name <- paste0("CHG.profile.axi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist160[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 1000,
                                       context = "CHG")
  
  object_name <- paste0("CHG.profile.epi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CHG.profile.axi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CHG.profile.axi.160.", "", profile_names.axi))]

CHG.160.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CHG.profile.epi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CHG.profile.epi.160.", "", profile_names.epi))]

CHG.160.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CpG.160.axi.profiles)) {
  x <- CHG.160.axi.profiles[[i]]
  y <- CHG.160.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CHG.160.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CHG.160.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CHG.160.profile.final.", "", final_profile_names))]

CHG.160.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary

CHG.160.axi.profiles.tot <- Reduce(c, CHG.160.axi.profiles)
CHG.160.epi.profiles.tot <- unlist(GRangesList(CHG.160.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.1k, CHG.160.axi.profiles.tot)

mcols(windows.1k)$Proportion<-NA

mcols(windows.1k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHG.160.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.1k))[colnames(mcols(windows.1k)) == "Proportion"] <- "160_CHG_axi_Proportion"

#epicormic


CHG.160.epi.profiles.tot <- Reduce(c, CHG.160.epi.profiles)
CHG.160.epi.profiles.tot <- unlist(GRangesList(CHG.160.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.1k, CHG.160.epi.profiles.tot)

mcols(windows.1k)$Proportion<-NA

mcols(windows.1k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHG.160.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.1k))[colnames(mcols(windows.1k)) == "Proportion"] <- "160_CHG_epi_Proportion"


#CHH

#100k loops

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist160[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 100000,
                                       context = "CHH")
  
  object_name <- paste0("CHH.profile.axi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist160[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 100000,
                                       context = "CHH")
  
  object_name <- paste0("CHH.profile.epi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CHH.profile.axi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CHH.profile.axi.160.", "", profile_names.axi))]

CHH.160.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CHH.profile.epi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CHH.profile.epi.160.", "", profile_names.epi))]

CHH.160.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CHH.160.axi.profiles)) {
  x <- CHH.160.axi.profiles[[i]]
  y <- CHH.160.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CHH.160.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CHH.160.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CHH.160.profile.final.", "", final_profile_names))]

CHH.160.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary
CHH.160.axi.profiles.tot <- Reduce(c, CHH.160.axi.profiles)
CHH.160.axi.profiles.tot <- unlist(GRangesList(CHH.160.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.100k, CHH.160.axi.profiles.tot)

mcols(windows.100k)$Proportion<-NA

mcols(windows.100k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHH.160.axi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.100k))[colnames(mcols(windows.100k)) == "Proportion"] <- "160_CHH_axi_Proportion"


#epicormic


CHH.160.epi.profiles.tot <- Reduce(c, CHH.160.epi.profiles)
CHH.160.epi.profiles.tot <- unlist(GRangesList(CHH.160.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.100k, CHH.160.epi.profiles.tot)

mcols(windows.100k)$Proportion<-NA

mcols(windows.100k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHH.160.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.100k))[colnames(mcols(windows.100k)) == "Proportion"] <- "160_CHH_epi_Proportion"


#10k loop

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist160[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 10000,
                                       context = "CHH")
  
  object_name <- paste0("CHH.profile.axi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist160[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 10000,
                                       context = "CHH")
  
  object_name <- paste0("CHH.profile.epi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CHH.profile.axi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CHH.profile.axi.160.", "", profile_names.axi))]

CHH.160.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CHH.profile.epi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CHH.profile.epi.160.", "", profile_names.epi))]

CHH.160.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CHG.160.axi.profiles)) {
  x <- CHH.160.axi.profiles[[i]]
  y <- CHH.160.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CHH.160.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CHH.160.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CHH.160.profile.final.", "", final_profile_names))]

CHH.160.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary
CHH.160.axi.profiles.tot <- Reduce(c, CHH.160.axi.profiles)
CHH.160.epi.profiles.tot <- unlist(GRangesList(CHH.160.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.10k, CHH.160.axi.profiles.tot)

mcols(windows.10k)$Proportion<-NA

mcols(windows.10k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHH.160.axi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.10k))[colnames(mcols(windows.10k)) == "Proportion"] <- "160_CHH_axi_Proportion"


#epicormic


CHH.160.epi.profiles.tot <- Reduce(c, CHH.160.epi.profiles)
CHH.160.epi.profiles.tot <- unlist(GRangesList(CHH.160.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.10k, CHH.160.epi.profiles.tot)

mcols(windows.10k)$Proportion<-NA

mcols(windows.10k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHH.160.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.10k))[colnames(mcols(windows.10k)) == "Proportion"] <- "160_CHH_epi_Proportion"

#1k loop

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist160[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 1000,
                                       context = "CHH")
  
  object_name <- paste0("CHH.profile.axi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist160[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 1000,
                                       context = "CHH")
  
  object_name <- paste0("CHH.profile.epi.160.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CHH.profile.axi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CHH.profile.axi.160.", "", profile_names.axi))]

CHH.160.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CHH.profile.epi.160.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CHH.profile.epi.160.", "", profile_names.epi))]

CHH.160.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CpG.160.axi.profiles)) {
  x <- CHH.160.axi.profiles[[i]]
  y <- CHH.160.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CHH.160.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CHH.160.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CHH.160.profile.final.", "", final_profile_names))]

CHH.160.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary

CHH.160.axi.profiles.tot <- Reduce(c, CHH.160.axi.profiles)
CHH.160.epi.profiles.tot <- unlist(GRangesList(CHH.160.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.1k, CHH.160.axi.profiles.tot)

mcols(windows.1k)$Proportion<-NA

mcols(windows.1k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHH.160.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.1k))[colnames(mcols(windows.1k)) == "Proportion"] <- "160_CHH_axi_Proportion"

#epicormic


CHH.160.epi.profiles.tot <- Reduce(c, CHH.160.epi.profiles)
CHH.160.epi.profiles.tot <- unlist(GRangesList(CHH.160.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.1k, CHH.160.epi.profiles.tot)

mcols(windows.1k)$Proportion<-NA

mcols(windows.1k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHH.160.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.1k))[colnames(mcols(windows.1k)) == "Proportion"] <- "160_CHH_epi_Proportion"


#8_201

#100k loops

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist201[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 100000,
                                       context = "CG")
  
  object_name <- paste0("CpG.profile.axi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist201[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 100000,
                                       context = "CG")
  
  object_name <- paste0("CpG.profile.epi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CpG.profile.axi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CpG.profile.axi.201.", "", profile_names.axi))]

CpG.201.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CpG.profile.epi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CpG.profile.epi.201.", "", profile_names.epi))]

CpG.201.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CpG.201.axi.profiles)) {
  x <- CpG.201.axi.profiles[[i]]
  y <- CpG.201.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CpG.201.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CpG.201.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CpG.201.profile.final.", "", final_profile_names))]

CpG.201.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary
CpG.201.axi.profiles.tot <- Reduce(c, CpG.201.axi.profiles)
CpG.201.epi.profiles.tot <- unlist(GRangesList(CpG.201.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.100k, CpG.201.axi.profiles.tot)

mcols(windows.100k)$Proportion<-NA

mcols(windows.100k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CpG.201.axi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.100k))[colnames(mcols(windows.100k)) == "Proportion"] <- "201_CpG_axi_Proportion"


#epicormic


CpG.201.epi.profiles.tot <- Reduce(c, CpG.201.epi.profiles)
CpG.201.epi.profiles.tot <- unlist(GRangesList(CpG.201.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.100k, CpG.201.epi.profiles.tot)

mcols(windows.100k)$Proportion<-NA

mcols(windows.100k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CpG.201.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.100k))[colnames(mcols(windows.100k)) == "Proportion"] <- "201_CpG_epi_Proportion"


#10k loop

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist201[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 10000,
                                       context = "CG")
  
  object_name <- paste0("CpG.profile.axi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist201[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 10000,
                                       context = "CG")
  
  object_name <- paste0("CpG.profile.epi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CpG.profile.axi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CpG.profile.axi.201.", "", profile_names.axi))]

CpG.201.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CpG.profile.epi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CpG.profile.epi.201.", "", profile_names.epi))]

CpG.201.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CpG.201.axi.profiles)) {
  x <- CpG.201.axi.profiles[[i]]
  y <- CpG.201.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CpG.201.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CpG.201.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CpG.201.profile.final.", "", final_profile_names))]

CpG.201.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary
CpG.201.axi.profiles.tot <- Reduce(c, CpG.201.axi.profiles)
CpG.201.epi.profiles.tot <- unlist(GRangesList(CpG.201.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.10k, CpG.201.axi.profiles.tot)

mcols(windows.10k)$Proportion<-NA

mcols(windows.10k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CpG.201.axi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.10k))[colnames(mcols(windows.10k)) == "Proportion"] <- "201_CpG_axi_Proportion"


#epicormic


CpG.201.epi.profiles.tot <- Reduce(c, CpG.201.epi.profiles)
CpG.201.epi.profiles.tot <- unlist(GRangesList(CpG.201.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.10k, CpG.201.epi.profiles.tot)

mcols(windows.10k)$Proportion<-NA

mcols(windows.10k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CpG.201.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.10k))[colnames(mcols(windows.10k)) == "Proportion"] <- "201_CpG_epi_Proportion"

#1k loop

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist201[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 1000,
                                       context = "CG")
  
  object_name <- paste0("CpG.profile.axi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist201[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 1000,
                                       context = "CG")
  
  object_name <- paste0("CpG.profile.epi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CpG.profile.axi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CpG.profile.axi.201.", "", profile_names.axi))]

CpG.201.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CpG.profile.epi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CpG.profile.epi.201.", "", profile_names.epi))]

CpG.201.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CpG.201.axi.profiles)) {
  x <- CpG.201.axi.profiles[[i]]
  y <- CpG.201.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CpG.201.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CpG.201.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CpG.201.profile.final.", "", final_profile_names))]

CpG.201.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary

CpG.201.axi.profiles.tot <- Reduce(c, CpG.201.axi.profiles)
CpG.201.epi.profiles.tot <- unlist(GRangesList(CpG.201.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.1k, CpG.201.axi.profiles.tot)

mcols(windows.1k)$Proportion<-NA

mcols(windows.1k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CpG.201.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.1k))[colnames(mcols(windows.1k)) == "Proportion"] <- "201_CpG_axi_Proportion"

#epicormic


CpG.201.epi.profiles.tot <- Reduce(c, CpG.201.epi.profiles)
CpG.201.epi.profiles.tot <- unlist(GRangesList(CpG.201.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.1k, CpG.201.epi.profiles.tot)

mcols(windows.1k)$Proportion<-NA

mcols(windows.1k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CpG.201.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.1k))[colnames(mcols(windows.1k)) == "Proportion"] <- "201_CpG_epi_Proportion"



#CHG

#100k loops

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist201[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 100000,
                                       context = "CHG")
  
  object_name <- paste0("CHG.profile.axi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist201[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 100000,
                                       context = "CHG")
  
  object_name <- paste0("CHG.profile.epi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CHG.profile.axi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CHG.profile.axi.201.", "", profile_names.axi))]

CHG.201.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CHG.profile.epi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CHG.profile.epi.201.", "", profile_names.epi))]

CHG.201.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CHG.201.axi.profiles)) {
  x <- CHG.201.axi.profiles[[i]]
  y <- CHG.201.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CHG.201.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CHG.201.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CHG.201.profile.final.", "", final_profile_names))]

CHG.201.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary
CHG.201.axi.profiles.tot <- Reduce(c, CHG.201.axi.profiles)
CHG.201.axi.profiles.tot <- unlist(GRangesList(CHG.201.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.100k, CHG.201.axi.profiles.tot)

mcols(windows.100k)$Proportion<-NA

mcols(windows.100k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHG.201.axi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.100k))[colnames(mcols(windows.100k)) == "Proportion"] <- "201_CHG_axi_Proportion"


#epicormic


CHG.201.epi.profiles.tot <- Reduce(c, CHG.201.epi.profiles)
CHG.201.epi.profiles.tot <- unlist(GRangesList(CHG.201.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.100k, CHG.201.epi.profiles.tot)

mcols(windows.100k)$Proportion<-NA

mcols(windows.100k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHG.201.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.100k))[colnames(mcols(windows.100k)) == "Proportion"] <- "201_CHG_epi_Proportion"


#10k loop

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist201[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 10000,
                                       context = "CHG")
  
  object_name <- paste0("CHG.profile.axi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist201[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 10000,
                                       context = "CHG")
  
  object_name <- paste0("CHG.profile.epi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CHG.profile.axi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CHG.profile.axi.201.", "", profile_names.axi))]

CHG.201.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CHG.profile.epi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CHG.profile.epi.201.", "", profile_names.epi))]

CHG.201.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CHG.201.axi.profiles)) {
  x <- CHG.201.axi.profiles[[i]]
  y <- CHG.201.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CHG.201.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CHG.201.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CHG.201.profile.final.", "", final_profile_names))]

CHG.201.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary
CHG.201.axi.profiles.tot <- Reduce(c, CHG.201.axi.profiles)
CHG.201.epi.profiles.tot <- unlist(GRangesList(CHG.201.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.10k, CHG.201.axi.profiles.tot)

mcols(windows.10k)$Proportion<-NA

mcols(windows.10k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHG.201.axi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.10k))[colnames(mcols(windows.10k)) == "Proportion"] <- "201_CHG_axi_Proportion"


#epicormic


CHG.201.epi.profiles.tot <- Reduce(c, CHG.201.epi.profiles)
CHG.201.epi.profiles.tot <- unlist(GRangesList(CHG.201.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.10k, CHG.201.epi.profiles.tot)

mcols(windows.10k)$Proportion<-NA

mcols(windows.10k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHG.201.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.10k))[colnames(mcols(windows.10k)) == "Proportion"] <- "201_CHG_epi_Proportion"

#1k loop

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist201[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 1000,
                                       context = "CHG")
  
  object_name <- paste0("CHG.profile.axi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist201[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 1000,
                                       context = "CHG")
  
  object_name <- paste0("CHG.profile.epi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CHG.profile.axi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CHG.profile.axi.201.", "", profile_names.axi))]

CHG.201.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CHG.profile.epi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CHG.profile.epi.201.", "", profile_names.epi))]

CHG.201.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CpG.201.axi.profiles)) {
  x <- CHG.201.axi.profiles[[i]]
  y <- CHG.201.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CHG.201.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CHG.201.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CHG.201.profile.final.", "", final_profile_names))]

CHG.201.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary

CHG.201.axi.profiles.tot <- Reduce(c, CHG.201.axi.profiles)
CHG.201.epi.profiles.tot <- unlist(GRangesList(CHG.201.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.1k, CHG.201.axi.profiles.tot)

mcols(windows.1k)$Proportion<-NA

mcols(windows.1k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHG.201.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.1k))[colnames(mcols(windows.1k)) == "Proportion"] <- "201_CHG_axi_Proportion"

#epicormic


CHG.201.epi.profiles.tot <- Reduce(c, CHG.201.epi.profiles)
CHG.201.epi.profiles.tot <- unlist(GRangesList(CHG.201.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.1k, CHG.201.epi.profiles.tot)

mcols(windows.1k)$Proportion<-NA

mcols(windows.1k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHG.201.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.1k))[colnames(mcols(windows.1k)) == "Proportion"] <- "201_CHG_epi_Proportion"


#CHH

#100k loops

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist201[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 100000,
                                       context = "CHH")
  
  object_name <- paste0("CHH.profile.axi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist201[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 100000,
                                       context = "CHH")
  
  object_name <- paste0("CHH.profile.epi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CHH.profile.axi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CHH.profile.axi.201.", "", profile_names.axi))]

CHH.201.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CHH.profile.epi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CHH.profile.epi.201.", "", profile_names.epi))]

CHH.201.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CHH.201.axi.profiles)) {
  x <- CHH.201.axi.profiles[[i]]
  y <- CHH.201.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CHH.201.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CHH.201.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CHH.201.profile.final.", "", final_profile_names))]

CHH.201.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary
CHH.201.axi.profiles.tot <- Reduce(c, CHH.201.axi.profiles)
CHH.201.axi.profiles.tot <- unlist(GRangesList(CHH.201.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.100k, CHH.201.axi.profiles.tot)

mcols(windows.100k)$Proportion<-NA

mcols(windows.100k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHH.201.axi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.100k))[colnames(mcols(windows.100k)) == "Proportion"] <- "201_CHH_axi_Proportion"


#epicormic


CHH.201.epi.profiles.tot <- Reduce(c, CHH.201.epi.profiles)
CHH.201.epi.profiles.tot <- unlist(GRangesList(CHH.201.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.100k, CHH.201.epi.profiles.tot)

mcols(windows.100k)$Proportion<-NA

mcols(windows.100k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHH.201.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.100k))[colnames(mcols(windows.100k)) == "Proportion"] <- "201_CHH_epi_Proportion"


#10k loop

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist201[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 10000,
                                       context = "CHH")
  
  object_name <- paste0("CHH.profile.axi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist201[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 10000,
                                       context = "CHH")
  
  object_name <- paste0("CHH.profile.epi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CHH.profile.axi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CHH.profile.axi.201.", "", profile_names.axi))]

CHH.201.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CHH.profile.epi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CHH.profile.epi.201.", "", profile_names.epi))]

CHH.201.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CHG.201.axi.profiles)) {
  x <- CHH.201.axi.profiles[[i]]
  y <- CHH.201.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CHH.201.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CHH.201.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CHH.201.profile.final.", "", final_profile_names))]

CHH.201.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary
CHH.201.axi.profiles.tot <- Reduce(c, CHH.201.axi.profiles)
CHH.201.epi.profiles.tot <- unlist(GRangesList(CHH.201.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.10k, CHH.201.axi.profiles.tot)

mcols(windows.10k)$Proportion<-NA

mcols(windows.10k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHH.201.axi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.10k))[colnames(mcols(windows.10k)) == "Proportion"] <- "201_CHH_axi_Proportion"


#epicormic


CHH.201.epi.profiles.tot <- Reduce(c, CHH.201.epi.profiles)
CHH.201.epi.profiles.tot <- unlist(GRangesList(CHH.201.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.10k, CHH.201.epi.profiles.tot)

mcols(windows.10k)$Proportion<-NA

mcols(windows.10k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHH.201.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.10k))[colnames(mcols(windows.10k)) == "Proportion"] <- "201_CHH_epi_Proportion"

#1k loop

#this loop builds axillary methylation profiles for each chromosome and then assigns them to an object
for (x in seq_along(chr_list)){
  profile <- computeMethylationProfile(epiaxilist201[["Axillary"]],
                                       chr_list[[x]],
                                       windowSize = 1000,
                                       context = "CHH")
  
  object_name <- paste0("CHH.profile.axi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}

#this loop builds epicormic methylation profiles for each chromosome and then assigns them to an objec
for (x in seq_along(chr_list)) {
  profile <- computeMethylationProfile(epiaxilist201[["Epicormic"]],
                                       chr_list[[x]],
                                       windowSize = 1000,
                                       context = "CHH")
  
  object_name <- paste0("CHH.profile.epi.201.", chr_names[x])
  
  assign(object_name, profile, envir = .GlobalEnv)
}


#the profiles are then merged into a single list containing all the profiles

#axillary
profile_names.axi <- grep("CHH.profile.axi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.axi <- profile_names.axi[match(chr_names, gsub("CHH.profile.axi.201.", "", profile_names.axi))]

CHH.201.axi.profiles <- mget(ordered_profile_names.axi, envir = .GlobalEnv)


#epicormic
profile_names.epi <- grep("CHH.profile.epi.201.", names(.GlobalEnv), value = TRUE)

ordered_profile_names.epi <- profile_names.epi[match(chr_names, gsub("CHH.profile.epi.201.", "", profile_names.epi))]

CHH.201.epi.profiles <- mget(ordered_profile_names.epi, envir = .GlobalEnv)


#this loop takes the profile from the epicormic and axillary lists and puts them into a single profile for each chromsome  
for (i in seq_along(CpG.201.axi.profiles)) {
  x <- CHH.201.axi.profiles[[i]]
  y <- CHH.201.epi.profiles[[i]]
  
  profiles <- GRangesList("Axillary" = x, "Epicormic" = y)
  
  object_name <- paste0("CHH.201.profile.final.", chr_names[i])
  
  assign(object_name, profiles, envir = .GlobalEnv)
}


#this grabs all those profiles and puts them in a final list
final_profile_names <- grep("CHH.201.profile.final.", names(.GlobalEnv), value = TRUE)

ordered_final_profile_names <- final_profile_names[match(chr_names, gsub("CHH.201.profile.final.", "", final_profile_names))]

CHH.201.final.profiles <- mget(ordered_final_profile_names, envir = .GlobalEnv)


#add the methylation proportions from granges list to the window gene counts for the 100k windows

#axillary

CHH.201.axi.profiles.tot <- Reduce(c, CHH.201.axi.profiles)
CHH.201.epi.profiles.tot <- unlist(GRangesList(CHH.201.axi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.1k, CHH.201.axi.profiles.tot)

mcols(windows.1k)$Proportion<-NA

mcols(windows.1k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHH.201.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.1k))[colnames(mcols(windows.1k)) == "Proportion"] <- "201_CHH_axi_Proportion"

#epicormic


CHH.201.epi.profiles.tot <- Reduce(c, CHH.201.epi.profiles)
CHH.160.epi.profiles.tot <- unlist(GRangesList(CHH.201.epi.profiles), use.names = FALSE)

overlaps<- findOverlaps(windows.1k, CHH.201.epi.profiles.tot)

mcols(windows.1k)$Proportion<-NA

mcols(windows.1k)[queryHits(overlaps), "Proportion"] <- 
  mcols(CHH.201.epi.profiles.tot)[subjectHits(overlaps), "Proportion"]


colnames(mcols(windows.1k))[colnames(mcols(windows.1k)) == "Proportion"] <- "201_CHH_epi_Proportion"


save(windows.100k, windows.10k, windows.1k, file = "gene_density_windows.Rdata")


