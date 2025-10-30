library(DMRcaller)
library(GenomicRanges)
library(genomation)
library(tidyverse)

load("/users/PAS1755/whazzard/Epicormic_Axillary/8_160_pooled_samples_fixed.Rdata")

load("/users/PAS1755/whazzard/Epicormic_Axillary/8_201_pooled_samples_fixed.Rdata")

#make the granges objects for each chromosome
chr1 <- GRanges(seqnames = Rle("chr001"), ranges = IRanges(1,43996934))
chr2 <- GRanges(seqnames = Rle("chr002"), ranges = IRanges(1,26138581))
chr3 <- GRanges(seqnames = Rle("chr003"), ranges = IRanges(1,24073526))
chr4 <- GRanges(seqnames = Rle("chr004"), ranges = IRanges(1,24375383))
chr5 <- GRanges(seqnames = Rle("chr005"), ranges = IRanges(1,18233718))
chr6 <- GRanges(seqnames = Rle("chr006"), ranges = IRanges(1,29596931))
chr7 <- GRanges(seqnames = Rle("chr007"), ranges = IRanges(1,21340587))
chr8 <- GRanges(seqnames = Rle("chr008"), ranges = IRanges(1,20428636))
mit <- GRanges(seqnames = Rle("mit"), ranges = IRanges(1,1157723))

#create a granges list using those object
chr_list<- list(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, mit)

#create a seperate vector of strings containing the names of the chromosomes
chr_names<- list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "mit")


samples_8.160<- list(axillary.8160.billings,
                     axillary.8160.chico,
                     axillary.8160.chowchilla,
                     axillary.8160.salida,
                     epicormic.8160.billings,
                     epicormic.8160.chico,
                     epicormic.8160.chowchilla,
                     epicormic.8160.salida)

names(samples_8.160) <- c("axillary.8160.billings",
                          "axillary.8160.chico",
                          "axillary.8160.chowchilla",
                          "axillary.8160.salida",
                          "epicormic.8160.billings",
                          "epicormic.8160.chico",
                          "epicormic.8160.chowchilla",
                          "epicormic.8160.salida")

for (x in seq_along(chr_list)) {
  for (y in seq_along(samples_8.160)) {

    profile <- computeMethylationProfile(
      samples_8.160[[y]],
      chr_list[[x]],
      windowSize = 10000,
      context = "CG"
    )

    object_name <- paste0(names(samples_8.160)[y], ".", chr_names[x])

    assign(object_name, profile, envir = .GlobalEnv)
  }
}


sample_names <- c("axillary.8160.billings",
                  "axillary.8160.chico",
                  "axillary.8160.chowchilla",
                  "axillary.8160.salida",
                  "epicormic.8160.billings",
                  "epicormic.8160.chico",
                  "epicormic.8160.chowchilla",
                  "epicormic.8160.salida")

# Extract tissue and location from sample names for adding metadata later
tissue_types <- c("axillary", "axillary", "axillary", "axillary",
                  "epicormic", "epicormic", "epicormic", "epicormic")
locations <- c("billings", "chico", "chowchilla", "salida",
               "billings", "chico", "chowchilla", "salida")

# Loop through each sample in the samples_8.160 list
for (i in seq_along(samples_8.160)) {

  # Combine the chromosomes for the current sample
  sample_object <- do.call(c, lapply(chr_names, function(chr) {
    get(paste0(sample_names[i], ".", chr))
  }))

  # Convert to data frame
  CpG_df <- as.data.frame(sample_object)

  # Add treatment group metadata (tissue, location, variety)
  CpG_df <- CpG_df %>%
    mutate(tissue = tissue_types[i],
           location = locations[i],
           variety = "8_160")

  # Assign the data frame to a global variable with a dynamic name
  object_name <- paste0(sample_names[i], ".df")
  assign(object_name, CpG_df, envir = .GlobalEnv)

}


samples_8.201<- list(axillary.8201.billings,
                     axillary.8201.chico,
                     axillary.8201.chowchilla,
                     axillary.8201.salida,
                     epicormic.8201.billings,
                     epicormic.8201.chico,
                     epicormic.8201.chowchilla,
                     epicormic.8201.salida)

names(samples_8.201) <- c("axillary.8201.billings",
                          "axillary.8201.chico",
                          "axillary.8201.chowchilla",
                          "axillary.8201.salida",
                          "epicormic.8201.billings",
                          "epicormic.8201.chico",
                          "epicormic.8201.chowchilla",
                          "epicormic.8201.salida")

for (x in seq_along(chr_list)) {
  for (y in seq_along(samples_8.201)) {

    profile <- computeMethylationProfile(
      samples_8.201[[y]],
      chr_list[[x]],
      windowSize = 10000,
      context = "CG"
    )

    object_name <- paste0(names(samples_8.201)[y], ".", chr_names[x])

    assign(object_name, profile, envir = .GlobalEnv)
  }
}


sample_names <- c("axillary.8201.billings",
                  "axillary.8201.chico",
                  "axillary.8201.chowchilla",
                  "axillary.8201.salida",
                  "epicormic.8201.billings",
                  "epicormic.8201.chico",
                  "epicormic.8201.chowchilla",
                  "epicormic.8201.salida")

# Extract tissue and location from sample names for adding metadata later
tissue_types <- c("axillary", "axillary", "axillary", "axillary",
                  "epicormic", "epicormic", "epicormic", "epicormic")
locations <- c("billings", "chico", "chowchilla", "salida",
               "billings", "chico", "chowchilla", "salida")

# Loop through each sample in the samples_8.160 list
for (i in seq_along(samples_8.201)) {

  # Combine the chromosomes for the current sample
  sample_object <- do.call(c, lapply(chr_names, function(chr) {
    get(paste0(sample_names[i], ".", chr))
  }))

  # Convert to data frame
  CpG_df <- as.data.frame(sample_object)

  # Add treatment group metadata (tissue, location, variety)
  CpG_df <- CpG_df %>%
    mutate(tissue = tissue_types[i],
           location = locations[i],
           variety = "8_201")

  # Assign the data frame to a global variable with a dynamic name
  object_name <- paste0(sample_names[i], ".df")
  assign(object_name, CpG_df, envir = .GlobalEnv)

}

methylation_df_all_samples<-bind_rows(axillary.8160.billings.df,
                                      epicormic.8160.billings.df,
                                      axillary.8160.chico.df,
                                      epicormic.8160.chico.df,
                                      axillary.8160.chowchilla.df,
                                      epicormic.8160.chowchilla.df,
                                      axillary.8160.salida.df,
                                      epicormic.8160.salida.df,
                                      axillary.8201.billings.df,
                                      epicormic.8201.billings.df,
                                      axillary.8201.chico.df,
                                      epicormic.8201.chico.df,
                                      axillary.8201.chowchilla.df,
                                      epicormic.8201.chowchilla.df,
                                      axillary.8201.salida.df,
                                      epicormic.8201.salida.df)

save(axillary.8160.billings.df,
                                      epicormic.8160.billings.df,
                                      axillary.8160.chico.df,
                                      epicormic.8160.chico.df,
                                      axillary.8160.chowchilla.df,
                                      epicormic.8160.chowchilla.df,
                                      axillary.8160.salida.df,
                                      epicormic.8160.salida.df,
                                      axillary.8201.billings.df,
                                      epicormic.8201.billings.df,
                                      axillary.8201.chico.df,
                                      epicormic.8201.chico.df,
                                      axillary.8201.chowchilla.df,
                                      epicormic.8201.chowchilla.df,
                                      axillary.8201.salida.df,
                                      epicormic.8201.salida.df,
     methylation_df_all_samples,
     file="/users/PAS1755/whazzard/10000bp_cpg_sample_profiles_all.Rdata")
