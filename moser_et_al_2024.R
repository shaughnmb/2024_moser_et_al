#!/bin/env Rscript

#### SCRIPT INFORMATION ####
  #### License Notice ####

##
# Copyright (c) 2024 Cedars-Sinai Medical Center
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
  
  #### Title and Authors ####
  
##
# Title: Human iPSC-derived mononuclear phagocytes restore cognition and neural 
#        health across multiple mouse models of aging and Alzheimer’s disease
#
# Authors:  V. Alexandra Moser 1, Rachel M. Lipman 1*, Luz Jovita Dimas-Harms 1*, 
#           Jake Inzalaco 1*, Shaughn Bell 1, George Lawless 1, Simion Kreimer 2, 
#           Tao Zhang 1, Sarah J. Parker 2, Helen S. Goodridge 1, 
#           Clive N. Svendsen 1,3
# 
# Affiliations: 1 Board of Governors Regenerative Medicine Institute, 
#                 Cedars-Sinai Medical Center, Los Angeles, CA, USA. 
#               2 Smidt Heart Institute, Department of Cardiology, 
#                 Cedars-Sinai Medical Center, Los Angeles, CA, USA.
#               3 Lead contact
#  
#  * These authors contributed equally. 
#  
#  Corresponding authors' email:  V. Alexandra Moser: alexandra.moser@cshs.org
#                                 Clive N. Svendsen: clive.svendsen@cshs.org 
# 
##

  #### Script Information ####
  
##
# R version 4.3.1
# R Script Title:  moser_et_al_2024.R
# R Script Author:  Shaughn Bell
# R Script Corresponding Email:  shaughn.bell@cshs.org
#
# Notes: 
#   A) Script makes use of variables set up under "project information" as 
#      well as additional "prefixes" throughout the script for ease of saving
#      files with a similar path and naming structure.  When reloading data 
#      (see note "B"), you must either use the full path or reload the prefixes.
#   B) Script saves intermediate steps at each major manipulation of the seurat
#      object via the "saveRDS" function.  If needed, these RDS objects can then 
#      be reloaded to easily restart at one of these save points without needing 
#      to start from scratch.  However, these are not required for analysis, and
#      they can be skipped to save time and disk space.
##

  #### Sample Information #### 
    #### RNAseq Data ####
##
# Single Nucleus RNAseq data
#
# GEO Accession Number:  GSE220548
#
# aligned to mus musculus mm10 via CellRanger v7.1.0
# used the 10x Genomics supplied reference file "refdata-gex-mm10-2020-A.tar.gz"
#
# fastq files used:
#
# 102YVEH_S2_L001_I1_001.fastq.gz
# 102YVEH_S2_L001_R1_001.fastq.gz
# 102YVEH_S2_L001_R2_001.fastq.gz
# 102YVEH_S2_L002_I1_001.fastq.gz
# 102YVEH_S2_L002_R1_001.fastq.gz
# 102YVEH_S2_L002_R2_001.fastq.gz
# 104YVEH_S1_L001_I1_001.fastq.gz
# 104YVEH_S1_L001_R1_001.fastq.gz
# 104YVEH_S1_L001_R2_001.fastq.gz
# 104YVEH_S1_L002_I1_001.fastq.gz
# 104YVEH_S1_L002_R1_001.fastq.gz
# 104YVEH_S1_L002_R2_001.fastq.gz
# 106YVEH_S5_L001_I1_001.fastq.gz
# 106YVEH_S5_L001_R1_001.fastq.gz
# 106YVEH_S5_L001_R2_001.fastq.gz
# 106YVEH_S5_L002_I1_001.fastq.gz
# 106YVEH_S5_L002_R1_001.fastq.gz
# 106YVEH_S5_L002_R2_001.fastq.gz
# 111YVEH_S3_L001_I1_001.fastq.gz
# 111YVEH_S3_L001_R1_001.fastq.gz
# 111YVEH_S3_L001_R2_001.fastq.gz
# 111YVEH_S3_L002_I1_001.fastq.gz
# 111YVEH_S3_L002_R1_001.fastq.gz
# 111YVEH_S3_L002_R2_001.fastq.gz
# 116YVEH_S7_L001_I1_001.fastq.gz
# 116YVEH_S7_L001_R1_001.fastq.gz
# 116YVEH_S7_L001_R2_001.fastq.gz
# 116YVEH_S7_L002_I1_001.fastq.gz
# 116YVEH_S7_L002_R1_001.fastq.gz
# 116YVEH_S7_L002_R2_001.fastq.gz
# 011AVEH_S4_L001_I1_001.fastq.gz
# 011AVEH_S4_L001_R1_001.fastq.gz
# 011AVEH_S4_L001_R2_001.fastq.gz
# 011AVEH_S4_L002_I1_001.fastq.gz
# 011AVEH_S4_L002_R1_001.fastq.gz
# 011AVEH_S4_L002_R2_001.fastq.gz
# 023AVEH_S5_L001_I1_001.fastq.gz
# 023AVEH_S5_L001_R1_001.fastq.gz
# 023AVEH_S5_L001_R2_001.fastq.gz
# 023AVEH_S5_L002_I1_001.fastq.gz
# 023AVEH_S5_L002_R1_001.fastq.gz
# 023AVEH_S5_L002_R2_001.fastq.gz
# 026AVEH_S6_L001_I1_001.fastq.gz
# 026AVEH_S6_L001_R1_001.fastq.gz
# 026AVEH_S6_L001_R2_001.fastq.gz
# 026AVEH_S6_L002_I1_001.fastq.gz
# 026AVEH_S6_L002_R1_001.fastq.gz
# 026AVEH_S6_L002_R2_001.fastq.gz
# 029AVEH_S3_L001_I1_001.fastq.gz
# 029AVEH_S3_L001_R1_001.fastq.gz
# 029AVEH_S3_L001_R2_001.fastq.gz
# 029AVEH_S3_L002_I1_001.fastq.gz
# 029AVEH_S3_L002_R1_001.fastq.gz
# 029AVEH_S3_L002_R2_001.fastq.gz
# 033AVEH_S2_L001_I1_001.fastq.gz
# 033AVEH_S2_L001_R1_001.fastq.gz
# 033AVEH_S2_L001_R2_001.fastq.gz
# 033AVEH_S2_L002_I1_001.fastq.gz
# 033AVEH_S2_L002_R1_001.fastq.gz
# 033AVEH_S2_L002_R2_001.fastq.gz
# 017AIMP_S1_L001_I1_001.fastq.gz
# 017AIMP_S1_L001_R1_001.fastq.gz
# 017AIMP_S1_L001_R2_001.fastq.gz
# 017AIMP_S1_L002_I1_001.fastq.gz
# 017AIMP_S1_L002_R1_001.fastq.gz
# 017AIMP_S1_L002_R2_001.fastq.gz
# 025AIMP_S8_L001_I1_001.fastq.gz
# 025AIMP_S8_L001_R1_001.fastq.gz
# 025AIMP_S8_L001_R2_001.fastq.gz
# 025AIMP_S8_L002_I1_001.fastq.gz
# 025AIMP_S8_L002_R1_001.fastq.gz
# 025AIMP_S8_L002_R2_001.fastq.gz
# 031AIMP_S4_L001_I1_001.fastq.gz
# 031AIMP_S4_L001_R1_001.fastq.gz
# 031AIMP_S4_L001_R2_001.fastq.gz
# 031AIMP_S4_L002_I1_001.fastq.gz
# 031AIMP_S4_L002_R1_001.fastq.gz
# 031AIMP_S4_L002_R2_001.fastq.gz
# 045AIMP_S6_L001_I1_001.fastq.gz
# 045AIMP_S6_L001_R1_001.fastq.gz
# 045AIMP_S6_L001_R2_001.fastq.gz
# 045AIMP_S6_L002_I1_001.fastq.gz
# 045AIMP_S6_L002_R1_001.fastq.gz
# 045AIMP_S6_L002_R2_001.fastq.gz
# 046AIMP_S7_L001_I1_001.fastq.gz
# 046AIMP_S7_L001_R1_001.fastq.gz
# 046AIMP_S7_L001_R2_001.fastq.gz
# 046AIMP_S7_L002_I1_001.fastq.gz
# 046AIMP_S7_L002_R1_001.fastq.gz
# 046AIMP_S7_L002_R2_001.fastq.gz
#
# Cellranger outs used:  
#
# 01_102YVEH_filtered_feature_bc_matrix.h5
# 02_104YVEH_filtered_feature_bc_matrix.h5
# 03_106YVEH_filtered_feature_bc_matrix.h5
# 04_111YVEH_filtered_feature_bc_matrix.h5
# 05_116YVEH_filtered_feature_bc_matrix.h5
# 06_011AVEH_filtered_feature_bc_matrix.h5
# 07_023AVEH_filtered_feature_bc_matrix.h5
# 08_026AVEH_filtered_feature_bc_matrix.h5
# 09_029AVEH_filtered_feature_bc_matrix.h5
# 10_033AVEH_filtered_feature_bc_matrix.h5
# 11_017AIMP_filtered_feature_bc_matrix.h5
# 12_025AIMP_filtered_feature_bc_matrix.h5
# 13_031AIMP_filtered_feature_bc_matrix.h5
# 14_045AIMP_filtered_feature_bc_matrix.h5
# 15_046AIMP_filtered_feature_bc_matrix.h5
#
##

    #### Proteomics Data ####
##
# Proteomics Data (GitHub https://github.com/shaughnmb/2024_moser_et_al)
# 
# Raw data and calculations in excel format:  
#
#   moser_et_al_2024_iMP_Plasma_Proteomics_Raw_and_Calculations.xlsx
#
# Z Score spreadsheet from Raw and Calculations excel with sample order manually curated:
# 
#   moser_et_al_2024_iMP_Plasma_Proteomics_ZScore_YV_AV_AI.csv
#
# DEPs reformatted and filtered with outliers removed for volcano plots:
#
#   moser_et_al_2024_iMP_Plasma_Proteomics_For_Volcano.csv
#
##

#### SETUP ####
  #### General and Path variables ####

date <- "20240529"
project <- "iMPs_15_samples"
datadir <- "/Users/bells/Library/CloudStorage/Box-Box/GEO Submissions/GEO_Submission_Moser_et_al_2024/20240529_iMPs_15_samples"
sourcedir <- "/Users/bells/Library/CloudStorage/Box-Box/GEO Submissions/GEO_Submission_Moser_et_al_2024/20240529_iMPs_15_samples.h5"

  #### Load required packages ####

library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(irlba)
library(gridExtra)
library(scCustomize)
library(nichenetr)
library(EnhancedVolcano)
library(Vennerable)
library(glmnet)
library(ggpubr)
library(vip)

  #### Set working dir and save session info ####
setwd(datadir)

sessionInfo()

sink(paste0(date,"_",project,"_devtools_sessionInfo.txt"))
devtools::session_info()
sink()

sink(paste0(date,"_",project,"_sessionInfo.txt"))
sessionInfo()
sink()

  #### Functions for saving images ####
hirestiff <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 600,
    limitsize = TRUE,
    bg = "white"
  )
}

lowrestiff <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 100,
    limitsize = TRUE,
    bg = "white"
  )
}

hirestiffsquare <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 600,
    limitsize = TRUE,
    bg = "white",
    width = 9,
    height = 9,
    units = c("in")
  )
}

savepdf <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "pdf",
    scale = 1,
    dpi = 600,
    limitsize = TRUE
  )
}

savepdfsquare <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "pdf",
    scale = 1,
    dpi = 600,
    width = 9,
    height = 9,
    limitsize = TRUE
  )
}
  #### restore prefixes ####

# must enter at minimum the date, project, and datadir variables

# after determining original clustering, removing poor clusters, reclustering,
# and labeling cell types:
#
# save(datadir, sourcedir, date,
#      meaningful.PCs, prefixPC, prefixPCres,
#      project, res, ident_res,
#      file = paste0(datadir,"/",date,"_",project,"_reclust_","prefixes.rdata"))

load(paste0(datadir,"/",date,"_",project,"_reclust_","prefixes.rdata"))

# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))
# 
seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_cell_types_labeled.rds"))

  #### restore color pallets and functions #### 

# save(aging_palette, cell_cols, hirestiff,
#      lowrestiff, hirestiffsquare, savepdf, savepdfsquare,
#      file = paste0(datadir,"/",date,"_",project,"_reclust_","colors_functions.rdata"))

load(paste0(datadir,"/",date,"_",project,"_reclust_","colors_functions.rdata"))

  #### set color scheme ####

aging_palette <- c("#565DFF","#E13970","#175F12")

#### LOAD DATA ####
  #### Load data from Cellranger output ####

sample_list <- c(
  "01_102YVEH_filtered_feature_bc_matrix.h5",
  "02_104YVEH_filtered_feature_bc_matrix.h5",
  "03_106YVEH_filtered_feature_bc_matrix.h5",
  "04_111YVEH_filtered_feature_bc_matrix.h5",
  "05_116YVEH_filtered_feature_bc_matrix.h5",
  "06_011AVEH_filtered_feature_bc_matrix.h5",
  "07_023AVEH_filtered_feature_bc_matrix.h5",
  "08_026AVEH_filtered_feature_bc_matrix.h5",
  "09_029AVEH_filtered_feature_bc_matrix.h5",
  "10_033AVEH_filtered_feature_bc_matrix.h5",
  "11_017AIMP_filtered_feature_bc_matrix.h5",
  "12_025AIMP_filtered_feature_bc_matrix.h5",
  "13_031AIMP_filtered_feature_bc_matrix.h5",
  "14_045AIMP_filtered_feature_bc_matrix.h5",
  "15_046AIMP_filtered_feature_bc_matrix.h5"
)

# make path_list
filelist <- list.files(path = sourcedir)
filelist <- as.data.frame(filelist, files = list.files(path = sourcedir))
path_list <- paste0(sourcedir,"/",filelist[,1])

path_list

# check that all files exist
all(file.exists(path_list))

# check that the sample list and path list are in the same order
for (i in 1:length(sample_list)) {
  TF <- grepl(sample_list[i],path_list[i])
  print(paste0(sample_list[i]," ",TF))
}

# Pre-allocate a list to store the data
data_list <- vector("list", length = length(sample_list))

# Use lapply to read in cellranger output and assign to the list
data_list <- lapply(seq_along(sample_list), function(i) {
  # read in cellranger output
  cellranger <- Read10X_h5(path_list[i])
  return(cellranger)
})

updated_samps <- c(
  "YV_102",
  "YV_104",
  "YV_106",
  "YV_111",
  "YV_116",
  "AV_011",
  "AV_023",
  "AV_026",
  "AV_029",
  "AV_033",
  "AI_017",
  "AI_025",
  "AI_031",
  "AI_045",
  "AI_046"
  )

names(data_list) <- updated_samps

  #### Seurat Object ####

# Initialize the Seurat object with the raw (non-normalized data).
# Keep all genes expressed in >= 1 cell

# Preallocate a list to store the Seurat objects
seurat_list <- vector("list", length = length(data_list))

# Use lapply to create Seurat objects and assign to the list
seurat_list <- lapply(seq_along(data_list), function(i) {
  s.data <- CreateSeuratObject(counts = data_list[[i]],
                               project = updated_samps[i],
                               min.cells = 1,
                               min.features = 0)
  return(s.data)
})

names(seurat_list) <- updated_samps

seurat_list

  #### Add metadata ####

# Create lists for metadata attributes

Age <- c(
  "YV_102" = "Young",
  "YV_104" = "Young",
  "YV_106" = "Young",
  "YV_111" = "Young",
  "YV_116" = "Young",
  "AV_011" = "Aging",
  "AV_023" = "Aging",
  "AV_026" = "Aging",
  "AV_029" = "Aging",
  "AV_033" = "Aging",
  "AI_017" = "Aging",
  "AI_025" = "Aging",
  "AI_031" = "Aging",
  "AI_045" = "Aging",
  "AI_046" = "Aging"
)

Treatment <- c(
  "YV_102" = "Veh",
  "YV_104" = "Veh",
  "YV_106" = "Veh",
  "YV_111" = "Veh",
  "YV_116" = "Veh",
  "AV_011" = "Veh",
  "AV_023" = "Veh",
  "AV_026" = "Veh",
  "AV_029" = "Veh",
  "AV_033" = "Veh",
  "AI_017" = "iMPs",
  "AI_025" = "iMPs",
  "AI_031" = "iMPs",
  "AI_045" = "iMPs",
  "AI_046" = "iMPs"
)

Batch <- c(
  "YV_102" = "Batch_2",
  "YV_104" = "Batch_2",
  "YV_106" = "Batch_1",
  "YV_111" = "Batch_2",
  "YV_116" = "Batch_1",
  "AV_011" = "Batch_2",
  "AV_023" = "Batch_2",
  "AV_026" = "Batch_2",
  "AV_029" = "Batch_1",
  "AV_033" = "Batch_1",
  "AI_017" = "Batch_1",
  "AI_025" = "Batch_2",
  "AI_031" = "Batch_1",
  "AI_045" = "Batch_1",
  "AI_046" = "Batch_2"
)

orig.ident <- c(
  "YV_102" = "YV_102",
  "YV_104" = "YV_104",
  "YV_106" = "YV_106",
  "YV_111" = "YV_111",
  "YV_116" = "YV_116",
  "AV_011" = "AV_011",
  "AV_023" = "AV_023",
  "AV_026" = "AV_026",
  "AV_029" = "AV_029",
  "AV_033" = "AV_033",
  "AI_017" = "AI_017",
  "AI_025" = "AI_025",
  "AI_031" = "AI_031",
  "AI_045" = "AI_045",
  "AI_046" = "AI_046"
)

# Combine metadata into a data frame for easier processing
metadata_df <- data.frame(Age, Treatment, Batch, orig.ident)

# Use lapply to assign metadata to Seurat objects
seurat_list <- lapply(names(seurat_list), function(key) {
  seurat_obj <- seurat_list[[key]]
  metadata_row <- metadata_df[key, ]
  
  seurat_obj$Age <- metadata_row$Age
  seurat_obj$Treatment <- metadata_row$Treatment
  seurat_obj$Batch <- metadata_row$Batch
  seurat_obj$Type <- paste(metadata_row$Age, metadata_row$Treatment, sep = " ")
  seurat_obj$orig.ident <- metadata_row$orig.ident
  
  return(seurat_obj)
})

  #### Merge all data into a single seurat object ####

names(seurat_list)
length(names(seurat_list))

tenx1 <- merge(seurat_list[[1]], 
               y = c(
                 seurat_list[[2]],
                 seurat_list[[3]],
                 seurat_list[[4]],
                 seurat_list[[5]],
                 seurat_list[[6]],
                 seurat_list[[7]],
                 seurat_list[[8]],
                 seurat_list[[9]],
                 seurat_list[[10]],
                 seurat_list[[11]],
                 seurat_list[[12]],
                 seurat_list[[13]],
                 seurat_list[[14]],
                 seurat_list[[15]]
                 ),
add.cell.ids = c(
  names(seurat_list)[1],
  names(seurat_list)[2],
  names(seurat_list)[3],
  names(seurat_list)[4],
  names(seurat_list)[5],
  names(seurat_list)[6],
  names(seurat_list)[7],
  names(seurat_list)[8],
  names(seurat_list)[9],
  names(seurat_list)[10],
  names(seurat_list)[11],
  names(seurat_list)[12],
  names(seurat_list)[13],
  names(seurat_list)[14],
  names(seurat_list)[15]
  ),
project = project)

tenx1

rm(data_list,
   filelist,
   metadata_df,
   seurat_list,
   Batch,
   i,
   orig.ident,
   path_list,
   sample_list,
   TF,
   Treatment,
   updated_samps)

# save merged seurat object 
# saveRDS(tenx1, file = paste0("./",date,"_",project,"_merged_seurat_prefilter.rds"))

# tenx1 <- read_rds(file = paste0("./",date,"_",project,"_merged_seurat_prefilter.rds"))

#### QC FILTERING ####
  #### QC Setup ####
# copy tenx1 to a temporary data object just in case 
data <- tenx1

# check to make sure pattern is correct and not grepping other genes

grep("^mt-",rownames(data),value = T) # mitochondrial genes
grep("^Rp[sl]",rownames(data),value = T) # ribosomal genes
grep("^Mrp[sl]",rownames(data),value = T) # mitochondrial-ribosomal genes
grep("^Hb[abq]",rownames(data),value = T) # hemoglobin genes

# Add gene set info to metadata
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^Rp[sl]")
data[["percent.mt.ribo"]] <- PercentageFeatureSet(data, pattern = "^Mrp[sl]")
data[["percent.hb"]] <- PercentageFeatureSet(data, pattern = "^Hb[abq]")

  #### Prefiltered Plots ####
vp <- VlnPlot(data, features = "nCount_RNA", pt.size = 0) +
  NoLegend() +
  ggtitle("nCount_RNA prefiltered")
vp
pdf(paste0("./",date,"_",project,"_prefilter_","nCount_RNA",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0) +
  NoLegend() +
  ggtitle("nFeature_RNA prefiltered")
vp
pdf(paste0("./",date,"_",project,"_prefilter_","nFeature_RNA",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.mt", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.mt prefiltered")
vp
pdf(paste0("./",date,"_",project,"_prefilter_","percent.mt",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.ribo", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.ribo prefiltered")
vp
pdf(paste0("./",date,"_",project,"_prefilter_","percent.ribo",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.mt.ribo", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.mt.ribo prefiltered")
vp
pdf(paste0("./",date,"_",project,"_prefilter_","percent.mt.ribo",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.hb", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.hb prefiltered")
vp
pdf(paste0("./",date,"_",project,"_prefilter_","percent.hb",".pdf"))
print(vp)
dev.off()

p1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") 
p2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
p1 + p2

pdf(paste0("./",date,"_",project,"_prefilter_","scatterplots",".pdf"))
print(p1 + p2)
dev.off()

  #### Additional Prefiltered QC metrics and plots ####
qc.metrics <- as_tibble(data[[c("nCount_RNA",
                                "nFeature_RNA",
                                "percent.mt",
                                "percent.ribo",
                                "percent.mt.ribo",
                                "percent.hb")]],
                        rownames="Cell.Barcode")

sp <- qc.metrics %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Mitochondial DNA Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.mt_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  arrange(percent.ribo) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.ribo)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Ribosomal DNA Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.ribo_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  arrange(percent.mt.ribo) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt.ribo)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Mito-Ribosomal DNA Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.mt.ribo_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  arrange(percent.hb) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.hb)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Hemoglobin DNA Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.hb_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.mt)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Mitochondria Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_Dist_mito_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.ribo)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Ribosomes Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_Dist_ribo_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.mt.ribo)) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Mito-Ribosomes Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_Dist_mt_ribo_prefilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.hb)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Hemoglobin Prefiltered")
sp
pdf(paste0("./",date,"_",project,"_Dist_hb_prefilter",".pdf"))
print(sp)
dev.off()

  #### Add Z-scores  ####
data[["nUMI.z"]] <- scale(data$nCount_RNA)
data[["nGene.z"]] <- scale(data$nFeature_RNA)
data[["percent.mt.z"]] <- scale(data$percent.mt)
data[["percent.ribo.z"]] <- scale(data$percent.ribo)
data[["percent.mt.ribo.z"]] <- scale(data$percent.mt.ribo)
data[["percent.hb.z"]] <- scale(data$percent.hb)

  #### Prefiltered QC data  ####
length(data@meta.data$orig.ident)
mean(data@meta.data$percent.mt)
mean(data@meta.data$percent.ribo)
mean(data@meta.data$percent.mt.ribo)
mean(data@meta.data$percent.hb)
mean(data@meta.data$nCount_RNA)
median(data@meta.data$nCount_RNA)
mean(data@meta.data$nFeature_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

# save prefilter QC data
sink(paste0("./",date,"_",project,"_prefilter_QC_metrics.txt"))
cat("length")
length(data@meta.data$orig.ident)
cat("mean percent mito")
mean(data@meta.data$percent.mt)
cat("mean percent ribo")
mean(data@meta.data$percent.ribo)
cat("mean percent mito-ribo")
mean(data@meta.data$percent.mt.ribo)
cat("mean hemoglobin")
mean(data@meta.data$percent.hb)
cat("mean percent counts")
mean(data@meta.data$nCount_RNA)
cat("median counts")
median(data@meta.data$nCount_RNA)
cat("mean features")
mean(data@meta.data$nFeature_RNA)
cat("median features")
median(data@meta.data$nFeature_RNA)
cat("max counts")
max(data@meta.data$nCount_RNA)
sink()

  #### Filter cells based on Z-score  ####
data <- subset(data, subset = percent.mt.z < 3)
data <- subset(data, subset = percent.ribo.z < 3)
data <- subset(data, subset = percent.mt.ribo.z < 3)
data <- subset(data, subset = percent.hb.z < 3)

  #### Postfiltered QC data ####
length(data@meta.data$orig.ident)
mean(data@meta.data$percent.mt)
mean(data@meta.data$percent.ribo)
mean(data@meta.data$percent.hb)
mean(data@meta.data$nCount_RNA)
median(data@meta.data$nCount_RNA)
mean(data@meta.data$nFeature_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

# save postfilter QC data
sink(paste0("./",date,"_",project,"_postfilter_QC_metrics.txt"))
cat("length")
length(data@meta.data$orig.ident)
cat("mean percent mito")
mean(data@meta.data$percent.mt)
cat("mean percent ribo")
mean(data@meta.data$percent.ribo)
cat("mean percent mito-ribo")
mean(data@meta.data$percent.mt.ribo)
cat("mean hemoglobin")
mean(data@meta.data$percent.hb)
cat("mean percent counts")
mean(data@meta.data$nCount_RNA)
cat("median counts")
median(data@meta.data$nCount_RNA)
cat("mean features")
mean(data@meta.data$nFeature_RNA)
cat("median features")
median(data@meta.data$nFeature_RNA)
cat("max counts")
max(data@meta.data$nCount_RNA)
sink()

  #### Postfiltered Plots ####
vp <- VlnPlot(data, features = "nCount_RNA", pt.size = 0, group.by = "orig.ident")  + 
  NoLegend() + 
  ggtitle("nCount_RNA filtered")
vp
pdf(paste0("./",date,"_",project,"_filtered_","nCount_RNA",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") + 
  NoLegend() + 
  ggtitle("nFeature_RNA filtered")
vp
pdf(paste0("./",date,"_",project,"_filtered_","nFeature_RNA",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.mt", pt.size = 0, group.by = "orig.ident") +
  NoLegend() +
  ggtitle("percent.mt filtered")
vp
pdf(paste0("./",date,"_",project,"_filtered_","percent.mt",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.ribo", pt.size = 0, group.by = "orig.ident") +
  NoLegend() +
  ggtitle("percent.ribo filtered")
vp
pdf(paste0("./",date,"_",project,"_filtered_","percent.ribo",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.mt.ribo", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.mt.ribo filtered")
vp
pdf(paste0("./",date,"_",project,"_filtered_","percent.mt.ribo",".pdf"))
print(vp)
dev.off()

vp <- VlnPlot(data, features = "percent.hb", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.hb filtered")
vp
pdf(paste0("./",date,"_",project,"_filtered_","percent.hb",".pdf"))
print(vp)
dev.off()

p1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") 
p2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
p1 + p2

pdf(paste0("./",date,"_",project,"_filtered_","scatterplots",".pdf"))
print(p1 + p2)
dev.off()

  #### Additional postfiltered QC metrics and plots ####
qc.metrics <- as_tibble(data[[c("nCount_RNA",
                                "nFeature_RNA",
                                "percent.mt",
                                "percent.ribo",
                                "percent.mt.ribo",
                                "percent.hb")]],
                        rownames="Cell.Barcode")

sp <- qc.metrics %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Mitochondial DNA Postfilter")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.mt_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  arrange(percent.ribo) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.ribo)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Ribosomal DNA Postfilter")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.ribo_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  arrange(percent.mt.ribo) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt.ribo)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Mito-Ribosomal DNA Postfilter")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.mt.ribo_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  arrange(percent.hb) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.hb)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics - percent Hemoglobin DNA Postfilter")
sp
pdf(paste0("./",date,"_",project,"_QC_Metrics_pct.hb_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.mt)) + 
  geom_histogram(binwidth = 0.25, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Mitochondria Postfilter")
sp
pdf(paste0("./",date,"_",project,"_Dist_mito_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.ribo)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Ribosomes Postfilter")
sp
pdf(paste0("./",date,"_",project,"_Dist_ribo_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.mt.ribo)) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Mito-Ribosomes Postfilter")
sp
pdf(paste0("./",date,"_",project,"_Dist_mt_ribo_postfilter",".pdf"))
print(sp)
dev.off()

sp <- qc.metrics %>%
  ggplot(aes(percent.hb)) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage of reads from Hemoglobin Postfilter")
sp
pdf(paste0("./",date,"_",project,"_Dist_hb_postfilter",".pdf"))
print(sp)
dev.off()

  #### post QC cleanup ####

tenx1
data

tenx1 <- data

rm(vp,sp,p1,p2,
   qc.metrics,data)

#### SCTRANSFORM, PCA, UMAP, CLUSTERING, and NORMALIZE ####
  #### set factor levels ####
seurat_sct <- tenx1

seurat_sct$orig.ident <- factor(x = seurat_sct$orig.ident, 
                                levels = c("YV_102","YV_104","YV_106","YV_111","YV_116",
                                           "AV_011","AV_023","AV_026","AV_029","AV_033",
                                           "AI_017","AI_025","AI_031","AI_045","AI_046"))
seurat_sct$Type <- factor(x = seurat_sct$Type, 
                          levels = c("Young Veh","Aging Veh","Aging iMPs"))
seurat_sct$Treatment <- factor(x = seurat_sct$Treatment, 
                               levels = c("Veh","iMPs"))
seurat_sct$Age <- factor(x = seurat_sct$Age, 
                         levels = c("Young","Aging"))
seurat_sct$Batch <- factor(x = seurat_sct$Batch, 
                         levels = c("Batch_1","Batch_2"))

rm(tenx1)

  #### join layers ####
seurat_sct[["RNA"]] <- JoinLayers(seurat_sct[["RNA"]]) # for seurat v5 + only

  #### SCTransform ####

seurat_sct <- SCTransform(seurat_sct,
                          vst.flavor = "v1", # for seurat v5 + only
                          ncells = 10595, #total cells / 10 (105951/10 = 10595)
                          vars.to.regress = c("percent.mt",
                                              "percent.ribo",
                                              "percent.mt.ribo",
                                              "percent.hb"),
                          verbose = TRUE)

# saveRDS(seurat_sct, file = paste0("./",date,"_",project,"_post_sct_pre_pca.rds"))

# seurat_sct <- readRDS(file = paste0("./",date,"_",project,"_post_sct_pre_pca.rds"))

  #### Run PCA ####

# this performs PCA on the seurat object
seurat_sct <- RunPCA(seurat_sct, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(seurat_sct@reductions$pca@cell.embeddings)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(seurat_sct@reductions$pca@feature.loadings)

# calculate eigenvalues for arrays

# generate squares of all sample coordinates
sq.xx.coord <- as.data.frame(xx.coord^2)
# create empty list for eigenvalues first
eig <- c()
# calculate the eigenvalue for each PC in sq.xx.coord by taking the sqrt of the sum of squares
for(i in 1:ncol(sq.xx.coord))
  eig[i] = sqrt(sum(sq.xx.coord[,i]))
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if they all contribute equally to the total variance
expected.contribution <- sum.eig/(length(xx.coord)-1)
# return the number of principal components with an eigenvalue greater than expected by equal variance
meaningful.PCs <- sum(eig > expected.contribution)

# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for(i in 1:length(eig))
  eig.percent[i] = 100*eig[i]/sum.eig
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for(i in 1:length(eig))
  if(i == 1) scree[i] = eig.percent[i] else scree[i] = scree[i-1] + eig.percent[i]

# create data frame for eigenvalue summaries
eigenvalues <- data.frame("PC" = colnames(xx.coord), "eig" = eig, "percent" = eig.percent, "scree" = scree)

# plot scree values
plot(eigenvalues$percent, ylim = c(0,100), type = "S", xlab = "PC", ylab = "Percent of variance",
     main = paste0(date,"_",project," scree plot all samples PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = meaningful.PCs, col = "blue")
text(meaningful.PCs, cut.off, label = paste("cutoff PC",meaningful.PCs),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0("./",date,"_",project,"_scree_plot.pdf"))
dev.off()

# meaningful.PCs <- 13

rm(eigenvalues,sq.xx.coord,xx.coord,xx.gload,cut.off,eig,
   eig.percent,expected.contribution,i,scree,sum.eig)

  #### Run UMAP, and look at UMAP plots ####

# update prefixed variable
prefixPC <- paste0("./",date,"_",project,"_",meaningful.PCs,"PCs")

# run UMAP
seurat_sct <- RunUMAP(seurat_sct, reduction = "pca", dims = 1:meaningful.PCs, verbose = TRUE)

# UMAP plot by sample name ("orig.ident")
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "orig.ident",
        raster = F)
hirestiff(paste0(prefixPC,"_UMAP_by_orig.ident_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_orig.ident_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "orig.ident", 
        group.by = "orig.ident", 
        combine = TRUE,
        raster = F) + 
  NoLegend()

hirestiff(paste0(prefixPC,"_UMAP_split_by_orig.ident_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_orig.ident_lowres.tiff"))

Idents(seurat_sct) <- "orig.ident"

for (i in 1:length(levels(seurat_sct))){
  p <- (DimPlot(seurat_sct, 
          reduction = "umap", 
          label = FALSE, 
          pt.size = .25,
          raster = F, 
          cols = Hue_Pal(length(levels(seurat_sct)))[i],
          cells = c(WhichCells(seurat_sct, idents = levels(seurat_sct)[i]))) + 
    NoLegend() +
    ggtitle(paste0(levels(seurat_sct)[i])))
  
  print(p)
hirestiff(paste0(prefixPC,"_UMAP_split_by_orig.ident_",levels(seurat_sct)[i],"_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_orig.ident_",levels(seurat_sct)[i],"_lowres.tiff"))
}

# by Type 

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "Type",
        cols = aging_palette,
        raster = F)
hirestiff(paste0(prefixPC,"_UMAP_by_type_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_type_lowres.tiff"))

DimPlot(seurat_sct, reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "Type", 
        group.by = "Type",
        cols = aging_palette,
        raster = F) + NoLegend()
hirestiff(paste0(prefixPC,"_UMAP_split_by_type_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_type_lowres.tiff"))

Idents(seurat_sct) <- "Type"

for (i in 1:length(levels(seurat_sct))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .25, 
                cols = aging_palette[i],
                raster = F,
                cells = c(WhichCells(seurat_sct, idents = levels(seurat_sct)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct)[i])))
  
  print(p)
  hirestiff(paste0(prefixPC,"_UMAP_split_by_type_",levels(seurat_sct)[i],"_hires.tiff"))
  lowrestiff(paste0(prefixPC,"_UMAP_split_by_type_",levels(seurat_sct)[i],"_lowres.tiff"))
}

# by Batch 

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "Batch",
        cols = aging_palette,
        raster = F)
hirestiff(paste0(prefixPC,"_UMAP_by_Batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_Batch_lowres.tiff"))

DimPlot(seurat_sct, reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "Batch", 
        group.by = "Batch",
        cols = aging_palette,
        raster = F) + NoLegend()
hirestiff(paste0(prefixPC,"_UMAP_split_by_Batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_Batch_lowres.tiff"))

Idents(seurat_sct) <- "Batch"

for (i in 1:length(levels(seurat_sct))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .25, 
                cols = aging_palette[i],
                raster = F,
                cells = c(WhichCells(seurat_sct, idents = levels(seurat_sct)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct)[i])))
  
  print(p)
  hirestiff(paste0(prefixPC,"_UMAP_split_by_Batch_",levels(seurat_sct)[i],"_hires.tiff"))
  lowrestiff(paste0(prefixPC,"_UMAP_split_by_Batch_",levels(seurat_sct)[i],"_lowres.tiff"))
}

  #### Normalize RNA, find variable features, scale data for visualization ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(seurat_sct) <- "RNA"

# Normalize, find variable features, scale data 
seurat_sct <- NormalizeData(seurat_sct, verbose = FALSE)
seurat_sct <- FindVariableFeatures(seurat_sct)
all.genes <- rownames(seurat_sct)
seurat_sct <- ScaleData(seurat_sct, features = all.genes)

write.csv(all.genes, 
          file = "all_hippocampus_genes.csv", 
          quote = T, 
          row.names = F)

  #### Add Cell Cycle info ####

## convert human cc genes to mouse genes
mm.cc.genes <- list()

mm.cc.genes$s.genes <- convert_human_to_mouse_symbols(cc.genes.updated.2019$s.genes)
mm.cc.genes$g2m.genes <- nichenetr::convert_human_to_mouse_symbols(cc.genes.updated.2019$g2m.genes)

# add cell cycle module scores
seurat_sct <- CellCycleScoring(seurat_sct,
                               s.features = mm.cc.genes$s.genes,
                               g2m.features = mm.cc.genes$g2m.genes,
                               set.ident = TRUE)

as_tibble(seurat_sct[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "Phase",
        raster = F)
hirestiff(paste0(prefixPC,"_UMAP_by_phase_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_phase_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, pt.size = .25, 
        split.by = "Type", 
        group.by = "Phase",
        raster = F)
hirestiff(paste0(prefixPC,"_UMAP_by_phase_split_by_type_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_phase_split_by_type_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, pt.size = .25, 
        split.by = "Batch", 
        group.by = "Phase",
        raster = F)
hirestiff(paste0(prefixPC,"_UMAP_by_phase_split_by_batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_phase_split_by_batch_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, pt.size = .25, 
        split.by = "Phase", 
        group.by = "Phase",
        raster = F)
hirestiff(paste0(prefixPC,"_UMAP_by_phase_split_by_phase_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_phase_split_by_phase_lowres.tiff"))

Idents(seurat_sct) <- "Phase"

for (i in 1:length(levels(seurat_sct))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .25,
                raster = F, 
                cols = Hue_Pal(length(levels(seurat_sct)))[i],
                cells = c(WhichCells(seurat_sct, idents = levels(seurat_sct)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct)[i])))
  
  print(p)
  hirestiff(paste0(prefixPC,"_UMAP_split_by_phase_",levels(seurat_sct)[i],"_hires.tiff"))
  lowrestiff(paste0(prefixPC,"_UMAP_split_by_phase_",levels(seurat_sct)[i],"_lowres.tiff"))
}

rm(mm.cc.genes,p,i)

# saveRDS(seurat_sct, file = paste0(prefixPC,"_seurat_post_umap.rds"))

# seurat_sct <- read_rds(file = paste0(prefixPC,"_seurat_post_umap.rds"))

  #### Clustering and Resolution ####

DefaultAssay(seurat_sct) <- "SCT"

# Determine the K-nearest neighbor graph
seurat_sct <- FindNeighbors(object = seurat_sct, 
                            reduction = "pca", 
                            dims = 1:meaningful.PCs)

# Determine the clusters                              

seurat_sct <- FindClusters(object = seurat_sct,
                           resolution = c(0.6))

res <- "_res.0.6"
ident_res <- paste0("SCT_snn",res)

Idents(seurat_sct) <- ident_res
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        raster = F) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

#update prefix
prefixPCres <- paste0(prefixPC,res)

# Plot the UMAP
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        raster = F)  + NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","by","cluster","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","by","cluster","lowres.tiff",sep = "_"))

# UMAP of cells in each cluster by type without cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        split.by = "Type",
        raster = F)  + NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","split_by","type","without_cluster_labels","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","split_by","type","without_cluster_labels","lowres.tiff",sep = "_"))

for (i in 1:length(levels(seurat_sct$Type))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .4,
                raster = F, 
                cells = c(WhichCells(seurat_sct, expression = Type == levels(seurat_sct$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","split_by","type","no_cluster_labels",
                  levels(seurat_sct$Type)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","split_by","type","no_cluster_labels",
                   levels(seurat_sct$Type)[i],"lowres.tiff",sep = "_"))
}

# UMAP of cells in each cluster by type with cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5, 
        split.by = "Type",
        raster = F)  + NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","split_by","type","with_cluster_labels","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","split_by","type","with_cluster_labels","lowres.tiff",sep = "_"))

for (i in 1:length(levels(seurat_sct$Type))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = TRUE, 
                pt.size = .4,
                raster = F,
                label.size = 5, 
                cells = c(WhichCells(seurat_sct, expression = Type == levels(seurat_sct$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","split_by","type","with_cluster_labels",
                  levels(seurat_sct$Type)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","split_by","type","with_cluster_labels",
                   levels(seurat_sct$Type)[i],"lowres.tiff",sep = "_"))
}

# UMAP of cells in each cluster by batch with cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = F, 
        split.by = "Batch",
        raster = F,
        label.size = 5) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","split_by","batch","without_cluster_labels","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","split_by","batch","without_cluster_labels","lowres.tiff",sep = "_"))

rm(p)

  #### set RNA as active assay for visualization ####
DefaultAssay(seurat_sct) <- "RNA"

  #### Determine cell types ####

# 0 - Excitatory Neurons EXCN
VlnPlot(seurat_sct,
        features = c("Slc17a7","Neurod6","Camk2a","Nell2"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 0,1,3,8 - Excitatory Neurons")
hirestiff(paste(prefixPCres,"vln","Cluster","0_1_3_8","ExNeurons","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","0_1_3_8","ExNeurons","lowres.tiff", sep = "_"))

# 1 - EXCN

# 2 - DG Granule (which also express glut markers) GRN
VlnPlot(seurat_sct,
        features = c("Bdnf","Calb1","Calm1","Calm2","Calm3",
                     "Ppp3ca","Disc1","Gabrb1","Gabra4","Gria1",
                     "Gria2","Gria3","Grin1","Grin2b","Chrm1",
                     "Nrgn","Nrp2","Pcp4","Ncam1","Prox1",
                     "Vsnl1","Npy1r","Actn2","Grm1","Grm2",
                     "Grm3","Grm5","Creb1","Slc17a7"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 1 - DG Granule Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","2","DG_Granule","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","2","DG_Granule","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Bdnf","Calb1","Disc1","Gabrb1","Gabra4",
                     "Gria1","Gria2","Gria3","Grin1","Grin2b",
                     "Chrm1","Nrp2","Ncam1","Prox1","Npy1r",
                     "Actn2","Grm1","Grm2","Grm3","Grm5",
                     "Creb1","Slc17a7"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 2 - DG Granule Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","2","DG_Granule","refined_list","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","2","DG_Granule","refined_list","lowres.tiff", sep = "_"))

# 3 - EXCN

# 4 - PYR Pyramidal PYR
VlnPlot(seurat_sct,
        features = c("Amigo2","Disc1","Bdnf","Bok",
                     "Ppp3ca","Gabrb1","Gabra4","Gabra5",
                     "Gria1","Man1a","Neurod6","Chrm1","Chrm3",
                     "Nrp2","Pou3f1","Npy2r",
                     "Actn2","Grm1","Grm2","Grm3","Grm4",
                     "Grm5","Pcp4","Creb1"),
        stack = TRUE,
        flip = TRUE) +  
  NoLegend() + ggtitle("Clusters 4, 9, 10, 19 - Pyramidal Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","4_9_10_19","pyramidal","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","4_9_10_19","pyramidal","lowres.tiff", sep = "_"))

# CA2 cannonical
VlnPlot(seurat_sct,
        features = c("Accn1", #not found
                     "Adcy5",
                     "Cacng5",
                     "Ccdc3",
                     "Fam40b", #not found
                     "Gpr12",
                     "Gsto1",
                     "Map3k15",
                     "Ptpn5",
                     "Prss23",
                     "Rgs14",
                     "S100b",
                     "Sostdc1"),
        stack = TRUE,
        flip = TRUE,
        idents = c(4,9,10,19)) + 
  NoLegend()

# CA2 "new"
VlnPlot(seurat_sct,
        features = c("170024P16Rik", # not found
                     "Aldh1a1",
                     "Arg2",
                     "Crabp1",
                     "Ctsc",
                     "Dusp5",
                     "F2r",
                     "Fgf2",
                     "Fgf5",
                     "Glul",
                     "Gm13847", # not found
                     "Gm13921", # not found
                     "Gm20459", # not found
                     "Maob",
                     "Ntsr2",
                     "Plcb4",
                     "Pygo1",
                     "Rgs5",
                     "Scgn",
                     "Srgap2",
                     "Srl",
                     "Syce2",
                     "Tgfb1i1",
                     "Vcan",
                     "Vit"
        ),
        stack = TRUE,
        flip = TRUE,
        idents = c(4,9,10,19)) +  
  NoLegend()

# Broad Hippocampal populations (Cembrowski)
VlnPlot(seurat_sct,
        features = c("Prox1", # GRN
                     "Dkk3", # Non-GRN
                     "Calb2", # Mossy
                     "Ociad2", # All Pyrmidal
                     "Cacng5", # CA2 Pyr
                     "Fibcd1", # CA1 Pyr
                     "Pdzd2", # GRN dorsal
                     "Tox3", # GRN ventral
                     "Iyd", # CA3 dorsal
                     "Coch", # CA3 ventral
                     "Wfs1", # CA1 dorsal
                     "Dcn"), # CA1 ventral
        # idents = c(4,9,10,19),
        stack = TRUE,
        flip = TRUE
        ) +  
  NoLegend()

# 5 - Oligodendrocytes OLIG
VlnPlot(seurat_sct,
        features = c("Cnp","Cldn11","Mbp","Mog","Olig1","Sox10"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 5_6 - Oligodendrocytes")
hirestiff(paste(prefixPCres,"vln","Cluster","5_6","Oligodendrocytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","5_6","Oligodendrocytes","lowres.tiff", sep = "_"))

# 6 - Oligodendrocytes OLIG

# 7 - Inhibitory Neurons INHN
VlnPlot(seurat_sct,
        features = c("Gad1","Synpr","Gad2","Sst"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 7 - Inhibitory Neurons")
hirestiff(paste(prefixPCres,"vln","Cluster","7_12","Inhibitory","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","7_12","Inhibitory","lowres.tiff", sep = "_"))

# 8 - EXCN

# 9 - PYR

# 10 - PYR

# 11 - Astrocytes ASTR
VlnPlot(seurat_sct,
        features = c("Rgs20","Aqp4","Ntsr2","Gfap","S100b"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 11 - Astrocytes")
hirestiff(paste(prefixPCres,"vln","Cluster","11","Astrocytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","11","Astrocytes","lowres.tiff", sep = "_"))

# 12 - DG Basket (inhibitory via Zhong) BSKT
VlnPlot(seurat_sct,
        features = c("Dlx6os1","Maf","Erbb4","Gad2","Grik1",
                     "Nxph1","Slc6a1","Gad1","Gm13629","Dlx1",
                     "Dlx6","Col19a1","Utrn","Reln","Npas1",
                     "Sox6","Crhbp","Npy","Sst",
                     "Gabra1","Pvalb"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 12 - DG Basket Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","12","DG","Basket","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","12","DG","Basket","lowres.tiff", sep = "_"))

# 13 - Microglia MG
VlnPlot(seurat_sct,
        features = c("Itgam","Ptprc","Tmem119","Cx3cr1","P2ry12","C1qa"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 13 - Microglia")
hirestiff(paste(prefixPCres,"vln","Cluster","13","Microglia","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","13","Microglia","lowres.tiff", sep = "_"))

# 14 - OPC
VlnPlot(seurat_sct,
        features = c("Pdgfra","Cspg4","Sox10","Olig1","Olig2"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 14 - OPCs")
hirestiff(paste(prefixPCres,"vln","Cluster","14","OPCs","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","14","OPCs","lowres.tiff", sep = "_"))

# 15 - DG Mossy MOSS
VlnPlot(seurat_sct,
        features = c("Gad1","Gad2","Calb2",
                     "Gria2","Gria3","Pcp4",
                     "Gabra6","Pvalb","Nos1"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 15 - DG Mossy Cells (neg for Gabra6, Pvalb & Nos1)")
hirestiff(paste(prefixPCres,"vln","Cluster","15","DG","Mossy","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","15","DG","Mossy","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Gria2","Gria3","Pcp4"),
        stack = TRUE,
        flip = TRUE,
        idents = c("15"),
        split.by = "Type",
        cols = aging_palette) + 
  ggtitle("Cluster 15 - DG Mossy Cells split by Type")
lowrestiff(paste(prefixPCres,"vln","Cluster","15","DG","Mossy","by","type","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Gria2"),
        split.by = "Type",
        idents = c("15"),
        flip = TRUE,
        cols = aging_palette,
        pt.size = 1)
lowrestiff(paste(prefixPCres,"vln","Cluster","15","DG","Mossy","by","type","Gria2","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Gria3"),
        split.by = "Type",
        idents = c("15"),
        flip = TRUE,
        cols = aging_palette,
        pt.size = 1)
lowrestiff(paste(prefixPCres,"vln","Cluster","15","DG","Mossy","by","type","Gria3","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Pcp4"),
        split.by = "Type",
        idents = c("15"),
        flip = TRUE,
        cols = aging_palette,
        pt.size = 1)
lowrestiff(paste(prefixPCres,"vln","Cluster","15","DG","Mossy","by","type","Pcp4","lowres.tiff", sep = "_"))

# 16 - CPCa

# CPC Choroid Plexus Cells
VlnPlot(seurat_sct,features = c("Ttr","Atp5g1","Clic6","Ndufv3","Uqcrb",
                                "Hemk1","Atp5l","Ndufa4","Ppp1r1b","Cox7b",
                                "Atp5mpl","Atp5md","Ndufa1","Uqcrq","Cox6c",
                                "Cox7c","Atp5j","Atp5e","Uqcr10","Cox8b",
                                "Acad8","Uqcrh","Atp5j2","Uqcr11","Cox5b",
                                "Cox8a","Cox6b1","Ldhb"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 16 - Choroid Plexus Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","16","CPC","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","16","CPC","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Ttr","Clic6","Hemk1","Ppp1r1b",
                     "Acad8","Folr1","Prlr","Kl"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 16 - Choroid Plexus Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","16","CPC","refined_list","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","16","CPC","refined_list","lowres.tiff", sep = "_"))

# 17 - Vascular Leptomeningeal Cells VLMC
VlnPlot(seurat_sct,
        features = c("Pdgfra","Gja1","Slc6a13"),
        stack = TRUE,flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 17 - VLMC")
hirestiff(paste(prefixPCres,"vln","Cluster","17","VLMC","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","17","VLMC","lowres.tiff", sep = "_"))

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Pdgfra","Gja1","Slc6a13","Nes"), 
            order = TRUE,
            # min.cutoff = 'q10', 
            label = TRUE,
            label.size = 3,
            pt.size = 0.3,
            raster = F)

# 18 - CPCb

# 19 - PYR

# 20 - Pericytes PERI
VlnPlot(seurat_sct,
        features = c("Kcnj8","Cspg4","Pdgfrb","Rgs5"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 20 - Pericytes")
hirestiff(paste(prefixPCres,"vln","Cluster","20","Pericytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","20","Pericytes","lowres.tiff", sep = "_"))

# 21 - Activated Microglia MG-A
VlnPlot(seurat_sct,features = c("Itgam","Tmem119","Cx3cr1",
                                "P2ry12","Mrc1",
                                "Jun","Junb","Jund","C3"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 21 - Activated Microglia")
hirestiff(paste(prefixPCres,"vln","Cluster","21","Activated","Microglia","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","21","Activated","Microglia","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,features = c("Itgam","Tmem119","Cx3cr1",
                                "Tlr4","Fosb","Jun","Junb","Jund"),
        idents = c(13,21),
        split.by = "Type",
        stack = TRUE,
        flip = TRUE) + 
  ggtitle("Cluster 21 - Activated Microglia")
hirestiff(paste(prefixPCres,"vln","Cluster","13_21","Microglia","by_type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","13_21","Microglia","by_type","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Mrc1"),
        split.by = "Type",
        idents = c("13","21"),
        flip = TRUE,
        pt.size = 1)
lowrestiff(paste(prefixPCres,"vln","Cluster","21","Activated","Microglia","Mrc1","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Tlr4"),
        split.by = "Type",
        idents = c("13","21"),
        flip = TRUE,
        pt.size = 1)
lowrestiff(paste(prefixPCres,"vln","Cluster","21","Activated","Microglia","Tlr4","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Cd163"),
        split.by = "Type",
        idents = c("13","21"),
        flip = TRUE,
        pt.size = 1)
lowrestiff(paste(prefixPCres,"vln","Cluster","21","Activated","Microglia","Cd163","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,features = c("Tmem119"),split.by = "Type",idents = c("13","21"),flip = TRUE)
VlnPlot(seurat_sct,features = c("Jund"),split.by = "Type",idents = c("13","21"),flip = TRUE)
VlnPlot(seurat_sct,features = c("Mapk14"),split.by = "Type",idents = c("13","21"),flip = TRUE)

VlnPlot(seurat_sct,
        features = c("C3"),
        split.by = "Type",
        idents = c("13","21"),
        flip = TRUE,
        cols = aging_palette,
        pt.size = 1)

hirestiff(paste(prefixPCres,"vln","MG","C3","by","Type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","MG","C3","by","Type","lowres.tiff", sep = "_"))

# 22 - A1 Astrocytes ASTR-A1
VlnPlot(seurat_sct,
        features = c("Rgs20","Aqp4","Ntsr2","Gfap","S100b","C3"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 22 - A1 Astrocytes")

hirestiff(paste(prefixPCres,"vln","Cluster","22","A1","Astrocytes","res.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","22","A1","Astrocytes","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("C3"),
        split.by = "Type",
        idents = c("11","22"),
        flip = TRUE,
        cols = aging_palette,
        pt.size = 1)

hirestiff(paste(prefixPCres,"vln","astro","A1","C3","by","Type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","astro","A1","C3","by","Type","lowres.tiff", sep = "_"))

# 23 - multiple
# 24 - multiple
# 25 - multiple
# 26 - multiple

  #### Remove poor clusters ####

# remove clusters with more than one cell type marker (multiplets)
seurat_sct <- subset(seurat_sct, idents = c(23,24,25,26), invert = TRUE)

DimPlot(seurat_sct,
        reduction = "umap",
        label = TRUE,
        label.size = 4,
        pt.size = 0.5,
        raster = F) + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","confounding","clusters","removed","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","confounding","clusters","removed","lowres.tiff", sep = "_"))
#### redo SCTRANSFORM, PCA, UMAP, and NORMALIZE ####
  #### SCTransform ####

seurat_sct

seurat_sct <- SCTransform(seurat_sct,
                          vst.flavor = "v1", # only for seurat v5 +
                          ncells = 10561, #total cells / 10 (105612/10 = 10561)
                          vars.to.regress = c("percent.mt",
                                              "percent.ribo",
                                              "percent.mt.ribo",
                                              "percent.hb"),
                          verbose = TRUE)

# saveRDS(seurat_sct, file = paste0("./",date,"_",project,"_reclust_post_sct_pre_pca.rds"))

# seurat_sct <- readRDS(file = paste0("./",date,"_",project,"_reclust_post_sct_pre_pca.rds"))

  #### Run PCA ####

# this performs PCA on the seurat object
seurat_sct <- RunPCA(seurat_sct, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(seurat_sct@reductions$pca@cell.embeddings)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(seurat_sct@reductions$pca@feature.loadings)

# calculate eigenvalues for arrays

# generate squares of all sample coordinates
sq.xx.coord <- as.data.frame(xx.coord^2)
# create empty list for eigenvalues first
eig <- c()
# calculate the eigenvalue for each PC in sq.xx.coord by taking the sqrt of the sum of squares
for(i in 1:ncol(sq.xx.coord))
  eig[i] = sqrt(sum(sq.xx.coord[,i]))
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if they all contribute equally to the total variance
expected.contribution <- sum.eig/(length(xx.coord)-1)
# return the number of principal components with an eigenvalue greater than expected by equal variance
meaningful.PCs <- sum(eig > expected.contribution)

# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for(i in 1:length(eig))
  eig.percent[i] = 100*eig[i]/sum.eig
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for(i in 1:length(eig))
  if(i == 1) scree[i] = eig.percent[i] else scree[i] = scree[i-1] + eig.percent[i]

# create data frame for eigenvalue summaries
eigenvalues <- data.frame("PC" = colnames(xx.coord), "eig" = eig, "percent" = eig.percent, "scree" = scree)

# plot scree values
plot(eigenvalues$percent, ylim = c(0,100), type = "S", xlab = "PC", ylab = "Percent of variance",
     main = paste0(date,"_",project," scree plot all samples PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = meaningful.PCs, col = "blue")
text(meaningful.PCs, cut.off, label = paste("cutoff PC",meaningful.PCs),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0("./",date,"_",project,"_scree_plot.pdf"))
dev.off()

# meaningful.PCs <- 13

rm(eigenvalues,sq.xx.coord,xx.coord,xx.gload,cut.off,eig,
   eig.percent,expected.contribution,i,scree,sum.eig)

  #### Run UMAP, and look at UMAP plots ####

# update prefixed variable
prefixPC <- paste0("./",date,"_",project,"_reclust_",meaningful.PCs,"PCs")

# run UMAP
seurat_sct <- RunUMAP(seurat_sct, reduction = "pca", dims = 1:meaningful.PCs, verbose = TRUE)

# UMAP plot by sample name ("orig.ident")
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "orig.ident",
        raster = F)
hirestiff(paste0(prefixPC,"_UMAP_by_orig.ident_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_orig.ident_lowres.tiff"))

# by Type 

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "Type",
        cols = aging_palette,
        raster = F)
hirestiff(paste0(prefixPC,"_UMAP_by_type_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_type_lowres.tiff"))

DimPlot(seurat_sct, reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "Type", 
        group.by = "Type",
        cols = aging_palette,
        raster = F) + NoLegend()
hirestiff(paste0(prefixPC,"_UMAP_split_by_type_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_type_lowres.tiff"))

Idents(seurat_sct) <- "Type"

for (i in 1:length(levels(seurat_sct))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .25, 
                cols = aging_palette[i],
                raster = F,
                cells = c(WhichCells(seurat_sct, idents = levels(seurat_sct)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct)[i])))
  
  print(p)
  hirestiff(paste0(prefixPC,"_UMAP_split_by_type_",levels(seurat_sct)[i],"_hires.tiff"))
  lowrestiff(paste0(prefixPC,"_UMAP_split_by_type_",levels(seurat_sct)[i],"_lowres.tiff"))
}

# by Batch 

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "Batch",
        cols = aging_palette,
        raster = F)
hirestiff(paste0(prefixPC,"_UMAP_by_Batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_Batch_lowres.tiff"))

DimPlot(seurat_sct, reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "Batch", 
        group.by = "Batch",
        cols = aging_palette,
        raster = F) + NoLegend()
hirestiff(paste0(prefixPC,"_UMAP_split_by_Batch_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_Batch_lowres.tiff"))

Idents(seurat_sct) <- "Batch"

for (i in 1:length(levels(seurat_sct))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .25, 
                cols = aging_palette[i],
                raster = F,
                cells = c(WhichCells(seurat_sct, idents = levels(seurat_sct)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct)[i])))
  
  print(p)
  hirestiff(paste0(prefixPC,"_UMAP_split_by_Batch_",levels(seurat_sct)[i],"_hires.tiff"))
  lowrestiff(paste0(prefixPC,"_UMAP_split_by_Batch_",levels(seurat_sct)[i],"_lowres.tiff"))
}

rm(p)
  #### Clustering and Resolution ####

DefaultAssay(seurat_sct) <- "SCT"

# Determine the K-nearest neighbor graph
seurat_sct <- FindNeighbors(object = seurat_sct, 
                            reduction = "pca", 
                            dims = 1:meaningful.PCs)

# Determine the clusters                              

seurat_sct <- FindClusters(object = seurat_sct,
                           resolution = c(0.6))

res <- "_res.0.6"
ident_res <- paste0("SCT_snn",res)

Idents(seurat_sct) <- ident_res
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        raster = F) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

#update prefix
prefixPCres <- paste0(prefixPC,res)

# Plot the UMAP
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        raster = F)  + NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","by","cluster","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","by","cluster","lowres.tiff",sep = "_"))

# UMAP of cells in each cluster by type without cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        split.by = "Type",
        raster = F)  + NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","split_by","type","without_cluster_labels","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","split_by","type","without_cluster_labels","lowres.tiff",sep = "_"))

for (i in 1:length(levels(seurat_sct$Type))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .4,
                raster = F, 
                cells = c(WhichCells(seurat_sct, expression = Type == levels(seurat_sct$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","split_by","type","no_cluster_labels",
                  levels(seurat_sct$Type)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","split_by","type","no_cluster_labels",
                   levels(seurat_sct$Type)[i],"lowres.tiff",sep = "_"))
}

# UMAP of cells in each cluster by type with cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE,
        label.size = 5, 
        split.by = "Type",
        raster = F)  + NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","split_by","type","with_cluster_labels","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","split_by","type","with_cluster_labels","lowres.tiff",sep = "_"))

for (i in 1:length(levels(seurat_sct$Type))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = TRUE, 
                pt.size = .4,
                raster = F,
                label.size = 5, 
                cells = c(WhichCells(seurat_sct, expression = Type == levels(seurat_sct$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","split_by","type","with_cluster_labels",
                  levels(seurat_sct$Type)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","split_by","type","with_cluster_labels",
                   levels(seurat_sct$Type)[i],"lowres.tiff",sep = "_"))
}

# UMAP of cells in each cluster by batch with cluster labels
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = F, 
        split.by = "Batch",
        raster = F,
        label.size = 5) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","split_by","batch","without_cluster_labels","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","split_by","batch","without_cluster_labels","lowres.tiff",sep = "_"))

rm(p)

  #### Normalize RNA, find variable features, scale data for visualization ####

# Select the RNA counts slot to be the default assay for visualization purposes

DefaultAssay(seurat_sct) <- "RNA"

# Normalize, find variable features, scale data 
seurat_sct <- NormalizeData(seurat_sct, verbose = FALSE)
seurat_sct <- FindVariableFeatures(seurat_sct)
all.genes <- rownames(seurat_sct)
seurat_sct <- ScaleData(seurat_sct, features = all.genes)

  #### Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(seurat_sct, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
n_cells
write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_cluster.csv"))

# number of cells per Type
n_cells <- FetchData(seurat_sct, vars = c("Type")) %>%
  dplyr::count(Type) 
n_cells
write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_type.csv"))


# number of cells per Type per cluster
n_cells <- FetchData(seurat_sct, vars = c("ident", "Type")) %>%
  dplyr::count(ident, Type) %>%
  tidyr::spread(ident, n)

n_cells
write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_type_per_cluster.csv"))

rm(n_cells)

# save RDS containing reduction and cluster idents
# saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_clustering.rds"))

# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_clustering.rds"))

  #### save RDS containing RNA normalized data ####
# saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

  

#### LABEL CELL TYPES ####
  #### Determine cell types ####

# 0 - Excitatory Neurons EXCN
VlnPlot(seurat_sct,
        features = c("Slc17a7","Neurod6","Camk2a","Nell2"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 0,2,3,6 - Excitatory Neurons")
hirestiff(paste(prefixPCres,"vln","Cluster","0_2_3_6","ExNeurons","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","0_2_3_6","ExNeurons","lowres.tiff", sep = "_"))

# 1 - DG Granule (which also express glut markers) GRN
VlnPlot(seurat_sct,
        features = c("Bdnf","Calb1","Calm1","Calm2","Calm3",
                     "Ppp3ca","Disc1","Gabrb1","Gabra4","Gria1",
                     "Gria2","Gria3","Grin1","Grin2b","Chrm1",
                     "Nrgn","Nrp2","Pcp4","Ncam1","Prox1",
                     "Vsnl1","Npy1r","Actn2","Grm1","Grm2",
                     "Grm3","Grm5","Creb1","Slc17a7"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 1,23 - DG Granule Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","1_23","DG_Granule","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","1_23","DG_Granule","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Bdnf","Calb1","Disc1","Gabrb1","Gabra4",
                     "Gria1","Gria2","Gria3","Grin1","Grin2b",
                     "Chrm1","Nrp2","Ncam1","Prox1","Npy1r",
                     "Actn2","Grm1","Grm2","Grm3","Grm5",
                     "Creb1","Slc17a7"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 1,23 - DG Granule Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","1_23","DG_Granule","refined_list","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","1_23","DG_Granule","refined_list","lowres.tiff", sep = "_"))

# 2 - EXCN

# 3 - EXCN

# 4 - Oligodendrocytes OLIG
VlnPlot(seurat_sct,
        features = c("Cnp","Cldn11","Mbp","Mog","Olig1","Sox10"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 4, 10, 20 - Oligodendrocytes")
hirestiff(paste(prefixPCres,"vln","Cluster","4_10_20","Oligodendrocytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","4_10_20","Oligodendrocytes","lowres.tiff", sep = "_"))

# 5 - PYR Pyramidal PYR
VlnPlot(seurat_sct,
        features = c("Amigo2","Disc1","Bdnf","Bok",
                     "Ppp3ca","Gabrb1","Gabra4","Gabra5",
                     "Gria1","Man1a","Neurod6","Chrm1","Chrm3",
                     "Nrp2","Pou3f1","Npy2r",
                     "Actn2","Grm1","Grm2","Grm3","Grm4",
                     "Grm5","Pcp4","Creb1"),
        stack = TRUE,
        flip = TRUE) +  
  NoLegend() + ggtitle("Clusters 5, 8, 9 - Pyramidal Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","5_8_9","pyramidal","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","5_8_9","pyramidal","lowres.tiff", sep = "_"))

# CA2 cannonical
VlnPlot(seurat_sct,
        features = c(
          # "Accn1", #not found
                     "Adcy5",
                     "Cacng5",
                     "Ccdc3",
                     # "Fam40b", #not found
                     "Gpr12",
                     "Gsto1",
                     "Map3k15",
                     "Ptpn5",
                     "Prss23",
                     "Rgs14",
                     "S100b",
                     "Sostdc1"),
        stack = TRUE,
        flip = TRUE,
        idents = c(5,8,9)) + 
  NoLegend()

# CA2 "new"
VlnPlot(seurat_sct,
        features = c(
          # "170024P16Rik", # not found
                     "Aldh1a1",
                     "Arg2",
                     "Crabp1",
                     "Ctsc",
                     "Dusp5",
                     "F2r",
                     "Fgf2",
                     "Fgf5",
                     "Glul",
                     # "Gm13847", # not found
                     # "Gm13921", # not found
                     # "Gm20459", # not found
                     "Maob",
                     "Ntsr2",
                     "Plcb4",
                     "Pygo1",
                     "Rgs5",
                     "Scgn",
                     "Srgap2",
                     "Srl",
                     "Syce2",
                     "Tgfb1i1",
                     "Vcan",
                     "Vit"
        ),
        stack = TRUE,
        flip = TRUE,
        idents = c(5,8,9)) +  
  NoLegend()

# Broad Hippocampal populations (Cembrowski)
VlnPlot(seurat_sct,
        features = c("Prox1", # GRN
                     "Dkk3", # Non-GRN
                     "Calb2", # Mossy
                     "Ociad2", # All Pyrmidal
                     "Cacng5", # CA2 Pyr
                     "Fibcd1", # CA1 Pyr
                     "Pdzd2", # GRN dorsal
                     "Tox3", # GRN ventral
                     "Iyd", # CA3 dorsal
                     "Coch", # CA3 ventral
                     "Wfs1", # CA1 dorsal
                     "Dcn"), # CA1 ventral
        # idents = c(5,8,9),
        sort = TRUE,
        stack = TRUE,
        flip = TRUE
        ) +  
  NoLegend()

# 6 - EXCN

# 7 - Inhibitory Neurons INHN
VlnPlot(seurat_sct,
        features = c("Gad1","Synpr","Gad2","Sst"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 7 - Inhibitory Neurons")
hirestiff(paste(prefixPCres,"vln","Cluster","7","Inhibitory","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","7","Inhibitory","lowres.tiff", sep = "_"))

# 8 - EXCN

# 9 - PYR

# 10 - OLIG

# 11 - Astrocytes ASTR
VlnPlot(seurat_sct,
        features = c("Rgs20","Aqp4","Ntsr2","Gfap","S100b"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 11 - Astrocytes")
hirestiff(paste(prefixPCres,"vln","Cluster","11","Astrocytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","11","Astrocytes","lowres.tiff", sep = "_"))

# 12 - DG Basket BSKT
VlnPlot(seurat_sct,
        features = c("Dlx6os1","Maf","Erbb4","Gad2","Grik1",
                     "Nxph1","Slc6a1","Gad1","Gm13629","Dlx1",
                     "Dlx6","Col19a1","Utrn","Reln","Npas1",
                     "Sox6","Crhbp","Npy","Sst",
                     "Gabra1","Pvalb"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 12 - DG Basket Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","12","DG","Basket","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","12","DG","Basket","lowres.tiff", sep = "_"))

# 13 - Microglia MG
VlnPlot(seurat_sct,
        features = c("Itgam","Ptprc","Tmem119","Cx3cr1","P2ry12","C1qa"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 13 - Microglia")
hirestiff(paste(prefixPCres,"vln","Cluster","13","Microglia","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","13","Microglia","lowres.tiff", sep = "_"))

# 14 - OPC
VlnPlot(seurat_sct,
        features = c("Pdgfra","Cspg4","Sox10","Olig1","Olig2"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 14 - OPCs")
hirestiff(paste(prefixPCres,"vln","Cluster","14","OPCs","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","14","OPCs","lowres.tiff", sep = "_"))

# 15 - DG Mossy MOSS
VlnPlot(seurat_sct,
        features = c("Gad1","Gad2","Calb2",
                     "Gria2","Gria3","Pcp4",
                     "Gabra6","Pvalb","Nos1"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 15 - DG Mossy Cells (neg for Gabra6, Pvalb & Nos1)")
hirestiff(paste(prefixPCres,"vln","Cluster","15","DG","Mossy","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","15","DG","Mossy","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Gria2","Gria3","Pcp4"),
        stack = TRUE,
        flip = TRUE,
        idents = c("15"),
        split.by = "Type",
        cols = aging_palette) + 
  ggtitle("Cluster 15 - DG Mossy Cells split by Type")
lowrestiff(paste(prefixPCres,"vln","Cluster","15","DG","Mossy","by","type","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Gria2"),
        split.by = "Type",
        idents = c("15"),
        flip = TRUE,
        cols = aging_palette,
        pt.size = 1)
lowrestiff(paste(prefixPCres,"vln","Cluster","15","DG","Mossy","by","type","Gria2","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Gria3"),
        split.by = "Type",
        idents = c("15"),
        flip = TRUE,
        cols = aging_palette,
        pt.size = 1)
lowrestiff(paste(prefixPCres,"vln","Cluster","15","DG","Mossy","by","type","Gria3","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Pcp4"),
        split.by = "Type",
        idents = c("15"),
        flip = TRUE,
        cols = aging_palette,
        pt.size = 1)
lowrestiff(paste(prefixPCres,"vln","Cluster","15","DG","Mossy","by","type","Pcp4","lowres.tiff", sep = "_"))

# 16 - CPCa

# CPC Choroid Plexus Cells
VlnPlot(seurat_sct,features = c("Ttr","Atp5g1","Clic6","Ndufv3","Uqcrb",
                                "Hemk1","Atp5l","Ndufa4","Ppp1r1b","Cox7b",
                                "Atp5mpl","Atp5md","Ndufa1","Uqcrq","Cox6c",
                                "Cox7c","Atp5j","Atp5e","Uqcr10","Cox8b",
                                "Acad8","Uqcrh","Atp5j2","Uqcr11","Cox5b",
                                "Cox8a","Cox6b1","Ldhb"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 16 - Choroid Plexus Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","16","CPC","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","16","CPC","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Ttr","Clic6","Hemk1","Ppp1r1b",
                     "Acad8","Folr1","Prlr","Kl"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 16 - Choroid Plexus Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","16","CPC","refined_list","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","16","CPC","refined_list","lowres.tiff", sep = "_"))

# 17 - Vascular Leptomeningeal Cells VLMC
VlnPlot(seurat_sct,
        features = c("Pdgfra","Gja1","Slc6a13"),
        stack = TRUE,flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 17 - VLMC")
hirestiff(paste(prefixPCres,"vln","Cluster","17","VLMC","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","17","VLMC","lowres.tiff", sep = "_"))

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Pdgfra","Gja1","Slc6a13","Nes"), 
            order = TRUE,
            # min.cutoff = 'q10', 
            label = TRUE,
            label.size = 3,
            pt.size = 0.3,
            raster = F)

# 18 - CPCb

# 19 - Pericytes PERI
VlnPlot(seurat_sct,
        features = c("Kcnj8","Cspg4","Pdgfrb","Rgs5"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 20 - Pericytes")
hirestiff(paste(prefixPCres,"vln","Cluster","20","Pericytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","20","Pericytes","lowres.tiff", sep = "_"))

# 20 - OLIG

# 21 - A1 Astrocytes ASTR-A1
VlnPlot(seurat_sct,
        features = c("Rgs20","Aqp4","Ntsr2","Gfap","S100b","C3"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 21 - A1 Astrocytes")

hirestiff(paste(prefixPCres,"vln","Cluster","21","A1","Astrocytes","res.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","21","A1","Astrocytes","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("C3"),
        split.by = "Type",
        idents = c("11","21"),
        flip = TRUE,
        cols = aging_palette,
        pt.size = 1)

hirestiff(paste(prefixPCres,"vln","astro","A1","C3","by","Type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","astro","A1","C3","by","Type","lowres.tiff", sep = "_"))

# 22 - Activated Microglia MG-A
VlnPlot(seurat_sct,features = c("Itgam","Tmem119","Cx3cr1",
                                "P2ry12","Mrc1",
                                "Jun","Junb","Jund","C3"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 22 - Activated Microglia")
hirestiff(paste(prefixPCres,"vln","Cluster","22","Activated","Microglia","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","22","Activated","Microglia","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,features = c("Itgam","Tmem119","Cx3cr1",
                                "Tlr4","Fosb","Jun","Junb","Jund"),
        idents = c(13,22),
        split.by = "Type",
        stack = TRUE,
        flip = TRUE) + 
  ggtitle("Cluster 21 - Activated Microglia")
hirestiff(paste(prefixPCres,"vln","Cluster","13_22","Microglia","by_type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","13_22","Microglia","by_type","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Mrc1"),
        split.by = "Type",
        idents = c("13","22"),
        flip = TRUE,
        pt.size = 1)
lowrestiff(paste(prefixPCres,"vln","Cluster","22","Activated","Microglia","Mrc1","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Tlr4"),
        split.by = "Type",
        idents = c("13","22"),
        flip = TRUE,
        pt.size = 1)
lowrestiff(paste(prefixPCres,"vln","Cluster","22","Activated","Microglia","Tlr4","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Cd163"),
        split.by = "Type",
        idents = c("13","22"),
        flip = TRUE,
        pt.size = 1)
lowrestiff(paste(prefixPCres,"vln","Cluster","22","Activated","Microglia","Cd163","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,features = c("Tmem119"),split.by = "Type",idents = c("13","22"),flip = TRUE)
VlnPlot(seurat_sct,features = c("Jund"),split.by = "Type",idents = c("13","22"),flip = TRUE)
VlnPlot(seurat_sct,features = c("Mapk14"),split.by = "Type",idents = c("13","22"),flip = TRUE)

VlnPlot(seurat_sct,
        features = c("C3"),
        split.by = "Type",
        idents = c("13","22"),
        flip = TRUE,
        cols = aging_palette,
        pt.size = 1)

hirestiff(paste(prefixPCres,"vln","MG","C3","by","Type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","MG","C3","by","Type","lowres.tiff", sep = "_"))

# 23 - GRN

  #### Label clusters as cell types ####

# assign cell types to cluster numbers
seurat_sct <- RenameIdents(object = seurat_sct, 
                           '0' = "EXCN",
                           '1' = "GRN", 
                           '2' = "EXCN", 
                           '3' = "EXCN",
                           '4' = "OLIG",
                           '5' = "PYR",
                           '6' = "EXCN",
                           '7' = "INHN",
                           '8' = "PYR",
                           '9' = "PYR",
                           '10' = "OLIG",
                           '11' = "ASTR",
                           '12' = "BSKT",
                           '13' = "MG",
                           '14' = "OPC",
                           '15' = "MOSS",
                           '16' = "CPCa",
                           '17' = "VLMC",
                           '18' = "CPCb",
                           '19' = "PERI",
                           '20' = "OLIG",
                           '21' = "ASTR-A1",
                           '22' = "MG",
                           '23' = "GRN")

# set the order of cell types
seurat_sct@active.ident <- factor(x = seurat_sct@active.ident, 
                                  levels = c("INHN","BSKT","MOSS",
                                             "EXCN",
                                             "PYR",
                                             "GRN",
                                             "MG",
                                             "ASTR","ASTR-A1",
                                             "OLIG","OPC",
                                             "VLMC","CPCa","CPCb","PERI"))

# set cell type color
cell_cols <- c("INHN"="#00BFC4","BSKT"="#B79F00","MOSS"="#00BA38",
               "EXCN"="#F8766D",
               "PYR"="#619CFF",
               "GRN"="#E88526",
               "MG"="#93AA00",
               "ASTR"="#00B9E3","ASTR-A1"="#F17D50",
               "OLIG"="#8E92FF","OPC"="#D39200",
               "VLMC"="#5EB300","CPCa"="#00C19F","CPCb"="#BE92FF","PERI"="#00ADFA")

  #### save cell type as a metadata column ####
seurat_sct[["Cell.Type"]] <- seurat_sct@active.ident

  #### switch between cluster numbers and cell type idents ####

Idents(seurat_sct) <- ident_res

Idents(seurat_sct) <- "Cell.Type"

  #### UMAPS with cell type labels ####
DimPlot(seurat_sct,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        pt.size = 0.5,
        cols = cell_cols,
        raster = F) + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","cell","types","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","cell","types","lowres.tiff", sep = "_"))
hirestiffsquare(paste(prefixPCres,"UMAP","cell","types","hires_square.tiff", sep = "_"))
savepdfsquare(paste(prefixPCres,"UMAP","cell","types","square.pdf", sep = "_"))

DimPlot(seurat_sct,
        reduction = "umap",
        label = F,
        label.size = 6,
        pt.size = 0.5,
        cols = cell_cols,
        raster = F,
        split.by = "Batch") + 
  NoLegend()

DimPlot(seurat_sct,
        reduction = "umap",
        label = F,
        label.size = 6,
        pt.size = 0.5,
        cols = cell_cols,
        raster = F,
        split.by = "Type") + 
  NoLegend()

DimPlot(seurat_sct,
        reduction = "umap",
        label = F,
        label.size = 6,
        pt.size = 0.5,
        cols = cell_cols,
        raster = F,
        split.by = "orig.ident") + 
  NoLegend()

DimPlot(seurat_sct,
        reduction = "umap",
        label = FALSE,
        label.size = 4,
        pt.size = 0.5,
        cols = cell_cols,
        raster = F) + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","cell","types","no","label","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","cell","types","no","label","lowres.tiff", sep = "_"))
hirestiffsquare(paste(prefixPCres,"UMAP","cell","types","no","label","hires_square.tiff", sep = "_"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 3, 
        pt.size = 0.3,
        split.by = "Type",
        cols = cell_cols,
        raster = F) + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","cell","types","by","type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","cell","types","by","type","lowres.tiff", sep = "_"))

for (i in 1:length(levels(seurat_sct$Type))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = TRUE,
                label.size = 6,
                pt.size = 0.5,
                cols = cell_cols,
                raster = F, 
                cells = c(WhichCells(seurat_sct, 
                                     expression = Type == levels(seurat_sct$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                  levels(seurat_sct$Type)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                   levels(seurat_sct$Type)[i],"lowres.tiff",sep = "_"))
  hirestiffsquare(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                        levels(seurat_sct$Type)[i],"hires_square.tiff",sep = "_"))
}

for (i in 1:length(levels(seurat_sct$Type))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = FALSE,
                label.size = 4,
                pt.size = 0.5,
                cols = cell_cols, 
                cells = c(WhichCells(seurat_sct, 
                                     expression = Type == levels(seurat_sct$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                  levels(seurat_sct$Type)[i],"no","label","hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                   levels(seurat_sct$Type)[i],"no","label","lowres.tiff",sep = "_"))
  hirestiffsquare(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                        levels(seurat_sct$Type)[i],"no","label","hires_square.tiff",sep = "_"))
}

for (i in 1:length(levels(seurat_sct$Type))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = FALSE,
                label.size = 4,
                pt.size = 0.5,
                cols = cell_cols, 
                cells = c(WhichCells(seurat_sct, 
                                     expression = Type == levels(seurat_sct$Type)[i]))) + 
          ggtitle(paste0(levels(seurat_sct$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                  levels(seurat_sct$Type)[i],"no","label","with","legend","hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                   levels(seurat_sct$Type)[i],"no","label","with","legend","lowres.tiff",sep = "_"))
  hirestiffsquare(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                        levels(seurat_sct$Type)[i],"no","label","with","legend","hires_square.tiff",sep = "_"))
}

rm(p,i)
  #### All cell types VLN ####

VlnPlot(seurat_sct,features = c(
                                "Gad1","Gad2",            ## INHN
                                "Dlx1","Col19a1",         ## BSKT
                                "Cnr1","Calb2",           ## MOSS
                                "Slc17a7","Neurod6",      ## EXCN
                                "Grm1","Satb2",           ## PYR
                                "Calb1","Prox1",          ## GRN
                                "Tmem119","Itgam",        ## MG
                                "Aqp4","Gfap",            ## ASTR
                                "S100b","C3",             ## ASTR-A1
                                "Mbp","Mog",              ## OLIG
                                "Olig1","Sox10",          ## OPC
                                "Pdgfra","Slc6a13",       ## VLMC
                                "Hemk1","Acad8",          ## CPCa
                                "Folr1","Kl",             ## CPCb
                                "Kcnj8","Rgs5"            ## PERI
                                             ),
        stack = TRUE,
        flip = TRUE,
        sort = FALSE) + 
  NoLegend() + 
  ggtitle("All Cell Types") &
  theme(
        axis.text.x = element_text(size = 20), 
        axis.title = element_text(size = 22),
        axis.text.y.left = element_text(size = 12),
        axis.text.y.right = element_text(size = 12),
        strip.text = element_text(size = 15,
                                  face = "bold.italic"))

hirestiff(paste(prefixPCres,"VLN","cell","type","markers","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"VLN","cell","type","markers","lowres.tiff", sep = "_"))
hirestiffsquare(paste(prefixPCres,"VLN","cell","type","markers","hires_square.tiff", sep = "_"))
savepdfsquare(paste(prefixPCres,"VLN","cell","type","markers","square.pdf", sep = "_"))

  #### Extract number of cells per cluster per Cell.Type ####
n_cells <- FetchData(seurat_sct, vars = c("orig.ident", "Cell.Type")) %>%
  dplyr::count(orig.ident, Cell.Type) %>%
  tidyr::spread(orig.ident, n)
n_cells
write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_cell.type.csv"))

n_cells <- FetchData(seurat_sct, vars = c("Type", "Cell.Type")) %>%
  dplyr::count(Type, Cell.Type) %>%
  tidyr::spread(Type, n)
n_cells
write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_type_per_cell.type.csv"))

rm(n_cells)

# number of cells per Type
FetchData(seurat_sct, vars = c("Type")) %>%
  dplyr::count(Type) 
  #### DEGs within cell type clusters ####
# for loop for DEGs within cell type clusters

# get the names of cell types
clusterlist <- unique(seurat_sct@meta.data$Cell.Type)

# list the groups to compare.  Since we're interested in 
# Young Veh vs Aged Veh, 
# Aged Veh vs Aged iMAC, 
# and Aged iMAC vs Young Veh,
# I added an additional Young Veh to the end of the list.
# The loop will compare the first and second terms,
# then the second and third terms,
# then the third and fourth terms and then will terminate
grouplist <- c("Young Veh","Aging Veh","Aging iMPs","Young Veh")

for (i in 1:(length(grouplist)-1)){
  for (k in 1:length(clusterlist)){
    
    nam <- paste("clusterDEG", grouplist[i],"vs",grouplist[i+1],"cluster",clusterlist[k], sep = "_")
    assign(nam, FindMarkers(seurat_sct, 
                            ident.1 = grouplist[i], 
                            ident.2 = grouplist[i+1],
                            group.by = "Type", 
                            subset.ident = clusterlist[k]))
    
    write.csv(get(nam), file = paste(prefixPCres, nam, "list.csv", sep = "_"), row.names = TRUE)
    
  } 
}

  #### Save cell type labeled seurat obj ####
saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_cell_types_labeled_mg_combined.rds"))

# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_cell_types_labeled_mg_combined.rds"))

#### AGING CLOCKS ####
  #### Bootstrap Cells ####

# Purpose:
# Take a seurat object and return a bootstrapcell dataframe
# Technique of cell type specific pseudo bulking to optimize
# trade off between cell count and cell complexity
# Bootstrap sampling rather than random partitions
# Equally weights samples rather than cells

# convert data to dataframe
data <- seurat_sct

Convert_to_Dataframe <- function(svz) {
  DefaultAssay(svz) <- "RNA"
  meta <- svz@meta.data
  meta <- meta[, c("orig.ident", "Type", "Cell.Type")]
  raw_counts <- t(as.matrix(svz[["RNA"]]@counts))
  raw_counts <- raw_counts[, colSums(raw_counts) > 0]
  df <- as_tibble(cbind(meta, raw_counts))
  return(df)
}

df <- Convert_to_Dataframe(data) %>%
  group_by(Cell.Type, Type, orig.ident) %>%
  nest()

head(df)
# A tibble: 6 × 4
# Groups:   Cell.Type, Type, orig.ident [6]
# orig.ident Type      Cell.Type pseudocell_all         
# <fct>      <fct>     <fct>     <list>                 
#   1 YV_106     Young Veh EXCN      <tibble [100 × 28,902]>
#   2 YV_106     Young Veh PYR       <tibble [100 × 28,902]>
#   3 YV_106     Young Veh GRN       <tibble [100 × 28,902]>
#   4 YV_106     Young Veh BSKT      <tibble [100 × 28,902]>
#   5 YV_106     Young Veh MOSS      <tibble [100 × 28,902]>
#   6 YV_106     Young Veh MG        <tibble [100 × 28,902]> 

rm(data)

# bootstrap pseudocells
bootstrap.pseudocells <- function(df, size=15, n=100, replace="dynamic") {
  pseudocells <- c()
  # If dynamic then only sample with replacement if required due to shortage of cells.
  if (replace == "dynamic") {
    if (nrow(df) <= size) {replace <- TRUE} else {replace <- FALSE}
  }
  for (i in c(1:n)) {
    batch <- df[sample(1:nrow(df), size = size, replace = replace), ]
    pseudocells <- rbind(pseudocells, colSums(batch))
  }
  colnames(pseudocells) <- colnames(df)
  return(as_tibble(pseudocells))
}

# Apply boostrap.pseudocells using map()
set.seed(42)

df2 <- df %>% mutate(pseudocell_all = map(data, bootstrap.pseudocells)) # ~30 minutes

# Remove single cell data; keep just key metadata and pseudocells
df2$data <- NULL

head(df2)

saveRDS(df2, paste(prefixPCres,"bootstrap_dataframe_all_samples.rds", sep = "_"))

  #### Age in Months ####
    #### Train with 4/5 YV and 4/5 AV samples ####

df1 <- readRDS(paste(prefixPCres,"bootstrap_dataframe_all_samples.rds", sep = "_"))

# Add numeric Months column dependent on orig.ident

# c( "YV_102","YV_104","YV_106","YV_111","YV_116",
#    "AV_011","AV_023","AV_026","AV_029","AV_033",
#    "AI_017","AI_025","AI_031","AI_045","AI_046")

# Ages in Days / Months
#
# 102	86 3
# 104	86 3
# 106	86 3
# 111	84 3
# 116	57 2
#
# 11 457 15
# 23 413 14
# 26 413 14
# 29 413 14
# 33 413 14
#
# 17 457 15
# 25 413 14
# 31 413 14
# 45 405 14
# 46 405 14

df1$Months <- ifelse(df1$orig.ident == "YV_102", 3,
  ifelse(df1$orig.ident == "YV_104", 3,
    ifelse(df1$orig.ident == "YV_106", 3,
      ifelse(df1$orig.ident == "YV_111", 3,
        ifelse(df1$orig.ident == "YV_116", 2,
          ifelse(df1$orig.ident == "AV_011", 15,
            ifelse(df1$orig.ident == "AV_023", 14,
              ifelse(df1$orig.ident == "AV_026", 14,
                ifelse(df1$orig.ident == "AV_029", 14,
                  ifelse(df1$orig.ident == "AV_033", 14,
                    ifelse(df1$orig.ident == "AI_017", 15,
                      ifelse(df1$orig.ident == "AI_025", 14,
                        ifelse(df1$orig.ident == "AI_031", 14,
                          ifelse(df1$orig.ident == "AI_045", 14,
                            ifelse(df1$orig.ident == "AI_046", 14, NA)
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)

# reorder columns so that Months is before the pseudocell data
column_order <- c("orig.ident","Type","Cell.Type","Months", 
                  setdiff(names(df1),
                          c("orig.ident","Type","Cell.Type","Months")))

df1 <- df1 %>%
  select(all_of(column_order))

# # subset out all Aging iMP samples for later testing
# df_AI <- df1 %>% filter(Type == "Aging iMPs")
# 
# df_AI <- droplevels(df_AI)

# randomly choose which samples to hold out on each round of cross validation
GroupYV <- c("YV_102","YV_104","YV_106","YV_111","YV_116")
GroupAV <- c("AV_011","AV_029","AV_033","AV_023","AV_026")

# Initialize a list to store the samples
sample_results <- list()

# Sample from the two groups 5 times without replacement
for (i in 1:5) {
  sample_YV <- sample(GroupYV, 1)
  sample_AV <- sample(GroupAV, 1)
  
  # Create a list to store samples from this round
  round_samples <- list(Sample_YV = sample_YV, Sample_AV = sample_AV)
  
  # Append the list to the sample_results list
  sample_results[[i]] <- round_samples
  
  # Remove the chosen elements from the group vectors
  GroupYV <- setdiff(GroupYV, sample_YV)
  GroupAV <- setdiff(GroupAV, sample_AV)
}

# Access the results from each round
for (i in 1:5) {
  cat("Round", i, "\n")
  cat("Sample from YV:", sample_results[[i]]$Sample_YV, "\n")
  cat("Sample from AV:", sample_results[[i]]$Sample_AV, "\n\n")
}

# Round 1 
# Sample from YV: YV_104 
# Sample from AV: AV_029 
# 
# Round 2 
# Sample from YV: YV_102 
# Sample from AV: AV_033 
# 
# Round 3 
# Sample from YV: YV_116 
# Sample from AV: AV_023 
# 
# Round 4 
# Sample from YV: YV_106 
# Sample from AV: AV_011 
# 
# Round 5 
# Sample from YV: YV_111 
# Sample from AV: AV_026 

# subset out one young and one aging veh for validation testing
# make the five sets of hold outs (one YV, one AV)

for (i in 1:5) {
  nam <- paste("df_test",i, sep = "_")
  assign(nam, df1 %>% 
           filter(orig.ident %in% c(sample_results[[i]]$Sample_YV,
                                    sample_results[[i]]$Sample_AV)))
  
} 

df_test_1 <- droplevels(df_test_1)
df_test_2 <- droplevels(df_test_2)
df_test_3 <- droplevels(df_test_3)
df_test_4 <- droplevels(df_test_4)
df_test_5 <- droplevels(df_test_5)

# subset out all remaining samples
df_train_1 <- df1 %>% filter(!orig.ident %in% c(sample_results[[1]]$Sample_YV,
                                                sample_results[[1]]$Sample_AV))
df_train_2 <- df1 %>% filter(!orig.ident %in% c(sample_results[[2]]$Sample_YV,
                                                sample_results[[2]]$Sample_AV))
df_train_3 <- df1 %>% filter(!orig.ident %in% c(sample_results[[3]]$Sample_YV,
                                                sample_results[[3]]$Sample_AV))
df_train_4 <- df1 %>% filter(!orig.ident %in% c(sample_results[[4]]$Sample_YV,
                                                sample_results[[4]]$Sample_AV))
df_train_5 <- df1 %>% filter(!orig.ident %in% c(sample_results[[5]]$Sample_YV,
                                                sample_results[[5]]$Sample_AV))

df_train_1 <- droplevels(df_train_1)
df_train_2 <- droplevels(df_train_2)
df_train_3 <- droplevels(df_train_3)
df_train_4 <- droplevels(df_train_4)
df_train_5 <- droplevels(df_train_5)

# drop the orig.ident column since we're calculating by celltype within Type
df_train_1 <- df_train_1 %>% ungroup %>%
  select(-c(orig.ident))
df_train_2 <- df_train_2 %>% ungroup %>%
  select(-c(orig.ident))
df_train_3 <- df_train_3 %>% ungroup %>%
  select(-c(orig.ident))
df_train_4 <- df_train_4 %>% ungroup %>%
  select(-c(orig.ident))
df_train_5 <- df_train_5 %>% ungroup %>%
  select(-c(orig.ident))

# function for log normalizing
lognorm <- function(input) {
  norm <- sweep(input, MARGIN = 1, FUN = "/", STATS = rowSums(input))
  log1p(norm * 10000)
}

# log normalize, drop pseudocell data, unnest the tibble, renest by celltype
df_train_LN_1 <- df_train_1 %>% mutate(lognorm = map(pseudocell_all, lognorm))
df_train_LN_1$pseudocell_all <- NULL
df_train_LN_1 <- unnest(df_train_LN_1, lognorm)
df_train_LN_1 <- df_train_LN_1 %>% select(-Type)
by_celltype_1 <- df_train_LN_1 %>% group_by(Cell.Type) %>% nest()
colnames(by_celltype_1)[2] <- "lognormalized"
rm(df_train_1,df_train_LN_1)

df_train_LN_2 <- df_train_2 %>% mutate(lognorm = map(pseudocell_all, lognorm))
df_train_LN_2$pseudocell_all <- NULL
df_train_LN_2 <- unnest(df_train_LN_2, lognorm)
df_train_LN_2 <- df_train_LN_2 %>% select(-Type)
by_celltype_2 <- df_train_LN_2 %>% group_by(Cell.Type) %>% nest()
colnames(by_celltype_2)[2] <- "lognormalized"
rm(df_train_2,df_train_LN_2)

df_train_LN_3 <- df_train_3 %>% mutate(lognorm = map(pseudocell_all, lognorm))
df_train_LN_3$pseudocell_all <- NULL
df_train_LN_3 <- unnest(df_train_LN_3, lognorm)
df_train_LN_3 <- df_train_LN_3 %>% select(-Type)
by_celltype_3 <- df_train_LN_3 %>% group_by(Cell.Type) %>% nest()
colnames(by_celltype_3)[2] <- "lognormalized"
rm(df_train_3,df_train_LN_3)

df_train_LN_4 <- df_train_4 %>% mutate(lognorm = map(pseudocell_all, lognorm))
df_train_LN_4$pseudocell_all <- NULL
df_train_LN_4 <- unnest(df_train_LN_4, lognorm)
df_train_LN_4 <- df_train_LN_4 %>% select(-Type)
by_celltype_4 <- df_train_LN_4 %>% group_by(Cell.Type) %>% nest()
colnames(by_celltype_4)[2] <- "lognormalized"
rm(df_train_4,df_train_LN_4)

df_train_LN_5 <- df_train_5 %>% mutate(lognorm = map(pseudocell_all, lognorm))
df_train_LN_5$pseudocell_all <- NULL
df_train_LN_5 <- unnest(df_train_LN_5, lognorm)
df_train_LN_5 <- df_train_LN_5 %>% select(-Type)
by_celltype_5 <- df_train_LN_5 %>% group_by(Cell.Type) %>% nest()
colnames(by_celltype_5)[2] <- "lognormalized"
rm(df_train_5,df_train_LN_5)

# Model function
celltype_model <- function(input) {
  cv.glmnet(x = as.matrix(input[, -1]), 
            y = as.matrix(input[, 1]),
            type.measure = "mae", 
            standardize = F, 
            relax = F, 
            nfolds = 5)
}

# Apply model function using map()
models_1 <- by_celltype_1 %>% 
  mutate(model = map(lognormalized, celltype_model))
models_1$lognormalized <- NULL

models_2 <- by_celltype_2 %>% 
  mutate(model = map(lognormalized, celltype_model))
models_2$lognormalized <- NULL

models_3 <- by_celltype_3 %>% 
  mutate(model = map(lognormalized, celltype_model))
models_3$lognormalized <- NULL

models_4 <- by_celltype_4 %>% 
  mutate(model = map(lognormalized, celltype_model))
models_4$lognormalized <- NULL

models_5 <- by_celltype_5 %>% 
  mutate(model = map(lognormalized, celltype_model))
models_5$lognormalized <- NULL

# Save models
saveRDS(models_1, paste(prefixPCres,"models_1_by_month.rds", sep = "_"))
saveRDS(models_2, paste(prefixPCres,"models_2_by_month.rds", sep = "_"))
saveRDS(models_3, paste(prefixPCres,"models_3_by_month.rds", sep = "_"))
saveRDS(models_4, paste(prefixPCres,"models_4_by_month.rds", sep = "_"))
saveRDS(models_5, paste(prefixPCres,"models_5_by_month.rds", sep = "_"))

    #### predict held out samples ####

models_1 <- readRDS(paste(prefixPCres,"models_1_by_month.rds", sep = "_"))
models_2 <- readRDS(paste(prefixPCres,"models_2_by_month.rds", sep = "_"))
models_3 <- readRDS(paste(prefixPCres,"models_3_by_month.rds", sep = "_"))
models_4 <- readRDS(paste(prefixPCres,"models_4_by_month.rds", sep = "_"))
models_5 <- readRDS(paste(prefixPCres,"models_5_by_month.rds", sep = "_"))

    #### Round 1 ####
# Sample from YV: YV_104 
# Sample from AV: AV_029 

models <- models_1
samples <- df_test_1

samples <- samples %>% unnest(pseudocell_all)

meta <- as.data.frame(samples[, c(1:4)])
umi <- samples[, -c(1:4)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)

by_celltype <- as_tibble(cbind(meta,logged)) %>%
  group_by(Cell.Type) %>%
  nest()

by_celltype <- dplyr::inner_join(by_celltype, models, by = "Cell.Type")

by_celltype[1,2][[1]][[1]][1:5,1:5]

custom_add_predictions <- function(data, model) {
  predictions <- predict(model, 
                         newx = as.matrix(data[,-(1:3)]), 
                         s = "lambda.min")
  p <- as_tibble(data.frame("Pred" = predictions[,1],
                            "orig.ident" = data[,1],
                            "Type" = data[,2],
                            "Months" = data[,3]))
  return(p)
}

# Add predictions 
by_celltype <- by_celltype %>%
  mutate(Predictions = map2(data, model, custom_add_predictions))

output <- by_celltype %>% select(-c(data, model)) %>% unnest(Predictions)

range(output$Pred) # 0.09698989 16.53956264
head(output)
# A tibble: 6 × 5
# Groups:   Cell.Type [1]
# Cell.Type  Pred orig.ident Type      Months
# <fct>     <dbl> <fct>      <fct>      <dbl>
# 1 EXCN      11.8  AV_029     Aging Veh     14
# 2 EXCN      11.8  AV_029     Aging Veh     14
# 3 EXCN      11.1  AV_029     Aging Veh     14
# 4 EXCN      12.8  AV_029     Aging Veh     14
# 5 EXCN      11.4  AV_029     Aging Veh     14
# 6 EXCN       9.94 AV_029     Aging Veh     14

saveRDS(output, paste(prefixPCres,
                      "TestGroup1",
                      "YV_104",
                      "AV_029",
                      "predictions_by_month.rds",
                      sep = "_"))

write.csv(output, paste(prefixPCres,
                        "TestGroup1",
                        "YV_104",
                        "AV_029",
                        "predictions_by_month.csv",
                        sep = "_"))

    #### Round 2 ####
# Sample from YV: YV_102 
# Sample from AV: AV_033 

models <- models_2
samples <- df_test_2

samples <- samples %>% unnest(pseudocell_all)

meta <- as.data.frame(samples[, c(1:4)])
umi <- samples[, -c(1:4)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)

by_celltype <- as_tibble(cbind(meta,logged)) %>%
  group_by(Cell.Type) %>%
  nest()

by_celltype <- dplyr::inner_join(by_celltype, models, by = "Cell.Type")

by_celltype[1,2][[1]][[1]][1:5,1:5]

custom_add_predictions <- function(data, model) {
  predictions <- predict(model, 
                         newx = as.matrix(data[,-(1:3)]), 
                         s = "lambda.min")
  p <- as_tibble(data.frame("Pred" = predictions[,1],
                            "orig.ident" = data[,1],
                            "Type" = data[,2],
                            "Months" = data[,3]))
  return(p)
}

# Add predictions 
by_celltype <- by_celltype %>%
  mutate(Predictions = map2(data, model, custom_add_predictions))

output <- by_celltype %>% select(-c(data, model)) %>% unnest(Predictions)

range(output$Pred) # 0.6293697 15.6199859
head(output)
# A tibble: 6 × 5
# Groups:   Cell.Type [1]
# Cell.Type  Pred orig.ident Type      Months
# <fct>     <dbl> <fct>      <fct>      <dbl>
# 1 PYR       11.2  AV_033     Aging Veh     14
# 2 PYR        8.72 AV_033     Aging Veh     14
# 3 PYR       11.1  AV_033     Aging Veh     14
# 4 PYR       11.1  AV_033     Aging Veh     14
# 5 PYR        9.28 AV_033     Aging Veh     14
# 6 PYR       11.5  AV_033     Aging Veh     14

saveRDS(output, paste(prefixPCres,
                      "TestGroup2",
                      "YV_102",
                      "AV_033",
                      "predictions_by_month.rds",
                      sep = "_"))

write.csv(output, paste(prefixPCres,
                        "TestGroup2",
                        "YV_102",
                        "AV_033",
                        "predictions_by_month.csv",
                        sep = "_"))
    #### Round 3 ####
# Sample from YV: YV_116 
# Sample from AV: AV_023 

models <- models_3
samples <- df_test_3

samples <- samples %>% unnest(pseudocell_all)

meta <- as.data.frame(samples[, c(1:4)])
umi <- samples[, -c(1:4)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)

by_celltype <- as_tibble(cbind(meta,logged)) %>%
  group_by(Cell.Type) %>%
  nest()

by_celltype <- dplyr::inner_join(by_celltype, models, by = "Cell.Type")

by_celltype[1,2][[1]][[1]][1:5,1:5]

custom_add_predictions <- function(data, model) {
  predictions <- predict(model, 
                         newx = as.matrix(data[,-(1:3)]), 
                         s = "lambda.min")
  p <- as_tibble(data.frame("Pred" = predictions[,1],
                            "orig.ident" = data[,1],
                            "Type" = data[,2],
                            "Months" = data[,3]))
  return(p)
}

# Add predictions 
by_celltype <- by_celltype %>%
  mutate(Predictions = map2(data, model, custom_add_predictions))

output <- by_celltype %>% select(-c(data, model)) %>% unnest(Predictions)

range(output$Pred) # 3.896581 17.395502
head(output)
# A tibble: 6 × 5
# Groups:   Cell.Type [1]
# Cell.Type  Pred orig.ident Type      Months
# <fct>     <dbl> <fct>      <fct>      <dbl>
# 1 PERI       12.0 YV_116     Young Veh      2
# 2 PERI       13.3 YV_116     Young Veh      2
# 3 PERI       12.9 YV_116     Young Veh      2
# 4 PERI       13.4 YV_116     Young Veh      2
# 5 PERI       12.8 YV_116     Young Veh      2
# 6 PERI       13.0 YV_116     Young Veh      2

saveRDS(output, paste(prefixPCres,
                      "TestGroup3",
                      "YV_116",
                      "AV_023",
                      "predictions_by_month.rds",
                      sep = "_"))

write.csv(output, paste(prefixPCres,
                        "TestGroup3",
                        "YV_116",
                        "AV_023",
                        "predictions_by_month.csv",
                        sep = "_"))

    #### Round 4 ####
# Sample from YV: YV_106 
# Sample from AV: AV_011 

models <- models_4
samples <- df_test_4

samples <- samples %>% unnest(pseudocell_all)

meta <- as.data.frame(samples[, c(1:4)])
umi <- samples[, -c(1:4)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)

by_celltype <- as_tibble(cbind(meta,logged)) %>%
  group_by(Cell.Type) %>%
  nest()

by_celltype <- dplyr::inner_join(by_celltype, models, by = "Cell.Type")

by_celltype[1,2][[1]][[1]][1:5,1:5]

custom_add_predictions <- function(data, model) {
  predictions <- predict(model, 
                         newx = as.matrix(data[,-(1:3)]), 
                         s = "lambda.min")
  p <- as_tibble(data.frame("Pred" = predictions[,1],
                            "orig.ident" = data[,1],
                            "Type" = data[,2],
                            "Months" = data[,3]))
  return(p)
}

# Add predictions 
by_celltype <- by_celltype %>%
  mutate(Predictions = map2(data, model, custom_add_predictions))

output <- by_celltype %>% select(-c(data, model)) %>% unnest(Predictions)

range(output$Pred) # 1.11135 18.94375
head(output)
# A tibble: 6 × 5
# Groups:   Cell.Type [1]
# Cell.Type  Pred orig.ident Type      Months
# <fct>     <dbl> <fct>      <fct>      <dbl>
# 1 EXCN      10.8  YV_106     Young Veh      3
# 2 EXCN       9.22 YV_106     Young Veh      3
# 3 EXCN      10.7  YV_106     Young Veh      3
# 4 EXCN       8.08 YV_106     Young Veh      3
# 5 EXCN       9.57 YV_106     Young Veh      3
# 6 EXCN      10.3  YV_106     Young Veh      3

saveRDS(output, paste(prefixPCres,
                      "TestGroup4",
                      "YV_106",
                      "AV_011",
                      "predictions_by_month.rds",
                      sep = "_"))

write.csv(output, paste(prefixPCres,
                        "TestGroup4",
                        "YV_106",
                        "AV_011",
                        "predictions_by_month.csv",
                        sep = "_"))

    #### Round 5 ####
# Sample from YV: YV_111 
# Sample from AV: AV_026 

models <- models_5
samples <- df_test_5

samples <- samples %>% unnest(pseudocell_all)

meta <- as.data.frame(samples[, c(1:4)])
umi <- samples[, -c(1:4)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)

by_celltype <- as_tibble(cbind(meta,logged)) %>%
  group_by(Cell.Type) %>%
  nest()

by_celltype <- dplyr::inner_join(by_celltype, models, by = "Cell.Type")

by_celltype[1,2][[1]][[1]][1:5,1:5]

custom_add_predictions <- function(data, model) {
  predictions <- predict(model, 
                         newx = as.matrix(data[,-(1:3)]), 
                         s = "lambda.min")
  p <- as_tibble(data.frame("Pred" = predictions[,1],
                            "orig.ident" = data[,1],
                            "Type" = data[,2],
                            "Months" = data[,3]))
  return(p)
}

# Add predictions 
by_celltype <- by_celltype %>%
  mutate(Predictions = map2(data, model, custom_add_predictions))

output <- by_celltype %>% select(-c(data, model)) %>% unnest(Predictions)

range(output$Pred) # -0.2998924 19.6494343
head(output)
# A tibble: 6 × 5
# Groups:   Cell.Type [1]
# Cell.Type  Pred orig.ident Type      Months
# <fct>     <dbl> <fct>      <fct>      <dbl>
# 1 EXCN      15.0  AV_026     Aging Veh     14
# 2 EXCN      13.8  AV_026     Aging Veh     14
# 3 EXCN      14.3  AV_026     Aging Veh     14
# 4 EXCN      15.5  AV_026     Aging Veh     14
# 5 EXCN       9.52 AV_026     Aging Veh     14
# 6 EXCN      13.4  AV_026     Aging Veh     14
 
saveRDS(output, paste(prefixPCres,
                      "TestGroup5",
                      "YV_111",
                      "AV_026",
                      "predictions_by_month.rds",
                      sep = "_"))

write.csv(output, paste(prefixPCres,
                        "TestGroup5",
                        "YV_111",
                        "AV_026",
                        "predictions_by_month.csv",
                        sep = "_"))

    #### Batch Cross Validation Viz and Stats ####
# load data and put all rounds into a single tibble
b1_val <- readRDS(paste(prefixPCres,
                        "TestGroup1",
                        "YV_104",
                        "AV_029",
                        "predictions_by_month.rds",
                        sep = "_"))
b2_val <- readRDS(paste(prefixPCres,
                        "TestGroup2",
                        "YV_102",
                        "AV_033",
                        "predictions_by_month.rds",
                        sep = "_"))
b3_val <- readRDS(paste(prefixPCres,
                        "TestGroup3",
                        "YV_116",
                        "AV_023",
                        "predictions_by_month.rds",
                        sep = "_"))
b4_val <- readRDS(paste(prefixPCres,
                        "TestGroup4",
                        "YV_106",
                        "AV_011",
                        "predictions_by_month.rds",
                        sep = "_"))
b5_val <- readRDS(paste(prefixPCres,
                        "TestGroup5",
                        "YV_111",
                        "AV_026",
                        "predictions_by_month.rds",
                        sep = "_"))

full_df <- rbind(b1_val,
                 b2_val,
                 b3_val,
                 b4_val,
                 b5_val)

d <- full_df
sumStats <- c()
for (cell in unique(d$Cell.Type)) {
  print(cell)
  cell_df <- filter(d, Cell.Type == cell)
  mae <- round(median(abs(cell_df$Months - cell_df$Pred)), 3)
  meanae <- round(mean(abs(cell_df$Months - cell_df$Pred)), 3)
  rho <- round(cor(cell_df$Months, cell_df$Pred, method = "spearman"), 3)
  r2 <- round(cor(cell_df$Months, cell_df$Pred, method = "pearson") ** 2, 3)
  print(paste0("MedianAbsErr : ", mae))
  print(paste0("MeanAbsErr : ", meanae))
  print(paste0("r2 : ", r2))
  print(paste0("spearman's rho : ", rho))
  print("------------------------- ")
  sumStats <- rbind(sumStats, c(cell, mae, meanae, rho, r2))
}

sumStats <- as.data.frame(sumStats)
colnames(sumStats) <- c("Celltype", "MedianAbsErr", "MeanAbsErr", "rho", "r2")
sumStats

#    Celltype MedianAbsErr MeanAbsErr    rho    r2
# 1      EXCN        2.465      3.462  0.599 0.449
# 2       OPC        2.271      2.605  0.799 0.798
# 3      INHN        2.784       3.74   0.54 0.403
# 4        MG        1.912      2.504  0.763 0.742
# 5      ASTR        1.341      1.815  0.785 0.857
# 6       PYR        1.748      2.765  0.659 0.581
# 7      CPCa         2.43      3.125  0.851 0.703
# 8       GRN        1.172      1.551  0.791  0.89
# 9      OLIG         1.09      1.625  0.773 0.863
# 10     CPCb        5.512      6.152 -0.311 0.173
# 11     BSKT        2.211      2.906  0.683   0.6
# 12     VLMC        3.744      4.097  0.572   0.5
# 13     MOSS        4.304      4.332  0.487 0.387
# 14     PERI        4.068      4.917  0.192 0.173
# 15  ASTR-A1        4.587      4.279  0.546 0.389

# Stats on the median values

d2 <- d %>% group_by(Cell.Type, Months) %>% mutate(med = median(Pred)) %>% select(-Pred) %>% distinct()
sumStats <- c()
for (cell in unique(d2$Cell.Type)) {
  print(cell)
  cell_df <- filter(d2, Cell.Type == cell)
  mae <- round(median(abs(cell_df$Months - cell_df$med)), 3)
  meanae <- round(mean(abs(cell_df$Months - cell_df$med)), 3)
  rho <- round(cor(cell_df$Months, cell_df$med, method = "spearman"), 3)
  r2 <- round(cor(cell_df$Months, cell_df$med, method = "pearson") ** 2, 3)
  print(paste0("MedianAbsErr : ", mae))
  print(paste0("MeanAbsErr : ", meanae))
  print(paste0("r2 : ", r2))
  print(paste0("spearman's rho : ", rho))
  print("------------------------- ")
  sumStats <- rbind(sumStats, c(cell, mae, meanae, rho, r2))
}

sumStats <- as.data.frame(sumStats)
colnames(sumStats) <- c("Celltype", "MedianAbsErr", "MeanAbsErr", "rho", "r2")
sumStats

#    Celltype MedianAbsErr MeanAbsErr    rho    r2
# 1      EXCN        2.463      2.972  0.862 0.692
# 2       OPC        2.166      2.404  0.862 0.937
# 3      INHN        2.874      3.379  0.862  0.68
# 4        MG        1.734      2.175  0.862  0.88
# 5      ASTR        1.201       1.59  0.862 0.914
# 6       PYR        1.384      2.135  0.862 0.697
# 7      CPCa        0.558      2.448      1 0.974
# 8       GRN        0.799      1.082  0.862 0.968
# 9      OLIG        0.629      1.068  0.862 0.935
# 10     CPCb        3.324      6.102 -0.547 0.642
# 11     BSKT        2.097      2.676  0.862  0.74
# 12     VLMC        3.741      4.082  0.862 0.667
# 13     MOSS        5.102      4.104  0.724 0.765
# 14     PERI        5.206      4.954  0.241 0.323
# 15  ASTR-A1        5.368      4.561  0.862 0.881

celltype_order <- c("INHN","BSKT","MOSS",
                    "EXCN",
                    "PYR",
                    "GRN",
                    "MG",
                    "ASTR","ASTR-A1",
                    "OLIG","OPC",
                    "VLMC","CPCa","CPCb","PERI")
d$Cell.Type <- factor(d$Cell.Type, levels = celltype_order)
d2$Cell.Type <- factor(d2$Cell.Type, levels = celltype_order)

ggplot(d, aes(x = Months, y = Pred)) +
  facet_wrap(Cell.Type~.) +
  theme_classic() +
  geom_point(data = d2, aes(x=Months, y=med), 
             size = 2, color = "black") +
  geom_smooth(data = d2, aes(x=Months, y=med), 
              method = "lm", color = "black", linewidth = 1) +
  stat_cor(data = d2, aes(x=Months, y=med), 
           label.x = 8, label.y = 2, size = 2.5) # pearson r not rho

hirestiff(paste(prefixPCres,"Cross","Validation","Pearson","PredictedAges",
                "BootstrapChronological","Months","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"Cross","Validation","Pearson","PredictedAges",
                 "BootstrapChronological","Months","lowres.tiff", sep = "_"))

    #### retrain with all 5 YV and AV ####

# subset out all YV & AV samples
df_YV_AV <- df1 %>% filter(Type %in% c("Young Veh", "Aging Veh"))

df_YV_AV <- droplevels(df_YV_AV) # remove unused factors

# drop the orig.ident column since we're calculating by Cell.Type within Days
df_YV_AV <- df_YV_AV %>% 
  ungroup %>%
  select(-c(orig.ident))

# function for log normalizing
lognorm <- function(input) {
  norm <- sweep(input, MARGIN = 1, FUN = "/", STATS = rowSums(input))
  log1p(norm * 10000)
}

# log normalize, drop pseudocell data, unnest the nested tibble
df2 <- df_YV_AV %>% 
  mutate(lognorm = map(pseudocell_all, lognorm))
df2$pseudocell_all <- NULL
df2 <- unnest(df2, lognorm)

df2 <- df2[ ,-1] # drop Type column

# nest by Cell.Type
by_celltype <- df2 %>% group_by(Cell.Type) %>% nest()

colnames(by_celltype)[2] <- "lognormalized"

# Model function
celltype_model <- function(input) {
  cv.glmnet(x = as.matrix(input[, -1]), 
            y = as.matrix(input[, 1]),
            type.measure = "mae", 
            standardize = F, 
            relax = F, 
            nfolds = 5)
}

# Apply model function using map()
models <- by_celltype %>% 
  mutate(model = map(lognormalized, celltype_model))

models$lognormalized <- NULL

saveRDS(models, paste(prefixPCres,
                      "models_YV_AV_All_by_month.rds",
                      sep = "_"))

models <- readRDS(paste(prefixPCres,
                        "models_YV_AV_All_by_month.rds",
                        sep = "_"))

    #### Get gene weights ####

gene_importances <- tibble()

for (i in c(1:length(models$Cell.Type))) {
  # get Cell.Type
  ct <- as.character(models[i,1][[1]][[1]])
  
  # get glmnet object
  lasso <- models[i,2][[1]][[1]]
  
  # get each gene importance and filter out zeros
  c1 = vip::vi_model(lasso$glmnet.fit, s = lasso$lambda.min, method = "shap")
  d <- c1 %>% filter(Importance != 0)
  
  # make a column to add the correct sign to the importance value
  d$ImportanceSigned <- d$Importance
  d$ImportanceSigned[d$Sign == 'NEG'] <- -d$ImportanceSigned[d$Sign == 'NEG']
  
  d$Cell.Type <- ct

  # add celltype-specific object to cumulative object
  
  gene_importances <- rbind(gene_importances, d)
}

write.csv(gene_importances, 
          paste(prefixPCres,
                "gene_importances",
                "by_month.csv",
                sep = "_"))

    #### Overlapping / Cell.Type Specific Clock Genes ####
      #### sub out important genes based on Cell.Type and Sign ####

INHN_aging_genes <- c()
INHN_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="INHN" & 
                                           gene_importances$Sign=="POS",]$Variable
INHN_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="INHN" & 
                                           gene_importances$Sign=="NEG",]$Variable
BSKT_aging_genes <- c()
BSKT_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="BSKT" & 
                                       gene_importances$Sign=="POS",]$Variable
BSKT_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="BSKT" & 
                                       gene_importances$Sign=="NEG",]$Variable
MOSS_aging_genes <- c()
MOSS_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="MOSS" & 
                                       gene_importances$Sign=="POS",]$Variable
MOSS_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="MOSS" & 
                                       gene_importances$Sign=="NEG",]$Variable
EXCN_aging_genes <- c()
EXCN_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="EXCN" & 
                                       gene_importances$Sign=="POS",]$Variable
EXCN_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="EXCN" & 
                                       gene_importances$Sign=="NEG",]$Variable
PYR_aging_genes <- c()
PYR_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="PYR" & 
                                      gene_importances$Sign=="POS",]$Variable
PYR_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="PYR" & 
                                      gene_importances$Sign=="NEG",]$Variable
GRN_aging_genes <- c()
GRN_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="GRN" & 
                                      gene_importances$Sign=="POS",]$Variable
GRN_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="GRN" & 
                                      gene_importances$Sign=="NEG",]$Variable
MG_aging_genes <- c()
MG_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="MG" & 
                                     gene_importances$Sign=="POS",]$Variable
MG_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="MG" & 
                                     gene_importances$Sign=="NEG",]$Variable
ASTR_aging_genes <- c()
ASTR_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="ASTR" & 
                                       gene_importances$Sign=="POS",]$Variable
ASTR_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="ASTR" & 
                                       gene_importances$Sign=="NEG",]$Variable
ASTR_A1_aging_genes <- c()
ASTR_A1_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="ASTR-A1" & 
                                          gene_importances$Sign=="POS",]$Variable
ASTR_A1_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="ASTR-A1" & 
                                          gene_importances$Sign=="NEG",]$Variable
OLIG_aging_genes <- c()
OLIG_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="OLIG" & 
                                       gene_importances$Sign=="POS",]$Variable
OLIG_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="OLIG" & 
                                       gene_importances$Sign=="NEG",]$Variable
OPC_aging_genes <- c()
OPC_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="OPC" & 
                                      gene_importances$Sign=="POS",]$Variable
OPC_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="OPC" & 
                                      gene_importances$Sign=="NEG",]$Variable
VLMC_aging_genes <- c()
VLMC_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="VLMC" & 
                                       gene_importances$Sign=="POS",]$Variable
VLMC_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="VLMC" & 
                                       gene_importances$Sign=="NEG",]$Variable
CPCa_aging_genes <- c()
CPCa_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="CPCa" & 
                                       gene_importances$Sign=="POS",]$Variable
CPCa_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="CPCa" & 
                                       gene_importances$Sign=="NEG",]$Variable
CPCb_aging_genes <- c()
CPCb_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="CPCb" & 
                                       gene_importances$Sign=="POS",]$Variable
CPCb_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="CPCb" & 
                                       gene_importances$Sign=="NEG",]$Variable
PERI_aging_genes <- c()
PERI_aging_genes$POS <- gene_importances[gene_importances$Cell.Type=="PERI" & 
                                       gene_importances$Sign=="POS",]$Variable
PERI_aging_genes$NEG <- gene_importances[gene_importances$Cell.Type=="PERI" & 
                                       gene_importances$Sign=="NEG",]$Variable

      #### VENN for overlapping genes ####

# venn diagrams of overlapping genes

# BSKT vs INHN vs MOSS
# up
UP_NEURO <- list(BSKT_aging_genes$POS,INHN_aging_genes$POS,MOSS_aging_genes$POS)
names(UP_NEURO) <- c("BSKT_UP","INHN_UP","MOSS_UP")

vennup <- Venn(UP_NEURO)

plot(vennup, doWeights = TRUE, type = "circles")

overlappingup <- vennup@IntersectionSets$`111`
overlappingup <- sort(overlappingup)
overlappingup

write.csv(overlappingup, file = paste(prefixPCres,"BSKT_INHN_MOSS","overlapping","UP","DEGs.csv", sep = "_"))

# down
DN_NEURO <- list(BSKT_aging_genes$NEG,INHN_aging_genes$NEG,MOSS_aging_genes$NEG)
names(DN_NEURO) <- c("BSKT_DN","INHN_DN","MOSS_DN")

venndn <- Venn(DN_NEURO)

plot(venndn, doWeights = TRUE, type = "circles")

overlappingdn <- venndn@IntersectionSets$`111`
overlappingdn <- sort(overlappingdn)
overlappingdn

write.csv(overlappingdn, file = paste(prefixPCres,"BSKT_INHN_MOSS","overlapping","DN","DEGs.csv", sep = "_"))


# ASTR vs ASTR-A1 vs MG
# up
UP_AS_A1_MG <- list(ASTR_aging_genes$POS,ASTR_A1_aging_genes$POS,MG_aging_genes$POS)
names(UP_AS_A1_MG) <- c("ASTR_UP","ASTR_A1_UP","MG_UP")

vennup <- Venn(UP_AS_A1_MG)

plot(vennup, doWeights = TRUE, type = "circles")

overlappingup <- vennup@IntersectionSets$`111`
overlappingup <- sort(overlappingup)
overlappingup

write.csv(overlappingup, file = paste(prefixPCres,"MOSS","MG","overlapping","UP","DEGs.csv", sep = "_"))

vennup@IntersectionSets

# down
DN_AS_A1_MG <- list(ASTR_aging_genes$NEG,ASTR_A1_aging_genes$NEG,MG_aging_genes$NEG)
names(DN_AS_A1_MG) <- c("ASTR_DN","ASTR_A1_DN","MG_DN")

venndn <- Venn(DN_AS_A1_MG)

plot(venndn, doWeights = TRUE, type = "circles")

overlappingdn <- venndn@IntersectionSets$`011`
overlappingdn <- sort(overlappingdn)
overlappingdn

write.csv(overlappingdn, file = paste(prefixPCres,"MOSS","MG","overlapping","DN","DEGs.csv", sep = "_"))



dn_vs_aged <- list(y_v_o_MOSS_dn$symbol,i_v_o_MOSS_dn$symbol)
names(dn_vs_aged) <- c("YoungV_vs_AgedV_MOSS","iMACS_vs_AgedV_MOSS")

venndn <- Venn(dn_vs_aged)


plot(venndn, doWeights = TRUE)



overlappingdn <- venndn@IntersectionSets$`11`
overlappingdn <- sort(overlappingdn)
overlappingdn
write.csv(overlappingdn, file = paste(prefixPCres,"MOSS","overlapping","down","DEGs.csv", sep = "_"))



    #### Predict all ####
df1 <- df1 %>% unnest(pseudocell_all)

meta <- as.data.frame(df1[, c(1:4)])
umi <- df1[, -c(1:4)]
normed <- sweep(umi, MARGIN = 1, FUN = "/", STATS = rowSums(umi))
logged <- log1p(normed * 10000)

by_celltype <- as_tibble(cbind(meta,logged)) %>%
  group_by(Cell.Type) %>%
  nest()

by_celltype <- dplyr::inner_join(by_celltype, models, by = "Cell.Type")

by_celltype[1,2][[1]][[1]][1:5,1:5]

# A tibble: 5 × 5
# orig.ident Type      Months  Xkr4 Gm1992
# <fct>      <fct>      <dbl> <dbl>  <dbl>
# 1 YV_106     Young Veh      3 0          0
# 2 YV_106     Young Veh      3 0          0
# 3 YV_106     Young Veh      3 0          0
# 4 YV_106     Young Veh      3 0.454      0
# 5 YV_106     Young Veh      3 0.458      0

custom_add_predictions <- function(data, model) {
  predictions <- predict(model, 
                         newx = as.matrix(data[,-(1:3)]), 
                         s = "lambda.min")
  p <- as_tibble(data.frame("Pred" = predictions[,1],
                            "orig.ident" = data[,1],
                            "Type" = data[,2],
                            "Months" = data[,3]))
  return(p)
}

# Add predictions 
by_celltype <- by_celltype %>%
  mutate(Predictions = map2(data, model, custom_add_predictions))

output <- by_celltype %>% select(-c(data, model)) %>% unnest(Predictions)

output$lognormalized <- NULL

range(output$Pred) # -0.1774673 21.2516134

head(output)

# A tibble: 6 × 5
# Groups:   Cell.Type [1]
# Cell.Type  Pred orig.ident Type      Months
# <fct>     <dbl> <fct>      <fct>      <dbl>
# 1 EXCN       5.42 YV_106     Young Veh      3
# 2 EXCN       3.04 YV_106     Young Veh      3
# 3 EXCN       4.91 YV_106     Young Veh      3
# 4 EXCN       2.85 YV_106     Young Veh      3
# 5 EXCN       3.66 YV_106     Young Veh      3
# 6 EXCN       4.84 YV_106     Young Veh      3

write.csv(output, paste(prefixPCres,"ALL_predictions_by_month.csv", sep = "_"))

saveRDS(output, paste(prefixPCres,"ALL_predictions_by_month.rds", sep = "_"))

    #### Heatmaps of predicted ages ####
      #### Load data ####
# the "...ALL_predictions_by_month.csv" file was modified slightly in Excel by adding
# a "Cell.Num" column numbered 1 - 100 for each cell type and saved as a csv
dir1 <- "/Users/bells/Library/CloudStorage/Box-Box/GEO Submissions/GEO_Submission_Moser_et_al_2024"
dir2 <- "/20240529_iMPs_15_samples/20240529_cross_validation_and_predictions/"
pred_file <- "20240529_iMPs_15_samples_reclust_13PCs_res.0.6_ALL_predictions_by_month_for_heatmap.csv"

pred_scores <- read.csv(paste0(dir1,dir2,pred_file))

# make a nested tibble to loop through
pred_tibble <- as_tibble(pred_scores) %>% 
  group_by(Cell.Type) %>% 
  nest()

# check nested tibble structure
pred_tibble[1,2][[1]][[1]][1:6]

# access list of cell types
pred_tibble[1]

      #### For one cell type ####

ggplot(
  pred_tibble[1,2][[1]][[1]][1:6],
  aes(
    x = orig.ident,
    y = Cell.Num,
    fill = Pred
  )
) +
  geom_tile() +
  labs(x = "", y = "") +
  guides(fill = guide_legend(title = "Predicted Ages")) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 6
    ),
    axis.text.y = element_text(size = 0),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_x_discrete(limits = rev(
    levels(
      as.factor(
        pred_tibble[1,2][[1]][[1]]$orig.ident
      )
    )
  )) +
  coord_fixed(ratio = 0.5) +
  ggtitle(paste0(pred_tibble[1,1])) +
  scale_fill_viridis_c(
    option = "A",
    direction = -1,
    limits = c(0,20),
    breaks = seq(0,20, by = 5)
  )

      #### A for loop for each cell type, individual plots ####

for (i in 1:length(pred_tibble[1]$Cell.Type)) {
  p1 <- ggplot(pred_tibble[i,2][[1]][[1]][1:6], 
               aes(x = orig.ident,
                   y = Cell.Num,
                   fill = Pred)) + 
    geom_tile() +
    labs(x="", y="") +
    guides(fill=guide_legend(title="Predicted Ages"),
    ) +
    theme(legend.position="right", 
          legend.direction="vertical",
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     size = 6),
          axis.text.y = element_text(size = 0),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()
    ) +
    scale_x_discrete(limits = rev(
      levels(
        as.factor(
          pred_tibble[i,2][[1]][[1]]$orig.ident)))) +
    coord_fixed(ratio = 0.5) +
    ggtitle(paste0(pred_tibble[i,1])) + 
    scale_fill_viridis_c(option = "A", 
                         direction = -1, 
                         limits=c(0, 20), 
                         breaks=seq(0,20,by=5))
  
  print(p1)
  
  hirestiff(paste(prefixPCres,pred_tibble[i,1],"hires.tiff", sep = "_"))
  lowrestiff(paste(prefixPCres,pred_tibble[i,1],"lowres.tiff", sep = "_"))
  
}

      #### All cell types ####

ggplot(
  pred_scores,
  aes(
    x = orig.ident,
    y = Cell.Num,
    fill = Pred
  )
) +
  geom_tile() +
  labs(x = "", y = "") +
  guides(fill = guide_legend(title = "Predicted Ages")) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 6
    ),
    axis.text.y = element_text(size = 0),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_x_discrete(limits = rev(
    levels(
      as.factor(
        pred_scores$orig.ident
      )
    )
  )) +
  scale_fill_viridis_c(
    option = "A",
    direction = -1
  ) + 
  facet_wrap(~ factor(Cell.Type, 
                      levels = c("INHN","BSKT","MOSS",
                                 "EXCN",
                                 "PYR",
                                 "GRN",
                                 "MG",
                                 "ASTR","ASTR-A1",
                                 "OLIG","OPC",
                                 "VLMC","CPCa","CPCb","PERI")
                      ), nrow = 3)

hirestiff(paste(prefixPCres,"heatmap","predicted","ages","ALL","cell","types","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"heatmap","predicted","ages","ALL","cell","types","lowres.tiff", sep = "_"))

      #### Neuronal_Glia cell types ####

pred_scores_sig <- pred_scores[pred_scores$Cell.Type %in% c("ASTR",
                                                            "ASTR-A1",
                                                            "BSKT",
                                                            "EXCN",
                                                            "GRN",
                                                            "INHN",
                                                            "MG",
                                                            "MOSS",
                                                            "PYR"),]

ggplot(
  pred_scores_sig,
  aes(
    x = orig.ident,
    y = Cell.Num,
    fill = Pred
  )
) +
  geom_tile() +
  labs(x = "", y = "") +
  guides(fill = guide_legend(title = "Predicted Ages")) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 6
    ),
    axis.text.y = element_text(size = 0),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_x_discrete(limits = rev(
    levels(
      as.factor(
        pred_scores_sig$orig.ident
      )
    )
  )) +
  scale_fill_viridis_c(
    option = "A",
    direction = -1
  ) + 
  facet_wrap(~ factor(Cell.Type, 
                      levels = c("INHN","BSKT","MOSS",
                                 "EXCN",
                                 "PYR",
                                 "GRN",
                                 "MG",
                                 "ASTR","ASTR-A1")
                      ), nrow = 3)

hirestiff(paste(prefixPCres,"heatmap","predicted","ages","Neuronal_Glia","cell","types","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"heatmap","predicted","ages","Neuronal_Glia","cell","types","lowres.tiff", sep = "_"))

      #### Other cell types ####

pred_scores_other <- pred_scores[pred_scores$Cell.Type %in% c("CPCa",
                                                              "CPCb",
                                                              "OLIG",
                                                              "OPC",
                                                              "PERI",
                                                              "VLMC"),]

ggplot(
  pred_scores_other,
  aes(
    x = orig.ident,
    y = Cell.Num,
    fill = Pred
  )
) +
  geom_tile() +
  labs(x = "", y = "") +
  guides(fill = guide_legend(title = "Predicted Ages")) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 6
    ),
    axis.text.y = element_text(size = 0),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_x_discrete(limits = rev(
    levels(
      as.factor(
        pred_scores_other$orig.ident)))
    ) +
  scale_fill_viridis_c(
    option = "A",
    direction = -1
    ) + 
  facet_wrap(~ factor(Cell.Type, 
                      levels = c("OLIG","OPC","CPCa",
                                 "VLMC","PERI","CPCb")
                      ), nrow = 2)

hirestiff(paste(prefixPCres,"heatmap","predicted","ages","other","cell","types","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"heatmap","predicted","ages","other","cell","types","lowres.tiff", sep = "_"))

#### SUBCLUSTERING ####
  #### sub out Mossy (MOSS) cluster ####
MOSS_cluster <- subset(seurat_sct, idents = c("MOSS"))

DimPlot(MOSS_cluster, 
        reduction = "umap", 
        label = FALSE, 
        group.by = "Type",
        cols = aging_palette) + 
  ggtitle(paste0(ident_res))

hirestiff(paste(prefixPCres,"MOSS","cluster","subset","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"MOSS","cluster","subset","lowres.tiff", sep = "_"))

DimPlot(MOSS_cluster, 
        reduction = "umap", 
        label = FALSE, 
        split.by = "Type", 
        group.by = "Type", 
        cols = aging_palette) + 
  NoLegend()

hirestiff(paste(prefixPCres,"MOSS","cluster","subset","by","type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"MOSS","cluster","subset","by","type","lowres.tiff", sep = "_"))

for (i in 1:length(levels(MOSS_cluster$Type))){
  p <- (DimPlot(MOSS_cluster, 
                reduction = "umap", 
                label = FALSE, 
                cols = aging_palette[i],
                cells = c(WhichCells(MOSS_cluster, expression = Type == levels(MOSS_cluster$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(MOSS_cluster$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","MOSS_cluster",levels(MOSS_cluster$Type)[i],
                  "split_by","type","hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","MOSS_cluster",levels(MOSS_cluster$Type)[i],
                   "split_by","type","lowres.tiff",sep = "_"))
}


    #### SCtransform the subset of interest ####

subprefix_MOSS <- "MOSS_cluster"

MOSS_prefix <- paste0("./",date,"_",project,"_",subprefix_MOSS)

# Use SCTranform for the next steps
MOSS_cluster <- SCTransform(MOSS_cluster,
                            vst.flavor = "v1", 
                            vars.to.regress = c("percent.mt"), 
                            verbose = TRUE)

    #### PCA on the seurat object ####
MOSS_cluster <- RunPCA(MOSS_cluster, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(MOSS_cluster@reductions$pca@cell.embeddings)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(MOSS_cluster@reductions$pca@feature.loadings)

# calculate eigenvalues for arrays

# generate squares of all sample coordinates
sq.xx.coord <- as.data.frame(xx.coord^2)
# create empty list for eigenvalues first
eig <- c()
# calculate the eigenvalue for each PC in sq.xx.coord by taking the sqrt of the sum of squares
for(i in 1:ncol(sq.xx.coord))
  eig[i] = sqrt(sum(sq.xx.coord[,i]))
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if they all contribute equally to the total variance
expected.contribution <- sum.eig/(length(xx.coord)-1)
# return the number of principal components with an eigenvalue greater than expected by equal variance
PCs_MOSS <- sum(eig > expected.contribution)

# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for(i in 1:length(eig))
  eig.percent[i] = 100*eig[i]/sum.eig
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for(i in 1:length(eig))
  if(i == 1) scree[i] = eig.percent[i] else scree[i] = scree[i-1] + eig.percent[i]

# create data frame for eigenvalue summaries
eigenvalues <- data.frame("PC" = colnames(xx.coord), "eig" = eig, "percent" = eig.percent, "scree" = scree)

# plot scree values
plot(eigenvalues$percent, ylim = c(0,100), type = "S", xlab = "PC", ylab = "Percent of variance",
     main = paste0(MOSS_prefix,"scree plot PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = PCs_MOSS, col = "blue")
text(PCs_MOSS, cut.off, label = paste("cutoff PC",PCs_MOSS),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0(subprefix_MOSS,"_scree_plot.pdf"))
dev.off()

# Number of PCs to use:
# PCs_MOSS <- 14

MOSS_prefixpc <- paste0(MOSS_prefix,"_",PCs_MOSS,"PCs")

    #### Run UMAP ####
MOSS_cluster <- RunUMAP(MOSS_cluster, dims = 1:PCs_MOSS, reduction = "pca")

# Plot UMAP                             
## UMAP plot by sample name ("orig.ident")
DimPlot(MOSS_cluster, reduction = "umap", label = FALSE, pt.size = 0.75, group.by = "orig.ident")
hirestiff(paste0(MOSS_prefixpc,"_UMAP_by_sample_hires.tiff"))
lowrestiff(paste0(MOSS_prefixpc,"_UMAP_by_sample_lowres.tiff"))

DimPlot(MOSS_cluster, reduction = "umap", 
        label = FALSE, 
        pt.size = 1.25, 
        group.by = "Type", 
        split.by = "Type", 
        cols = aging_palette) + 
  NoLegend()
hirestiff(paste0(MOSS_prefixpc,"_UMAP_by_type_hires.tiff"))
lowrestiff(paste0(MOSS_prefixpc,"_UMAP_by_type_lowres.tiff"))

DimPlot(MOSS_cluster, reduction = "umap", 
        label = FALSE, 
        pt.size = 1.25, 
        group.by = "Type", 
        split.by = "Batch", 
        cols = aging_palette) + 
  NoLegend()
hirestiff(paste0(MOSS_prefixpc,"_UMAP_by_type_batch_hires.tiff"))
lowrestiff(paste0(MOSS_prefixpc,"_UMAP_by_type_batch_lowres.tiff"))

for (i in 1:length(levels(MOSS_cluster$Type))){
  p <- (DimPlot(MOSS_cluster, 
                reduction = "umap", 
                label = FALSE, 
                cols = aging_palette[i],
                cells = c(WhichCells(MOSS_cluster, expression = Type == levels(MOSS_cluster$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(MOSS_cluster$Type)[i])))
  
  print(p)
  
  hirestiff(paste(MOSS_prefixpc,"UMAP","MOSS_cluster",levels(MOSS_cluster$Type)[i],
                  "split_by","type","hires.tiff",sep = "_"))
  lowrestiff(paste(MOSS_prefixpc,"UMAP","MOSS_cluster",levels(MOSS_cluster$Type)[i],
                   "split_by","type","lowres.tiff",sep = "_"))
}

# saveRDS(MOSS_cluster, file = paste0(MOSS_prefixpc,"_afterUMAP.rds"))

# MOSS_cluster <- read_rds(file = paste0(MOSS_prefixpc,"_afterUMAP.rds"))

    #### Clustering and Resolution ####

# Determine the K-nearest neighbor graph
MOSS_cluster <- FindNeighbors(object = MOSS_cluster, reduction = "pca",  dims = 1:PCs_MOSS)

# Determine the clusters                              
MOSS_cluster <- FindClusters(object = MOSS_cluster,
                             resolution = c(0.4))

subres_MOSS <- "_res.0.4"
ident_subres_MOSS <- paste0("SCT_snn",subres_MOSS)

Idents(MOSS_cluster) <- ident_subres_MOSS
DimPlot(MOSS_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5) + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MOSS))

MOSS_prefixpcres <- paste0(MOSS_prefixpc,subres_MOSS)

DimPlot(MOSS_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5) + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MOSS))

hirestiff(paste(MOSS_prefixpcres,"UMAP","by","cluster","hires.tiff",sep = "_"))
lowrestiff(paste(MOSS_prefixpcres,"UMAP","by","cluster","lowres.tiff",sep = "_"))

DimPlot(MOSS_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 6,
        pt.size = 0.75,
        split.by = "Type",
        repel = TRUE) + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MOSS," split by Type"))

hirestiff(paste(MOSS_prefixpcres,"UMAP","by","cluster","split","by","type","hires.tiff",sep = "_"))
lowrestiff(paste(MOSS_prefixpcres,"UMAP","by","cluster","split","by","type","lowres.tiff",sep = "_"))

DimPlot(MOSS_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        pt.size = 0.75,
        split.by = "orig.ident",
        repel = TRUE) + 
  NoLegend()

DimPlot(MOSS_cluster, 
        reduction = "umap", 
        label = FALSE, 
        label.size = 8,
        pt.size = 0.75,
        split.by = "Type",
        repel = TRUE) + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MOSS," split by Type"))

hirestiff(paste(MOSS_prefixpcres,"UMAP","by","cluster","split","by","type","NO_LABEL","hires.tiff",sep = "_"))
lowrestiff(paste(MOSS_prefixpcres,"UMAP","by","cluster","split","by","type","NO_LABEL","lowres.tiff",sep = "_"))

for (i in 1:length(levels(MOSS_cluster$Type))){
  p <- (DimPlot(MOSS_cluster, 
                reduction = "umap", 
                label = TRUE,
                label.size = 8, 
                pt.size = .75, 
                cells = c(WhichCells(MOSS_cluster, expression = Type == levels(MOSS_cluster$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(MOSS_cluster$Type)[i])))
  
  print(p)
  
  hirestiff(paste(MOSS_prefixpcres,"UMAP","split_by","type",
                  levels(MOSS_cluster$Type)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(MOSS_prefixpcres,"UMAP","split_by","type",
                   levels(MOSS_cluster$Type)[i],"lowres.tiff",sep = "_"))
}

# saveRDS(MOSS_cluster, file = paste0(MOSS_prefixpcres,"_after_reclustering.rds"))

# MOSS_cluster <- readRDS(file = paste0(MOSS_prefixpcres,"_after_reclustering.rds"))

    #### Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(MOSS_cluster, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

n_cells

write.csv(n_cells, file = paste0(MOSS_prefixpcres,"_cells_per_cluster.csv"))

    #### Normalize RNA ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(MOSS_cluster) <- "RNA"

# Normalize, find variable features, scale data 
MOSS_cluster <- NormalizeData(MOSS_cluster, verbose = FALSE)
MOSS_cluster <- FindVariableFeatures(MOSS_cluster)
all.genes <- rownames(MOSS_cluster)
MOSS_cluster <- ScaleData(MOSS_cluster, features = all.genes)

saveRDS(MOSS_cluster, file = paste0(MOSS_prefixpcres,"_after_rna_norm.rds"))

# MOSS_cluster <- readRDS(file = paste0(MOSS_prefixpcres,"_after_rna_norm.rds"))

write.csv(all.genes, 
          file = "20230920_all_MOSS_genes.csv", 
          quote = T, 
          row.names = F)

    #### Cluster markers ####
cluster_markers_MOSS <- FindAllMarkers(MOSS_cluster, 
                                           only.pos = FALSE, 
                                           min.pct = 0.25, 
                                           logfc.threshold = 0.4)
write.csv(cluster_markers_MOSS, 
          file = paste(MOSS_prefixpcres,"cluster_markers",
                       "minpct0.25_logfcthresh0.4.csv",sep = "_"))
    #### Select gene expression ####
expression_data <- FetchData(object = MOSS_cluster, 
                             vars = c("orig.ident","Type",
                                      "Cell.Type","seurat_clusters",
                                      "Vav3","Kcnma1","Robo1","Tcf712",
                                      "Kcnc2","Hdac9","Lef1","Fgf13",
                                      "Rnf220","Cit","Shisa9","Syt9"
                             ))

write.csv(expression_data, 
          file = paste(MOSS_prefixpcres,"select","gene","expression","by","type.csv", sep = "_"))
    #### save MOSS objects ####

# save(subprefix_MOSS,MOSS_cluster,MOSS_prefix,MOSS_prefixpc,
#      MOSS_prefixpcres,PCs_MOSS,subres_MOSS,ident_subres_MOSS,
#      file = paste0(MOSS_prefixpcres,"_MOSS_objects.rdata"))

MOSS_prefixpcres <- paste0("./",date,"_",project,"_","MOSS_cluster","_","14","PCs","_res.0.4")

load(paste0(MOSS_prefixpcres,"_MOSS_objects.rdata"))

#### PROTEOMICS ####
  #### Setup ####
dir <- "/Users/bells/Library/CloudStorage/Box-Box/GEO Submissions/Proteomics Data/moser_et_al_2024/"

library(ggplot2)
library(ggrepel)
library(ggcorrplot)
library(tidyverse)
library(EnhancedVolcano)
library(reshape2)

  #### Heatmap ####

# DEPs were obtained from the AV vs VY comparison in the mapDIA_Statistics tab in the excel file:
# moser_et_al_2024_iMP_Plasma_Proteomics_Raw_and_Calculations_AI_AV_YV.xlsx
# 
# zscores were calculated from the raw data tab "ProteinData" in same file, 
# filtered for DEPs from the AV vs VY comparison, and saved as a csv for loading into R
# 
# The zscore csv was then manually curated for sample order

zfile <- "moser_et_al_2024_iMP_Plasma_Proteomics_ZScore_YV_AV_AI.csv"

zscores <- read.csv(paste0(dir,zfile))

zdata <- reshape2::melt(zscores, id.vars = "Gene")

head(zdata)

ggplot(zdata, aes(x = variable,
                  y = Gene,
                  fill = value)) + 
  geom_tile(
    colour="white", linewidth=0.1
    ) +
  labs(x="", y="") +
  guides(fill=guide_legend(title="Z-Score")) +
  theme(legend.position="right", legend.direction="vertical",
        axis.text.x = element_text(angle = 45,hjust = 1,size = 6),
        axis.text.y = element_text(size = 5)
        ) +
  geom_vline(xintercept = c(20.5, 29.5), color = "black") +
  scale_y_discrete(limits = rev(levels(as.factor(zdata$Gene)))) +
  scale_fill_viridis_c(option = "A", direction = -1) + 
  coord_fixed(ratio = 0.5)

  #### Volcano Plots ####
    #### Read in data, set colors ####
volc <- read_csv(paste0(dir,"moser_et_al_2024_iMP_Plasma_Proteomics_For_Volcano.csv"))

volc_colors <- c("#707070", "#00ff00", "#00ffff", "#ff00ff")

    #### AV vs YV ####
EnhancedVolcano(volc,
                lab = volc$Gene,
                x = 'log2FC_AVvsY',
                y = 'FDR_AVvsY',
                xlim = c(min(volc$log2FC_AVvsY, na.rm = TRUE) - 0.1, 
                         max(volc$log2FC_AVvsY, na.rm = TRUE) + 0.1),
                ylim = c(0, max(-log10(volc$FDR_AVvsY), na.rm = TRUE) + 2),
                title = "Aging Veh vs Young Veh",
                titleLabSize = 30,
                caption = NULL,
                subtitle = NULL,
                xlab = bquote(bold(~Log["2"] ~ "fold change")),
                ylab = bquote(bold(~"-" ~Log["10"] ~ FDR)),
                axisLabSize = 30,
                col = volc_colors,
                colAlpha = 1/1,
                legendLabels = c("NS",
                                 expression(Log[2] ~ FC),
                                 "FDR",
                                 expression(FDR ~ and ~ log[2] ~ FC)),
                legendPosition = 'right',
                legendLabSize = 25,
                legendIconSize = 7,
                drawConnectors = TRUE,
                min.segment.length = 5,
                labSize = 8,
                max.overlaps = Inf,
                selectLab = c("Hba","Adipoq","Gapdh","Plin4",
                              "Cfb","Cfh","C4b","Itih3",
                              "Fgg","Apcs","Saa4","C1qc",
                              "C1qa","Kctd8","Cd5l","Tirap",
                              "Arhgap29","Rgs6","Elp1","Vsnl1",
                              "Chmp5","Ric1"),
                pCutoff = 0.05,
                FCcutoff = 0.58,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 2
                ) 

    #### AI vs YV ####
EnhancedVolcano(volc,
                lab = volc$Gene,
                x = 'log2FC_AMvsY',
                y = 'FDR_AMvsY',
                xlim = c(min(volc$log2FC_AMvsY, na.rm = TRUE) - 0.1, 
                         max(volc$log2FC_AMvsY, na.rm = TRUE) + 0.1),
                ylim = c(0, max(-log10(volc$FDR_AMvsY), na.rm = TRUE) + 2),
                title = "Aging iMPs vs Young Veh",
                titleLabSize = 30,
                caption = NULL,
                subtitle = NULL,
                xlab = bquote(bold(~Log["2"] ~ "fold change")),
                ylab = bquote(bold(~"-" ~Log["10"] ~ FDR)),
                axisLabSize = 30,
                col = volc_colors,
                colAlpha = 1/1,
                legendLabels = c("NS", 
                                 expression(Log[2] ~ FC), 
                                 "FDR", 
                                 expression(FDR ~ and ~ log[2] ~ FC)),
                legendPosition = 'right',
                legendLabSize = 25,
                legendIconSize = 7,
                drawConnectors = TRUE,
                min.segment.length = 5,
                labSize = 8,
                max.overlaps = Inf,
                selectLab = c("Hba","Plin4",
                              "Itih3","Abcb9","Rgs13","Frem2",
                              "Fgg","Apcs","Saa4","C1qc",
                              "C1qa","Kctd8","Cd5l","Tirap",
                              "Rgs6","Elp1",
                              "Chmp5","Ric1","Car1"),
                pCutoff = 0.05,
                FCcutoff = 0.58,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 2
                ) 

    #### AI vs AV ####
EnhancedVolcano(volc,
                lab = volc$Gene,
                x = 'log2FC_AMvsAV',
                y = 'FDR_AMvsAV',
                xlim = c(min(volc$log2FC_AMvsAV, na.rm = TRUE) - 0.1, 
                         max(volc$log2FC_AMvsAV, na.rm = TRUE) + 0.1),
                ylim = c(0, max(-log10(volc$FDR_AMvsAV), na.rm = TRUE) + 2),
                title = "Aging iMPs vs Aging Veh",
                titleLabSize = 30,
                caption = NULL,
                subtitle = NULL,
                xlab = bquote(bold(~Log["2"] ~ "fold change")),
                ylab = bquote(bold(~"-" ~Log["10"] ~ FDR)),
                axisLabSize = 30,
                col = volc_colors,
                colAlpha = 1/1,
                legendLabels = c("NS", 
                                 expression(Log[2] ~ FC), 
                                 "FDR", 
                                 expression(FDR ~ and ~ log[2] ~ FC)),
                legendPosition = 'right',
                legendLabSize = 25,
                legendIconSize = 7,
                drawConnectors = TRUE,
                min.segment.length = 5,
                labSize = 8,
                max.overlaps = Inf,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 2
                )

    #### AV vs AI ####
EnhancedVolcano(volc,
                lab = volc$Gene,
                x = 'log2FC_AVvsAM',
                y = 'FDR_AVvsAM',
                xlim = c(min(volc$log2FC_AVvsAM, na.rm = TRUE) - 0.1, 
                         max(volc$log2FC_AVvsAM, na.rm = TRUE) + 0.1),
                ylim = c(0, max(-log10(volc$FDR_AMvsAV), na.rm = TRUE) + 2),
                title = "Aging Veh vs Aging iMPs",
                titleLabSize = 30,
                caption = NULL,
                subtitle = NULL,
                xlab = bquote(bold(~Log["2"] ~ "fold change")),
                ylab = bquote(bold(~"-" ~Log["10"] ~ FDR)),
                axisLabSize = 30,
                col = volc_colors,
                colAlpha = 1/1,
                legendLabels = c("NS", 
                                 expression(Log[2] ~ FC), 
                                 "FDR", 
                                 expression(FDR ~ and ~ log[2] ~ FC)),
                legendPosition = 'right',
                legendLabSize = 25,
                legendIconSize = 7,
                drawConnectors = TRUE,
                min.segment.length = 5,
                labSize = 8,
                max.overlaps = Inf,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 2
)

#### FIGURES ####
  #### Fig 3A ####
# From the section "UMAPS with cell type labels" under "LABEL CELL TYPES"

DimPlot(seurat_sct,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        pt.size = 0.5,
        cols = cell_cols,
        raster = F) + 
  NoLegend()
hirestiffsquare(paste(prefixPCres,"UMAP","cell","types","hires_square.tiff", sep = "_"))

  #### Fig 3B ####
# From the section "All cell types VLN" under "LABEL CELL TYPES"

VlnPlot(seurat_sct,features = c(
                                "Gad1","Gad2",            ## INHN
                                "Dlx1","Col19a1",         ## BSKT
                                "Cnr1","Calb2",           ## MOSS
                                "Slc17a7","Neurod6",      ## EXCN
                                "Grm1","Satb2",           ## PYR
                                "Calb1","Prox1",          ## GRN
                                "Tmem119","Itgam",        ## MG
                                "Aqp4","Gfap",            ## ASTR
                                "S100b","C3",             ## ASTR-A1
                                "Mbp","Mog",              ## OLIG
                                "Olig1","Sox10",          ## OPC
                                "Pdgfra","Slc6a13",       ## VLMC
                                "Hemk1","Acad8",          ## CPCa
                                "Folr1","Kl",             ## CPCb
                                "Kcnj8","Rgs5"            ## PERI
                                             ),
        stack = TRUE,
        flip = TRUE,
        sort = FALSE) + 
  NoLegend() + 
  ggtitle("All Cell Types") &
  theme(
        axis.text.x = element_text(size = 20), 
        axis.title = element_text(size = 22),
        axis.text.y.left = element_text(size = 12),
        axis.text.y.right = element_text(size = 12),
        strip.text = element_text(size = 15,
                                  face = "bold.italic"))

hirestiffsquare(paste(prefixPCres,"VLN","cell","type","markers","hires_square.tiff", sep = "_"))

  #### Fig 3C ####
# From the section "Clustering and Resolution" under "SUBCLUSTERING"

DimPlot(MOSS_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 6,
        pt.size = 0.75,
        split.by = "Type",
        repel = TRUE) + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MOSS," split by Type"))

hirestiff(paste(MOSS_prefixpcres,"UMAP","by","cluster","split","by","type","hires.tiff",sep = "_"))

  #### Data for Fig 3D, 3E ####
# From the section "Cluster markers" under "SUBCLUSTERING" 
 # Data loaded into ImageGP:  https://www.bic.ac.cn/BIC/
  # Cluster markers used for Fig 3D & 3E

cluster_markers_MOSS <- FindAllMarkers(MOSS_cluster, 
                                           only.pos = FALSE, 
                                           min.pct = 0.25, 
                                           logfc.threshold = 0.4)
write.csv(cluster_markers_MOSS, 
          file = paste(MOSS_prefixpcres,"cluster_markers",
                       "minpct0.25_logfcthresh0.4.csv",sep = "_"))

  #### Data for Fig 3F ####
# From the section "Select gene expression" under "SUBCLUSTERING"
 # Data loaded into Prism

expression_data <- FetchData(object = MOSS_cluster, 
                             vars = c("orig.ident","Type",
                                      "Cell.Type","seurat_clusters",
                                      "Vav3","Kcnma1","Robo1","Tcf712",
                                      "Kcnc2","Hdac9","Lef1","Fgf13",
                                      "Rnf220","Cit","Shisa9","Syt9"
                             ))

write.csv(expression_data, 
          file = paste(MOSS_prefixpcres,"select","gene","expression","by","type.csv", sep = "_"))

  #### Fig 4B ####
# From the section "Heatmaps of predicted ages" -> "Neuronal_Glia cell types" under "AGING CLOCKS"

pred_scores_sig <- pred_scores[pred_scores$Cell.Type %in% c("ASTR",
                                                            "ASTR-A1",
                                                            "BSKT",
                                                            "EXCN",
                                                            "GRN",
                                                            "INHN",
                                                            "MG",
                                                            "MOSS",
                                                            "PYR"),]

ggplot(
  pred_scores_sig,
  aes(
    x = orig.ident,
    y = Cell.Num,
    fill = Pred
  )
) +
  geom_tile() +
  labs(x = "", y = "") +
  guides(fill = guide_legend(title = "Predicted Ages")) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 6
    ),
    axis.text.y = element_text(size = 0),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_x_discrete(limits = rev(
    levels(
      as.factor(
        pred_scores_sig$orig.ident
      )
    )
  )) +
  scale_fill_viridis_c(
    option = "A",
    direction = -1
  ) + 
  facet_wrap(~ factor(Cell.Type, 
                      levels = c("INHN","BSKT","MOSS",
                                 "EXCN",
                                 "PYR",
                                 "GRN",
                                 "MG",
                                 "ASTR","ASTR-A1")
                      ), nrow = 3)

hirestiff(paste(prefixPCres,"heatmap","predicted","ages","Neuronal_Glia","cell","types","hires.tiff", sep = "_"))

  #### Data for Fig 4C ####
# From the section "Age in Months" -> "Predict all" under "AGING CLOCKS"
 # Data loaded into Prism

read.csv(paste(prefixPCres,
               "ALL_predictions_by_month.csv", 
               sep = "_"))

  #### Data for Fig 4D ####
# From the section "Age in Months" -> "Get gene weights" under "AGING CLOCKS"
 # Data loaded into ImageGP: https://www.bic.ac.cn/BIC/
  # Overlapping genes in at least two cell types were manually curated in excel

read.csv(paste(prefixPCres,
               "gene_importances",
               "by_month.csv",
               sep = "_"))


  #### Fig 5A ####
# From the section "Volcano Plots" -> "AI vs AV" under "PROTEOMICS"
EnhancedVolcano(volc,
                lab = volc$Gene,
                x = 'log2FC_AMvsAV',
                y = 'FDR_AMvsAV',
                xlim = c(min(volc$log2FC_AMvsAV, na.rm = TRUE) - 0.1, 
                         max(volc$log2FC_AMvsAV, na.rm = TRUE) + 0.1),
                ylim = c(0, max(-log10(volc$FDR_AMvsAV), na.rm = TRUE) + 2),
                title = "Aging iMPs vs Aging Veh",
                titleLabSize = 30,
                caption = NULL,
                subtitle = NULL,
                xlab = bquote(bold(~Log["2"] ~ "fold change")),
                ylab = bquote(bold(~"-" ~Log["10"] ~ FDR)),
                axisLabSize = 30,
                col = volc_colors,
                colAlpha = 1/1,
                legendLabels = c("NS", 
                                 expression(Log[2] ~ FC), 
                                 "FDR", 
                                 expression(FDR ~ and ~ log[2] ~ FC)),
                legendPosition = 'right',
                legendLabSize = 25,
                legendIconSize = 7,
                drawConnectors = TRUE,
                min.segment.length = 5,
                labSize = 8,
                max.overlaps = Inf,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 2
                )
  #### Fig 5B ####
# From the section "Heatmap" under "PROTEOMICS"

ggplot(zdata, aes(x = variable,
                  y = Gene,
                  fill = value)) + 
  geom_tile(
    colour="white", linewidth=0.1
    ) +
  labs(x="", y="") +
  guides(fill=guide_legend(title="Z-Score")) +
  theme(legend.position="right", legend.direction="vertical",
        axis.text.x = element_text(angle = 45,hjust = 1,size = 6),
        axis.text.y = element_text(size = 5)
        ) +
  geom_vline(xintercept = c(20.5, 29.5), color = "black") +
  scale_y_discrete(limits = rev(levels(as.factor(zdata$Gene)))) +
  scale_fill_viridis_c(option = "A", direction = -1) + 
  coord_fixed(ratio = 0.5)
#### SUPPLEMENTAL FIGURES ####
  #### Supplemental Fig 4A ####
# From the section "UMAPS with cell type labels" under "LABEL CELL TYPES"

for (i in 1:length(levels(seurat_sct$Type))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = TRUE,
                label.size = 6,
                pt.size = 0.5,
                cols = cell_cols,
                raster = F, 
                cells = c(WhichCells(seurat_sct, 
                                     expression = Type == levels(seurat_sct$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct$Type)[i])))
  
  print(p)
  
  hirestiffsquare(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                        levels(seurat_sct$Type)[i],"hires_square.tiff",sep = "_"))
}

  #### Data for Supplemental Fig 4B ####
# From the section "DEGs within cell type clusters" under "LABEL CELL TYPES"
 # Data manually compiled and loaded into Prism

# for loop for DEGs within cell type clusters

# get the names of cell types
clusterlist <- unique(seurat_sct@meta.data$Cell.Type)

# list the groups to compare.  Since we're interested in 
# Young Veh vs Aged Veh, 
# Aged Veh vs Aged iMAC, 
# and Aged iMAC vs Young Veh,
# I added an additional Young Veh to the end of the list.
# The loop will compare the first and second terms,
# then the second and third terms,
# then the third and fourth terms and then will terminate
grouplist <- c("Young Veh","Aging Veh","Aging iMPs","Young Veh")

for (i in 1:(length(grouplist)-1)){
  for (k in 1:length(clusterlist)){
    
    nam <- paste("clusterDEG", grouplist[i],"vs",grouplist[i+1],"cluster",clusterlist[k], sep = "_")
    assign(nam, FindMarkers(seurat_sct, 
                            ident.1 = grouplist[i], 
                            ident.2 = grouplist[i+1],
                            group.by = "Type", 
                            subset.ident = clusterlist[k]))
    
    write.csv(get(nam), file = paste(prefixPCres, nam, "list.csv", sep = "_"), row.names = TRUE)
    
  } 
}

  #### Supplemental Fig 5A ####
# From the section "Batch Cross Validation Viz and Stats" under "AGING CLOCKS"

ggplot(d, aes(x = Months, y = Pred)) +
  facet_wrap(Cell.Type~.) +
  theme_classic() +
  geom_point(data = d2, aes(x=Months, y=med), 
             size = 2, color = "black") +
  geom_smooth(data = d2, aes(x=Months, y=med), 
              method = "lm", color = "black", linewidth = 1) +
  stat_cor(data = d2, aes(x=Months, y=med), 
           label.x = 8, label.y = 2, size = 2.5) # pearson r not rho

hirestiff(paste(prefixPCres,"Cross","Validation","Pearson","PredictedAges",
                "BootstrapChronological","Months","hires.tiff", sep = "_"))

  #### Supplemental Fig 5B ####
# From the section "Heatmaps of predicted ages" -> "Other cell types" under "AGING CLOCKS"

pred_scores_other <- pred_scores[pred_scores$Cell.Type %in% c("CPCa",
                                                              "CPCb",
                                                              "OLIG",
                                                              "OPC",
                                                              "PERI",
                                                              "VLMC"),]

ggplot(
  pred_scores_other,
  aes(
    x = orig.ident,
    y = Cell.Num,
    fill = Pred
  )
) +
  geom_tile() +
  labs(x = "", y = "") +
  guides(fill = guide_legend(title = "Predicted Ages")) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 6
    ),
    axis.text.y = element_text(size = 0),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_x_discrete(limits = rev(
    levels(
      as.factor(
        pred_scores_other$orig.ident)))
    ) +
  scale_fill_viridis_c(
    option = "A",
    direction = -1
    ) + 
  facet_wrap(~ factor(Cell.Type, 
                      levels = c("OLIG","OPC","CPCa",
                                 "VLMC","PERI","CPCb")
                      ), nrow = 2)

hirestiff(paste(prefixPCres,"heatmap","predicted","ages","other","cell","types","hires.tiff", sep = "_"))

  #### Data for Supplemental Fig 5C ####
# from the section "Age in Months" -> "Get gene weights" under "AGING CLOCKS"
 # Data loaded into prism

read.csv(paste(prefixPCres,
               "ALL_predictions_by_month.csv", 
               sep = "_"))

  #### Supplemental Fig 7A ####
# From the section "Volcano Plots" -> "AV vs YV" under "PROTEOMICS"

EnhancedVolcano(volc,
                lab = volc$Gene,
                x = 'log2FC_AVvsY',
                y = 'FDR_AVvsY',
                xlim = c(min(volc$log2FC_AVvsY, na.rm = TRUE) - 0.1, 
                         max(volc$log2FC_AVvsY, na.rm = TRUE) + 0.1),
                ylim = c(0, max(-log10(volc$FDR_AVvsY), na.rm = TRUE) + 2),
                title = "Aging Veh vs Young Veh",
                titleLabSize = 30,
                caption = NULL,
                subtitle = NULL,
                xlab = bquote(bold(~Log["2"] ~ "fold change")),
                ylab = bquote(bold(~"-" ~Log["10"] ~ FDR)),
                axisLabSize = 30,
                col = volc_colors,
                colAlpha = 1/1,
                legendLabels = c("NS",
                                 expression(Log[2] ~ FC),
                                 "FDR",
                                 expression(FDR ~ and ~ log[2] ~ FC)),
                legendPosition = 'right',
                legendLabSize = 25,
                legendIconSize = 7,
                drawConnectors = TRUE,
                min.segment.length = 5,
                labSize = 8,
                max.overlaps = Inf,
                selectLab = c("Hba","Adipoq","Gapdh","Plin4",
                              "Cfb","Cfh","C4b","Itih3",
                              "Fgg","Apcs","Saa4","C1qc",
                              "C1qa","Kctd8","Cd5l","Tirap",
                              "Arhgap29","Rgs6","Elp1","Vsnl1",
                              "Chmp5","Ric1"),
                pCutoff = 0.05,
                FCcutoff = 0.58,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 2
                ) 

  #### Supplemental Fig 7B ####
# From the section "Volcano Plots" -> "AI vs YV" under "PROTEOMICS"

EnhancedVolcano(volc,
                lab = volc$Gene,
                x = 'log2FC_AMvsY',
                y = 'FDR_AMvsY',
                xlim = c(min(volc$log2FC_AMvsY, na.rm = TRUE) - 0.1, 
                         max(volc$log2FC_AMvsY, na.rm = TRUE) + 0.1),
                ylim = c(0, max(-log10(volc$FDR_AMvsY), na.rm = TRUE) + 2),
                title = "Aging iMPs vs Young Veh",
                titleLabSize = 30,
                caption = NULL,
                subtitle = NULL,
                xlab = bquote(bold(~Log["2"] ~ "fold change")),
                ylab = bquote(bold(~"-" ~Log["10"] ~ FDR)),
                axisLabSize = 30,
                col = volc_colors,
                colAlpha = 1/1,
                legendLabels = c("NS", 
                                 expression(Log[2] ~ FC), 
                                 "FDR", 
                                 expression(FDR ~ and ~ log[2] ~ FC)),
                legendPosition = 'right',
                legendLabSize = 25,
                legendIconSize = 7,
                drawConnectors = TRUE,
                min.segment.length = 5,
                labSize = 8,
                max.overlaps = Inf,
                selectLab = c("Hba","Plin4",
                              "Itih3","Abcb9","Rgs13","Frem2",
                              "Fgg","Apcs","Saa4","C1qc",
                              "C1qa","Kctd8","Cd5l","Tirap",
                              "Rgs6","Elp1",
                              "Chmp5","Ric1","Car1"),
                pCutoff = 0.05,
                FCcutoff = 0.58,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 2
                ) 
