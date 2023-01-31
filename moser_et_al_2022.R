#!/bin/env Rscript

  #### License Notice ####

##
# Copyright (c) 2022 Cedars-Sinai Medical Center
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
# Title: Human iPSC-derived mononuclear phagocytes restore cognition, 
#        neural health, and a subpopulation of hippocampal mossy cells 
#        in aging mice
#
# Authors:  V. Alexandra Moser 1, Rachel M. Lipman 1, Shaughn Bell 1, 
#           George Lawless 1, Luz Jovita Dimas-Harms 1, Jake Inzalaco 1, 
#           Simion Kreimer 2, Sarah J. Parker 2, Helen S. Goodridge 1, 
#           Clive N. Svendsen**1
# 
# Affiliations: 1 Board of Governors Regenerative Medicine Institute, 
#                 Cedars-Sinai Medical Center, Los Angeles, CA, USA. 
#               2 Smidt Heart Institute, Department of Cardiology, 
#                 Cedars-Sinai Medical Center, Los Angeles, CA, USA.
#               
# **Corresponding author email:  Clive.Svendsen@cshs.org
##

  #### Script Information ####

##
# R version 4.2.2
# R Script Title:  moser_et_al_2022.R
# R Script Author:  Shaughn Bell
# R Script Corresponding Email:  shaughn.bell@cshs.org
#
# Notes: 
#   A) Script makes use of the variables set up under "project information" as 
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

##
# Single Nucleus RNAseq data
#
# GEO Accession Number:  GSE220548
#
# aligned to mus musculus mm10 via CellRanger v6.1.2
# used the 10x Genomics supplied reference file "refdata-gex-mm10-2020-A.tar.gz"
#
# fastq files used:
#
# 106HPC_S5_L001_I1_001.fastq.gz
# 106HPC_S5_L001_R1_001.fastq.gz
# 106HPC_S5_L001_R2_001.fastq.gz
# 106HPC_S5_L002_I1_001.fastq.gz
# 106HPC_S5_L002_R1_001.fastq.gz
# 106HPC_S5_L002_R2_001.fastq.gz
# 116HPC_S7_L001_I1_001.fastq.gz
# 116HPC_S7_L001_R1_001.fastq.gz
# 116HPC_S7_L001_R2_001.fastq.gz
# 116HPC_S7_L002_I1_001.fastq.gz
# 116HPC_S7_L002_R1_001.fastq.gz
# 116HPC_S7_L002_R2_001.fastq.gz
# 29HPC_S3_L001_I1_001.fastq.gz
# 29HPC_S3_L001_R1_001.fastq.gz
# 29HPC_S3_L001_R2_001.fastq.gz
# 29HPC_S3_L002_I1_001.fastq.gz
# 29HPC_S3_L002_R1_001.fastq.gz
# 29HPC_S3_L002_R2_001.fastq.gz
# 33HPC_S2_L001_I1_001.fastq.gz
# 33HPC_S2_L001_R1_001.fastq.gz
# 33HPC_S2_L001_R2_001.fastq.gz
# 33HPC_S2_L002_I1_001.fastq.gz
# 33HPC_S2_L002_R1_001.fastq.gz
# 33HPC_S2_L002_R2_001.fastq.gz
# 42HPC_S8_L001_I1_001.fastq.gz
# 42HPC_S8_L001_R1_001.fastq.gz
# 42HPC_S8_L001_R2_001.fastq.gz
# 42HPC_S8_L002_I1_001.fastq.gz
# 42HPC_S8_L002_R1_001.fastq.gz
# 42HPC_S8_L002_R2_001.fastq.gz
# 17HPC_S1_L001_I1_001.fastq.gz
# 17HPC_S1_L001_R1_001.fastq.gz
# 17HPC_S1_L001_R2_001.fastq.gz
# 17HPC_S1_L002_I1_001.fastq.gz
# 17HPC_S1_L002_R1_001.fastq.gz
# 17HPC_S1_L002_R2_001.fastq.gz
# 31HPC_S4_L001_I1_001.fastq.gz
# 31HPC_S4_L001_R1_001.fastq.gz
# 31HPC_S4_L001_R2_001.fastq.gz
# 31HPC_S4_L002_I1_001.fastq.gz
# 31HPC_S4_L002_R1_001.fastq.gz
# 31HPC_S4_L002_R2_001.fastq.gz
# 45HPC_S6_L001_I1_001.fastq.gz
# 45HPC_S6_L001_R1_001.fastq.gz
# 45HPC_S6_L001_R2_001.fastq.gz
# 45HPC_S6_L002_I1_001.fastq.gz
# 45HPC_S6_L002_R1_001.fastq.gz
# 45HPC_S6_L002_R2_001.fastq.gz
#
# Cellranger outs used:  
#
# mm_106HPC_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# mm_106HPC_outs/filtered_feature_bc_matrix/features.tsv.gz
# mm_106HPC_outs/filtered_feature_bc_matrix/matrix.mtx.gz
# 
# mm_116HPC_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# mm_116HPC_outs/filtered_feature_bc_matrix/features.tsv.gz
# mm_116HPC_outs/filtered_feature_bc_matrix/matrix.mtx.gz
# 
# mm_29HPC_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# mm_29HPC_outs/filtered_feature_bc_matrix/features.tsv.gz
# mm_29HPC_outs/filtered_feature_bc_matrix/matrix.mtx.gz
# 
# mm_33HPC_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# mm_33HPC_outs/filtered_feature_bc_matrix/features.tsv.gz
# mm_33HPC_outs/filtered_feature_bc_matrix/matrix.mtx.gz
# 
# mm_17HPC_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# mm_17HPC_outs/filtered_feature_bc_matrix/features.tsv.gz
# mm_17HPC_outs/filtered_feature_bc_matrix/matrix.mtx.gz
# 
# mm_31HPC_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# mm_31HPC_outs/filtered_feature_bc_matrix/features.tsv.gz
# mm_31HPC_outs/filtered_feature_bc_matrix/matrix.mtx.gz
# 
# mm_45HPC_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# mm_45HPC_outs/filtered_feature_bc_matrix/features.tsv.gz
# mm_45HPC_outs/filtered_feature_bc_matrix/matrix.mtx.gz
# 
# Cellranger outs are not included in the GEO record.
#
# Genesummed matricies of counts:  
# 
# 1_YV_106_filtered_genesummed_mtx.csv
# 2_YV_116_filtered_genesummed_mtx.csv
# 3_AV_029_filtered_genesummed_mtx.csv
# 4_AV_033_filtered_genesummed_mtx.csv
# 5_AV_042_filtered_genesummed_mtx.csv
# 6_AI_017_filtered_genesummed_mtx.csv
# 7_AI_031_filtered_genesummed_mtx.csv
# 8_AI_045_filtered_genesummed_mtx.csv
#
# Genesummed matrix files are included in the GEO record. 
##

##
# Proteomics Data (included in GitHub record)
# 
# Raw data and calculations:  
#
# moser_et_al_2022_iMP_Plasma_Proteomics_Raw_and_Calculations.xlsx
#
# Raw data in csv format:
#
# moser_et_al_2022_iMP_Plasma_Proteomics_Raw.csv
#
# Differentially expressed proteins as described in the proteomics section:
#
# moser_et_al_2022_iMP_Plasma_Proteomics_deps.csv
#
# Raw csv spreadsheet that only includes the differentially expressed proteins:
#
# moser_et_al_2022_iMP_Plasma_Proteomics_Raw_Filtered_DEPs.csv
#
# Formatted and filtered with outliers removed for volcano plots:
#
# moser_et_al_2022_iMP_Plasma_Proteomics_For_Volcano.csv
#
# Z Score spreadsheet from Raw and Caluclations excel:
# 
# moser_et_al_2022_iMP_Plasma_Proteomics_ZScore.csv
##

  #### project information ####

date <- "20220923"
project <- "iMPs_all"
datadir <- "/20220927_iMPs_all" # directory to save files generated
sourcedir <- "/20220523_mm_outs_copy" # save Cellranger outs here

  #### Load required packages ####

library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(ggcorrplot)
library(reshape2)
library(cowplot)
library(patchwork)
library(irlba)
library(gridExtra)
library(scCustomize)
library(nichenetr)
library(EnhancedVolcano)
library(Vennerable)

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

#### ####
  #### Gene Sum ####

# Some gene symbols refer to more than one ENSG ID.  
# In order to ensure all ENSG IDs are represented while using the much more
# user-friendly gene symbol, the cellranger outputs were processed as follows:
#    1) Cellranger output was loaded into a matrix object
#    2) All genes with more than one ENSG ID were put into a new sub-matrix
#    3) Counts from all instances of a given gene symbol were summed
#    4) The genesummed sub-matrix was rejoined to the main matrix
#    5) The ENSG ID column was removed
# Genesummed matrix files are included in the GEO record.  
# Cellranger outs are not included in the GEO record.

# Create matrix files from the cellranger outs

setwd(sourcedir)

filelist <- list.files(path = sourcedir)
filelist <- as.data.frame(filelist, files = list.files(path = sourcedir))
path_list <- paste0(sourcedir,"/",filelist[,1],"/filtered_feature_bc_matrix")
### IMPORTANT ###
# Files with numbers in their name may be sorted in a 
# different order depending on the way files are listed.  
# Some systems sort by number instead of numeric value.
# Make sure to use the same sort order as in the path_list object.
path_list

all(file.exists(path_list)) # check that all the files exist

#Read the cellranger output files and save expression matrix csv
for (k in 1:length(path_list)){
  
  # read in and transpose the barcode file
  bc <- as.data.frame(read.delim(file = paste0(path_list[k],"/barcodes.tsv.gz"), header = F))%>%t(.)
  # read in the features file and select just the first two columns
  feat <- as.data.frame(read.delim(file = paste0(path_list[k],"/features.tsv.gz"), header = F))%>%.[,1:2]
  # read in the counts
  matx <- as.data.frame(readMM(file = paste0(path_list[k],"/matrix.mtx.gz")))
  
  colnames(matx) <- bc[1,] #assign barcodes as colnames for the matrix
  colnames(feat) <- c("ensmbl.id", "symbol") #rename the columns of features table
  matx <- cbind(feat, matx) #cbind the features table and the matrix
  
  write.csv(matx, file = paste0(path_list[k],"/filtered_combined_mtx.csv"), row.names = FALSE)
  
}

rm(matx,bc,feat,filelist,k)

#### Load the data files in to "data" objects

# !!!!! check that the path list is in the same sample order !!!!!
path_list

data.106 <- read.csv(file=paste0(path_list[1],"/filtered_combined_mtx.csv"), header = TRUE)
data.116 <- read.csv(file=paste0(path_list[2],"/filtered_combined_mtx.csv"), header = TRUE)
data.017 <- read.csv(file=paste0(path_list[3],"/filtered_combined_mtx.csv"), header = TRUE)
data.029 <- read.csv(file=paste0(path_list[4],"/filtered_combined_mtx.csv"), header = TRUE)
data.031 <- read.csv(file=paste0(path_list[5],"/filtered_combined_mtx.csv"), header = TRUE)
data.033 <- read.csv(file=paste0(path_list[6],"/filtered_combined_mtx.csv"), header = TRUE)
data.042 <- read.csv(file=paste0(path_list[7],"/filtered_combined_mtx.csv"), header = TRUE)
data.045 <- read.csv(file=paste0(path_list[8],"/filtered_combined_mtx.csv"), header = TRUE)

# put all the loaded data into one list to iterate over
obj.list <- list(data.106,data.116,
                 data.017,data.029,
                 data.031,data.033,
                 data.042,data.045)

# iterate the gene sum loop over all samples in obj.list
for (j in 1:length(obj.list)){
  
  #create a list of duplicated genes
  gene.duplicates<-obj.list[[j]]$symbol[duplicated(obj.list[[j]]$symbol, incomparables = NA)]%>%as.list(.)%>%unique(.)
  
  #Create empty dataframe to rbind duplicated rows to
  dup.rows<-data.frame()
  
  #Create identical expression matrix to subtract duplicated rows from
  sub.rows<-obj.list[[j]]
  
  #Extract rows containing the duplicated genes. 
  for (i in 1:length(gene.duplicates)){
    dups<-subset(obj.list[[j]], subset = symbol==gene.duplicates[[i]])#subsets duplicated rows
    dup.rows<-rbind(dup.rows, dups)#binds duplicated rows together
    sub.rows<-subset(sub.rows, subset = symbol!=gene.duplicates[[i]])#deleted duplicated rows from expression matrix
  }
  
  #Remove ensembl ID columns
  dup.rows<-dup.rows[,2:length(colnames(dup.rows))]
  sub.rows<-sub.rows[,2:length(colnames(sub.rows))]
  
  #Find the gsum of duplicated rows
  dup.rows<-aggregate(. ~ symbol, data = dup.rows, sum)
  
  #rbind the gsum of duplicated rows with the expression matrix that had these genes removed
  new.mat<-rbind(sub.rows, dup.rows)
  
  if (nrow(new.mat) == length(unique(obj.list[[j]]$symbol))) {
    new.mat <- as.data.frame(new.mat)
    rownames(new.mat) <- new.mat[,1]
    new.mat <- subset(new.mat, select = -c(symbol))
    obj.list[[j]] <- new.mat
    message("Successfully merged!")
    } else {
      message("There's an error!")
      }
}

# assign the gene summed items in obj.list back to their original data object
data.106 <- obj.list[[1]]
data.116 <- obj.list[[2]]
data.017 <- obj.list[[3]]
data.029 <- obj.list[[4]]
data.031 <- obj.list[[5]]
data.033 <- obj.list[[6]]
data.042 <- obj.list[[7]]
data.045 <- obj.list[[8]]

# save gene summed matrix objects

# !!!!! check that the path list is in the same sample order !!!!!
path_list

write.csv(data.106, file = paste0(path_list[1],"/","YV_106","_filtered_genesummed_mtx.csv"), row.names = FALSE)
write.csv(data.116, file = paste0(path_list[2],"/","YV_116","_filtered_genesummed_mtx.csv"), row.names = FALSE)
write.csv(data.017, file = paste0(path_list[3],"/","AI_017","_filtered_genesummed_mtx.csv"), row.names = FALSE)
write.csv(data.029, file = paste0(path_list[4],"/","AV_029","_filtered_genesummed_mtx.csv"), row.names = FALSE)
write.csv(data.031, file = paste0(path_list[5],"/","AI_031","_filtered_genesummed_mtx.csv"), row.names = FALSE)
write.csv(data.033, file = paste0(path_list[6],"/","AV_033","_filtered_genesummed_mtx.csv"), row.names = FALSE)
write.csv(data.042, file = paste0(path_list[7],"/","AV_042","_filtered_genesummed_mtx.csv"), row.names = FALSE)
write.csv(data.045, file = paste0(path_list[8],"/","AI_045","_filtered_genesummed_mtx.csv"), row.names = FALSE)

# save R objects for easy reloading
save(data.106,data.116,
     data.017,data.029,
     data.031,data.033,
     data.042,data.045,
     file = paste0(datadir,"/",date,"_",project,"_genesummed.rdata"))

# load(file = paste0(datadir,"/",date,"_",project,"_genesummed.rdata"))

  #### Seurat Object ####

# set wd to datadir
setwd(datadir)

# Initialize the Seurat object with the raw (non-normalized data).
# Keep all genes expressed in >= 1 cell

s106 <- CreateSeuratObject(counts = obj.list[[1]], project = "YV-106", min.cells = 1, min.features = 0)
s116 <- CreateSeuratObject(counts = obj.list[[2]], project = "YV-116", min.cells = 1, min.features = 0)
s017 <- CreateSeuratObject(counts = obj.list[[3]], project = "AI-017", min.cells = 1, min.features = 0)
s029 <- CreateSeuratObject(counts = obj.list[[4]], project = "AV-029", min.cells = 1, min.features = 0)
s031 <- CreateSeuratObject(counts = obj.list[[5]], project = "AI-031", min.cells = 1, min.features = 0)
s033 <- CreateSeuratObject(counts = obj.list[[6]], project = "AV-033", min.cells = 1, min.features = 0)
s042 <- CreateSeuratObject(counts = obj.list[[7]], project = "AV-042", min.cells = 1, min.features = 0)
s045 <- CreateSeuratObject(counts = obj.list[[8]], project = "AI-045", min.cells = 1, min.features = 0)

rm(data.106,data.116,
   data.017,data.029,
   data.031,data.033,
   data.042,data.045,
   dup.rows,dups,gene.duplicates,new.mat,
   obj.list,sub.rows,i,j,path_list)

save(s106,s116,
     s017,s029,
     s031,s033,
     s042,s045,
     file = paste0(datadir,"/",date,"_",project,"_seuratobjs.rdata"))

# load(file = paste0(datadir,"/",date,"_",project,"_seuratobjs.rdata"))

  #### Add sample metadata ####
s106[["Age"]] <- "Young"
s116[["Age"]] <- "Young"
s017[["Age"]] <- "Aging"
s029[["Age"]] <- "Aging"
s031[["Age"]] <- "Aging"
s033[["Age"]] <- "Aging"
s042[["Age"]] <- "Aging"
s045[["Age"]] <- "Aging"

s106[["Treatment"]] <- "Veh"
s116[["Treatment"]] <- "Veh"
s017[["Treatment"]] <- "iMPs"
s029[["Treatment"]] <- "Veh"
s031[["Treatment"]] <- "iMPs"
s033[["Treatment"]] <- "Veh"
s042[["Treatment"]] <- "Veh"
s045[["Treatment"]] <- "iMPs"

s106[["Type"]] <- "Young Veh"
s116[["Type"]] <- "Young Veh"
s017[["Type"]] <- "Aging iMPs"
s029[["Type"]] <- "Aging Veh"
s031[["Type"]] <- "Aging iMPs"
s033[["Type"]] <- "Aging Veh"
s042[["Type"]] <- "Aging Veh"
s045[["Type"]] <- "Aging iMPs"

  #### Merge all data into a single seurat object ####

tenx1 <- merge(s106, y = c(s116,
                           s017,s029,
                           s031,s033,
                           s042,s045),
               add.cell.ids = c("YV-106","YV-116",
                                "AI-017","AV-029",
                                "AI-031","AV-033",
                                "AV-042","AI-045"), 
               project = project)

tenx1

rm(s106,s116,
   s017,s029,
   s031,s033,
   s042,s045)

# save merged seurat object 
saveRDS(tenx1, file = paste0("./",date,"_",project,"_merged_seurat_prefilter.rds"))

# tenx1 <- read_rds(file = paste0("./",date,"_",project,"_merged_seurat_prefilter.rds"))

  #### Initial QC Filtering ####

data <- tenx1 # save merged object in a temp "data" object

# Add mito and ribo info to metadata
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^Rp[sl]")

##Plots of UMI#, Gene#, %mito & %ribo
vp1 <- VlnPlot(data, features = "nCount_RNA", pt.size = 0) +
  NoLegend() +
  ggtitle("nCount_RNA prefiltered")
vp1
pdf(paste0("./",date,"_",project,"_prefilter_","nCount_RNA",".pdf"))
print(vp1)
dev.off()

vp2 <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0) +
  NoLegend() +
  ggtitle("nFeature_RNA prefiltered")
vp2
pdf(paste0("./",date,"_",project,"_prefilter_","nFeature_RNA",".pdf"))
print(vp2)
dev.off()

vp3 <- VlnPlot(data, features = "percent.mt", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.mt prefiltered")
vp3
pdf(paste0("./",date,"_",project,"_prefilter_","percent.mt",".pdf"))
print(vp3)
dev.off()

vp4 <- VlnPlot(data, features = "percent.ribo", pt.size = 0) +
  NoLegend() +
  ggtitle("percent.ribo prefiltered")
vp4
pdf(paste0("./",date,"_",project,"_prefilter_","percent.ribo",".pdf"))
print(vp4)
dev.off()

p1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") 
p2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
p1 + p2

pdf(paste0("./",date,"_",project,"_prefilter_","scatterplots",".pdf"))
print(p1 + p2)
dev.off()

# Additional QC and Z-score QC metrics
data[["nUMI.z"]] <- scale(data$nCount_RNA)
data[["nGene.z"]] <- scale(data$nFeature_RNA)
data[["percent.mt.z"]] <- scale(data$percent.mt)
data[["percent.ribo.z"]] <- scale(data$percent.ribo)

##Filter cells based on Z-score
length(data@meta.data$orig.ident)
mean(data@meta.data$percent.mt)
mean(data@meta.data$percent.ribo)
mean(data@meta.data$nCount_RNA)
median(data@meta.data$nCount_RNA)
mean(data@meta.data$nFeature_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

sink(paste0("./",date,"_",project,"_prefilter_QC_metrics.txt"))
cat("length")
length(data@meta.data$orig.ident)
cat("mean percent mito")
mean(data@meta.data$percent.mt)
cat("mean percent ribo")
mean(data@meta.data$percent.ribo)
cat("mean counts")
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

data <- subset(data, subset = percent.mt.z < 3)
data <- subset(data, subset = percent.ribo.z < 3)
data <- subset(data, subset = nGene.z < 3)
data <- subset(data, subset = nUMI.z < 3)

length(data@meta.data$orig.ident)
mean(data@meta.data$percent.mt)
mean(data@meta.data$percent.ribo)
mean(data@meta.data$nCount_RNA)
median(data@meta.data$nCount_RNA)
mean(data@meta.data$nFeature_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

sink(paste0("./",date,"_",project,"_postfilter_QC_metrics.txt"))
cat("length")
length(data@meta.data$orig.ident)
cat("mean percent mito")
mean(data@meta.data$percent.mt)
cat("mean percent ribo")
mean(data@meta.data$percent.ribo)
cat("mean counts")
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

##Plots of UMI#, Gene#, %mito & %ribo
vp5 <- VlnPlot(data, features = "nCount_RNA", pt.size = 0, group.by = "orig.ident")  + 
  NoLegend() + 
  ggtitle("nCount_RNA filtered")
vp5
pdf(paste0("./",date,"_",project,"_filtered_","nCount_RNA",".pdf"))
print(vp5)
dev.off()

vp6 <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") + 
  NoLegend() + 
  ggtitle("nFeature_RNA filtered")
vp6
pdf(paste0("./",date,"_",project,"_filtered_","nFeature_RNA",".pdf"))
print(vp6)
dev.off()

vp7 <- VlnPlot(data, features = "percent.mt", pt.size = 0, group.by = "orig.ident") +
  NoLegend() + 
  ggtitle("percent.mt filtered")
vp7
pdf(paste0("./",date,"_",project,"_filtered_","percent.mt",".pdf"))
print(vp7)
dev.off()

vp8 <- VlnPlot(data, features = "percent.ribo", pt.size = 0, group.by = "orig.ident") +
  NoLegend() + 
  ggtitle("percent.ribo filtered")
vp8
pdf(paste0("./",date,"_",project,"_filtered_","percent.ribo",".pdf"))
print(vp8)
dev.off()

p3 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") 
p4 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
p3 + p4

pdf(paste0("./",date,"_",project,"_filtered_","scatterplots",".pdf"))
print(p3 + p4)
dev.off()

# Additional QC metrics
qc.metrics <- as_tibble(data[[c("nCount_RNA",
                                "nFeature_RNA",
                                "percent.mt",
                                "percent.ribo")]],
                        rownames="Cell.Barcode")

p5 <- qc.metrics %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics")
p5
pdf(paste0("./",date,"_",project,"QC_Metrics",".pdf"))
print(p5)
dev.off()

p6 <- qc.metrics %>%
  ggplot(aes(percent.mt)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Mitochondrion")
p6
pdf(paste0("./",date,"_",project,"Dist_mito",".pdf"))
print(p6)
dev.off()

p7 <- qc.metrics %>%
  ggplot(aes(percent.ribo)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Ribosomal")
p7
pdf(paste0("./",date,"_",project,"Dist_ribo",".pdf"))
print(p7)
dev.off()

rm(vp1,vp2,vp3,vp4,vp5,vp6,vp7,vp8,
   p1,p2,p3,p4,p5,p6,p7,qc.metrics)

tenx1 <- data

rm(data)

  #### set factor levels ####
seurat_sct <- tenx1
seurat_sct$orig.ident <- factor(x = seurat_sct$orig.ident, levels = c("YV-106","YV-116",
                                                                      "AV-029","AV-033","AV-042",
                                                                      "AI-017","AI-031","AI-045"))
seurat_sct$Type <- factor(x = seurat_sct$Type, levels = c("Young Veh","Aging Veh","Aging iMPs"))
seurat_sct$Treatment <- factor(x = seurat_sct$Treatment, levels = c("Veh","iMPs"))
seurat_sct$Age <- factor(x = seurat_sct$Age, levels = c("Young","Aging"))

rm(tenx1)

  #### set color scheme ####

aging_palette <- c("#000080","#666633","#80003F")

saveRDS(seurat_sct, file = paste0("./",date,"_",project,"_seurat_post_qc"))

# seurat_sct <- read_rds(file = paste0("./",date,"_",project,"_seurat_post_qc"))

  #### SCTransform ####
seurat_sct <- SCTransform(seurat_sct, vars.to.regress = c("percent.mt"), verbose = TRUE)

saveRDS(seurat_sct, file = paste0("./",date,"_",project,"_post_sct_pre_pca.rds"))

# seurat_sct <- readRDS(file = paste0("./",date,"_",project,"_post_sct_pre_pca.rds"))

  #### Run PCA ####

# this performs PCA on the seurat object
seurat_sct <- RunPCA(seurat_sct, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(seurat_sct@reductions$pca@cell.embeddings)
# write.table(xx.coord, file = paste0("./",date,"_",project,"_PCAcoord_R.txt"), col.names = NA, sep = "\t", row.names = T)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(seurat_sct@reductions$pca@feature.loadings)
# write.table(xx.gload, file = paste0("./",date,"_",project,"_PCAgload_R.txt"), col.names = NA, sep = "\t", row.names = T)

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

# write csv for eigenvalues
# write.csv(eigenvalues, file = paste0("./",date,"_",project,"_PCA_eigenvalues.csv"), row.names = F)

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

seurat_sct <- RunUMAP(seurat_sct, 
                      reduction = "pca", 
                      dims = 1:meaningful.PCs, 
                      verbose = TRUE)

# update prefixed variable
prefixPC <- paste0("./",date,"_",project,"_",meaningful.PCs,"PCs")

## UMAP plot by sample name ("orig.ident")
DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "orig.ident")

hirestiff(paste0(prefixPC,"_UMAP_by_orig.ident_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_orig.ident_lowres.tiff"))

DimPlot(seurat_sct, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "orig.ident", 
        group.by = "orig.ident", 
        combine = TRUE) + 
  NoLegend()

hirestiff(paste0(prefixPC,"_UMAP_split_by_orig.ident_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_orig.ident_lowres.tiff"))

# These "for loops" are used to generate UMAPs with the specific identity group
# in one UMAP as opposed to the default "split plot" behavior where all identities
# are squished onto one UMAP

Idents(seurat_sct) <- "orig.ident"

for (i in 1:length(levels(seurat_sct))){
  p <- (DimPlot(seurat_sct, 
          reduction = "umap", 
          label = FALSE, 
          pt.size = .25, 
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
        cols = aging_palette)

hirestiff(paste0(prefixPC,"_UMAP_by_type_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_type_lowres.tiff"))

DimPlot(seurat_sct, reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "Type", 
        group.by = "Type",
        cols = aging_palette) + 
  NoLegend()

hirestiff(paste0(prefixPC,"_UMAP_split_by_type_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_split_by_type_lowres.tiff"))

Idents(seurat_sct) <- "Type"

for (i in 1:length(levels(seurat_sct))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .25, 
                cols = aging_palette[i],
                cells = c(WhichCells(seurat_sct, idents = levels(seurat_sct)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct)[i])))
  
  print(p)
  
  hirestiff(paste0(prefixPC,"_UMAP_split_by_type_",levels(seurat_sct)[i],"_hires.tiff"))
  lowrestiff(paste0(prefixPC,"_UMAP_split_by_type_",levels(seurat_sct)[i],"_lowres.tiff"))
}

  #### Add Cell cycle info ####

DefaultAssay(seurat_sct) <- "RNA"

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

DimPlot(seurat_sct, reduction = "umap", label = FALSE, pt.size = .25, group.by = "Phase")
hirestiff(paste0(prefixPC,"_UMAP_by_phase_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_phase_lowres.tiff"))

DimPlot(seurat_sct, reduction = "umap", label = FALSE, pt.size = .25, split.by = "Type", group.by = "Phase")
hirestiff(paste0(prefixPC,"_UMAP_by_phase_split_by_type_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_phase_split_by_type_lowres.tiff"))

DimPlot(seurat_sct, reduction = "umap", label = FALSE, pt.size = .25, split.by = "Phase", group.by = "Phase")
hirestiff(paste0(prefixPC,"_UMAP_by_phase_split_by_phase_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_phase_split_by_phase_lowres.tiff"))

Idents(seurat_sct) <- "Phase"

for (i in 1:length(levels(seurat_sct))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .25, 
                cols = Hue_Pal(length(levels(seurat_sct)))[i],
                cells = c(WhichCells(seurat_sct, idents = levels(seurat_sct)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct)[i])))
  
  print(p)
  hirestiff(paste0(prefixPC,"_UMAP_split_by_phase_",levels(seurat_sct)[i],"_hires.tiff"))
  lowrestiff(paste0(prefixPC,"_UMAP_split_by_phase_",levels(seurat_sct)[i],"_lowres.tiff"))
}

rm(mm.cc.genes)

saveRDS(seurat_sct, file = paste0(prefixPC,"_seurat_post_umap.rds"))

# seurat_sct <- read_rds(file = paste0(prefixPC,"_seurat_post_umap.rds"))

  #### Clustering and Resolution ####

DefaultAssay(seurat_sct) <- "SCT"

# Determine the K-nearest neighbor graph
seurat_sct <- FindNeighbors(object = seurat_sct, reduction = "pca", dims = 1:meaningful.PCs)

# Determine the clusters                              
seurat_sct <- FindClusters(object = seurat_sct,
                           resolution = c(0.3))

res <- "_res.0.3"
ident_res <- paste0("SCT_snn",res)

Idents(seurat_sct) <- ident_res
DimPlot(seurat_sct, reduction = "umap", label = TRUE, label.size = 5) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

#update prefix
prefixPCres <- paste0(prefixPC,res)

# Plot the UMAP
DimPlot(seurat_sct, reduction = "umap", label = TRUE, label.size = 5)  + NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","by","cluster","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","by","cluster","lowres.tiff",sep = "_"))

# UMAP of cells in each cluster by type without cluster labels
DimPlot(seurat_sct, reduction = "umap", label = FALSE, split.by = "Type")  + NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","split_by","type","without_cluster_labels","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","split_by","type","without_cluster_labels","lowres.tiff",sep = "_"))

for (i in 1:length(levels(seurat_sct$Type))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .4, 
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
DimPlot(seurat_sct, reduction = "umap", label = TRUE, split.by = "Type")  + NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste(prefixPCres,"UMAP","split_by","type","with_cluster_labels","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","split_by","type","with_cluster_labels","lowres.tiff",sep = "_"))

for (i in 1:length(levels(seurat_sct$Type))){
  p <- (DimPlot(seurat_sct, 
                reduction = "umap", 
                label = TRUE, 
                pt.size = .25, 
                cells = c(WhichCells(seurat_sct, expression = Type == levels(seurat_sct$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","split_by","type","with_cluster_labels",
                  levels(seurat_sct$Type)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","split_by","type","with_cluster_labels",
                   levels(seurat_sct$Type)[i],"lowres.tiff",sep = "_"))
}

  #### Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(seurat_sct, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
n_cells
write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_cluster.csv"))

rm(n_cells)

  #### Save RDS containing reduction and cluster idents ####
saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_clustering.rds"))

# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_clustering.rds"))

  #### Normalize RNA, find variable features, & scale data for visualization ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(seurat_sct) <- "RNA"

# Normalize, find variable features, scale data 
seurat_sct <- NormalizeData(seurat_sct, verbose = FALSE)
seurat_sct <- FindVariableFeatures(seurat_sct)
all.genes <- rownames(seurat_sct)
seurat_sct <- ScaleData(seurat_sct, features = all.genes)

# save RDS containing RNA normalized data
saveRDS(seurat_sct, paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

# seurat_sct <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

  #### Cell types ####

# 0 - Excitatory Neurons EXCN
VlnPlot(seurat_sct,features = c("Slc17a7","Neurod6","Camk2a","Nell2"),stack = TRUE,flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 0,2,4 - Excitatory Neurons")
hirestiff(paste(prefixPCres,"vln","Cluster","0_2_4","ExNeurons","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","0_2_4","ExNeurons","lowres.tiff", sep = "_"))

# 1 - Oligodendrocytes OLIG
VlnPlot(seurat_sct,features = c("Cnp","Cldn11","Mbp","Mog","Olig1","Sox10"),stack = TRUE,flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 1 - Oligodendrocytes")
hirestiff(paste(prefixPCres,"vln","Cluster","1","Oligodendrocytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","1","Oligodendrocytes","lowres.tiff", sep = "_"))

# 2 - EXCN

# 3 - DG Granule (which also express glut markers) GRN
VlnPlot(seurat_sct,
        features = c("Bdnf","Calb1","Calm1","Calm2","Calm3",
                     "Ppp3ca","Disc1","Gabrb1","Gabra4","Gria1",
                     "Gria2","Gria3","Grin1","Grin2b","Chrm1",
                     "Nrgn","Nrp2","Pcp4","Ncam1","Prox1",
                     "Vsnl1","Npy1r","Actn2","Grm1","Grm2",
                     "Grm3","Grm5","Creb1","Slc17a7"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 3 - DG Granule Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","3","DG_Granule","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","3","DG_Granule","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Bdnf","Calb1","Disc1","Gabrb1","Gabra4",
                     "Gria1","Gria2","Gria3","Grin1","Grin2b",
                     "Chrm1","Nrp2","Ncam1","Prox1","Npy1r",
                     "Actn2","Grm1","Grm2","Grm3","Grm5",
                     "Creb1","Slc17a7"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 3 - DG Granule Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","3","DG_Granule","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","3","DG_Granule","lowres.tiff", sep = "_"))

# 4 - EXCN

# 5 - CA3 Pyramidal PYR-CA3
VlnPlot(seurat_sct,
        features = c("Bdnf","Bok","Calm1","Calm2","Calm3",
                     "Ppp3ca","Gabrb1","Gabra4","Gabra5",
                     "Gria1","Man1a","Neurod6","Chrm1","Chrm3",
                     "Nrgn","Nrp2","Pou3f1","Vsnl1","Npy2r",
                     "Actn2","Grm1","Grm5","Creb1","Slc17a7"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 5 - CA3 Pyramidal Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","5","CA3","pyramidal","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","5","CA3","pyramidal","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Bdnf","Bok","Slc1a1","Gabrb1","Gabra4",
                     "Gabra5","Gria1","Man1a","Neurod6","Chrm1","Chrm3",
                     "Satb2","Nrp2","Pou3f1","Npy2r","Actn2",
                     "Grm1","Grm2","Grm3","Grm5","Creb1",
                     "Slc17a7"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 5 - CA3 Pyramidal Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","5","CA3","pyramidal","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","5","CA3","pyramidal","lowres.tiff", sep = "_"))

# 6 - CA1 Pyramidal PYR-CA1
VlnPlot(seurat_sct,
        features = c("Bdnf","Calm1","Calm2","Calm3","Ppp3ca",
                     "Slc1a1","Gabra4","Man1a","Neurod6",
                     "Chrm1","Chrm3","Nrgn","Satb2",
                     "Pou3f1","Vsnl1","Actn2","Grm1","Grm2",
                     "Grm3","Grm5","Creb1","Slc17a7"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 6 - CA1 Pyramidal Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","6","CA1","pyramidal","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","6","CA1","pyramidal","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Bdnf","Bok","Slc1a1","Gabrb1","Gabra4",
                     "Gabra5","Gria1","Man1a","Neurod6","Chrm1","Chrm3",
                     "Satb2","Nrp2","Pou3f1","Npy2r","Actn2",
                     "Grm1","Grm2","Grm3","Grm5","Creb1",
                     "Slc17a7"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 6 - CA1 Pyramidal Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","6","CA1","pyramidal","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","6","CA1","pyramidal","lowres.tiff", sep = "_"))

# 7 - Inhibitory Neurons INHN
VlnPlot(seurat_sct,
        features = c("Gad1","Synpr","Gad2","Sst"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 7 - Inhibitory Neurons")
hirestiff(paste(prefixPCres,"vln","Cluster","7","Inhibitory","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","7","Inhibitory","lowres.tiff", sep = "_"))

# 8 - Astrocytes ASTR
VlnPlot(seurat_sct,
        features = c("Rgs20","Aqp4","Ntsr2","Gfap","S100b"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 8 - Astrocytes")
hirestiff(paste(prefixPCres,"vln","Cluster","8","Astrocytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","8","Astrocytes","lowres.tiff", sep = "_"))

# 9 - Microglia MG
VlnPlot(seurat_sct,
        features = c("Itgam","Ptprc","Tmem119","Cx3cr1","P2ry12","C1qa"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 9 - Microglia")
hirestiff(paste(prefixPCres,"vln","Cluster","9","Microglia","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","9","Microglia","lowres.tiff", sep = "_"))

# 10 - DG Basket (inhibitory via Zhong) BSKT
VlnPlot(seurat_sct,
        features = c("Dlx6os1","Maf","Erbb4","Gad2","Grik1",
                     "Nxph1","Slc6a1","Gad1","Gm13629","Dlx1",
                     "Dlx6","Col19a1","Utrn","Reln","Npas1",
                     "Sox6","Crhbp","Npy","Sst",
                     "Gabra1","Pvalb"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 10 - DG Basket Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","10","DG","Basket","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","10","DG","Basket","lowres.tiff", sep = "_"))

# 11 - DG Mossy MOSS
VlnPlot(seurat_sct,
        features = c("Gria2","Gria3","Pcp4",
                     "Gabra6","Pvalb","Nos1"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 11 - DG Mossy Cells (neg for Gabra6, Pvalb & Nos1)")
hirestiff(paste(prefixPCres,"vln","Cluster","11","DG","Mossy","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","11","DG","Mossy","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Gria2","Gria3","Pcp4"),
        stack = TRUE,
        flip = TRUE,
        idents = 11,
        split.by = "Type",
        cols = aging_palette) + 
  ggtitle("Cluster 11 - DG Mossy Cells split by Type")
lowrestiff(paste(prefixPCres,"vln","Cluster","11","DG","Mossy","by","type","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,features = c("Gria2"),split.by = "Type",idents = c("11"),flip = TRUE)
lowrestiff(paste(prefixPCres,"vln","Cluster","11","DG","Mossy","by","type","Gria2","lowres.tiff", sep = "_"))
VlnPlot(seurat_sct,features = c("Gria3"),split.by = "Type",idents = c("11"),flip = TRUE)
lowrestiff(paste(prefixPCres,"vln","Cluster","11","DG","Mossy","by","type","Gria3","lowres.tiff", sep = "_"))
VlnPlot(seurat_sct,features = c("Pcp4"),split.by = "Type",idents = c("11"),flip = TRUE)
lowrestiff(paste(prefixPCres,"vln","Cluster","11","DG","Mossy","by","type","Pcp4","lowres.tiff", sep = "_"))

# 12 - OPC
VlnPlot(seurat_sct,
        features = c("Pdgfra","Cspg4","Sox10","Olig1","Olig2"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 12 - OPCs")
hirestiff(paste(prefixPCres,"vln","Cluster","12","OPCs","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","12","OPCs","lowres.tiff", sep = "_"))

# 13 - CPC Choroid Plexus
VlnPlot(seurat_sct,features = c("Ttr","Atp5g1","Clic6","Ndufv3","Uqcrb",
                                "Hemk1","Atp5l","Ndufa4","Ppp1r1b","Cox7b",
                                "Atp5mpl","Atp5md","Ndufa1","Uqcrq","Cox6c",
                                "Cox7c","Atp5j","Atp5e","Uqcr10","Cox8b",
                                "Acad8","Uqcrh","Atp5j2","Uqcr11","Cox5b",
                                "Cox8a","Cox6b1","Ldhb"),
        stack = TRUE,flip = TRUE) + NoLegend() + ggtitle("Cluster 13 - CPCs")
hirestiff(paste(prefixPCres,"vln","Cluster","13","CPC","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","13","CPC","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("Ttr","Folr1","Prlr","Kl"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 13 - Choroid Plexus Cells")
hirestiff(paste(prefixPCres,"vln","Cluster","13","CPC","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","13","CPC","lowres.tiff", sep = "_"))

# 14 - Vascular Leptomeningeal Cells VLMC
VlnPlot(seurat_sct,
        features = c("Pdgfra","Gja1","Slc6a13"),
        stack = TRUE,flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 14 - VLMC")
hirestiff(paste(prefixPCres,"vln","Cluster","14","VLMC","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","14","VLMC","lowres.tiff", sep = "_"))

FeaturePlot(seurat_sct, 
            reduction = "umap", 
            features = c("Pdgfra","Gja1","Slc6a13","Nes"), 
            order = TRUE,
            # min.cutoff = 'q10', 
            label = TRUE,
            label.size = 3,
            pt.size = 0.3)

# 15 - - Activated Microglia MG-A
VlnPlot(seurat_sct,features = c("Itgam","Tmem119","Cx3cr1",
                                "P2ry12","Mrc1",
                                "Jun","Junb"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 15 - Activated Microglia")
hirestiff(paste(prefixPCres,"vln","Cluster","15","Activated","Microglia","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","15","Activated","Microglia","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,features = c("Tlr4","Fosb","Jun","Junb","Jund"),
        # idents = c(9,15),
        # split.by = "Type",
        stack = TRUE,
        flip = TRUE)

VlnPlot(seurat_sct,features = c("Mrc1"),split.by = "Type",idents = c("10","16"),flip = TRUE)
lowrestiff(paste(prefixPCres,"vln","Cluster","15","Activated","Microglia","Mrc1","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,features = c("Tlr4"),split.by = "Type",idents = c("10","16"),flip = TRUE)
lowrestiff(paste(prefixPCres,"vln","Cluster","15","Activated","Microglia","Tlr4","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,features = c("Cd163"),split.by = "Type",idents = c("10","16"),flip = TRUE)
lowrestiff(paste(prefixPCres,"vln","Cluster","15","Activated","Microglia","Cd163","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,features = c("Tmem119"),split.by = "Type",idents = c("9","15"),flip = TRUE)
VlnPlot(seurat_sct,features = c("Jund"),split.by = "Type",idents = c("9","15"),flip = TRUE)
VlnPlot(seurat_sct,features = c("Mapk14"),split.by = "Type",idents = c("9","15"),flip = TRUE)

# 16 -  Pericytes PERI
VlnPlot(seurat_sct,
        features = c("Kcnj8","Cspg4","Pdgfrb","Rgs5"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 16 - Pericytes")
hirestiff(paste(prefixPCres,"vln","Cluster","16","Pericytes","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","16","Pericytes","lowres.tiff", sep = "_"))

# 19 - Multiple
VlnPlot(seurat_sct,features = c("Creb1","Lrrtm4","Lingo2","Tmem108","Actn2"),
        stack = TRUE,flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 19 - GRNL cells")

hirestiff(paste(prefixPCres,"vln","Cluster","19","CPGLcells","res.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","19","CPGLcells","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,features = c("Tubb1","Ncam1","Olig2","Pdgfra","Cspg4",
                                "Mbp","Cnp","Penk","Pcp4","Cck"),
        stack = TRUE,flip = TRUE) + 
  NoLegend() + ggtitle("Cluster 19 - CPGRNL cells")

# 17 - multiple

# 18 - A1 Astrocytes ASTR-A1
VlnPlot(seurat_sct,
        features = c("Rgs20","Aqp4","Ntsr2","Gfap","S100b","C3"),
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() + 
  ggtitle("Cluster 18 - A1 Astrocytes")

hirestiff(paste(prefixPCres,"vln","Cluster","18","A1","Astrocytes","res.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","Cluster","18","A1","Astrocytes","lowres.tiff", sep = "_"))

VlnPlot(seurat_sct,
        features = c("C3"),
        split.by = "Type",
        idents = c("8","18"),
        flip = TRUE,
        cols = aging_palette)

hirestiff(paste(prefixPCres,"vln","astro","A1","C3","by","Type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"vln","astro","A1","C3","by","Type","lowres.tiff", sep = "_"))

# 19 - multiple


  #### consolidate clusters into cell types####

# remove clusters with more than one cell type marker (multiplets)
seurat_sct <- subset(seurat_sct, idents = c(17,19), invert = TRUE)

DimPlot(seurat_sct,
        reduction = "umap",
        label = TRUE,
        label.size = 4,
        pt.size = 0.5) + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","confounding","clusters","removed","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","confounding","clusters","removed","lowres.tiff", sep = "_"))

## new cell counts after filtering

FetchData(seurat_sct, vars = c("Type")) %>%
  dplyr::count(Type) 

# save a new object and assign cell types to cluster numbers
seurat_sct_gen_clusters <- seurat_sct #make new object in case need to revert to numbers

seurat_sct_gen_clusters <- RenameIdents(object = seurat_sct_gen_clusters, 
                                        '0' = "EXCN",
                                        '1' = "OLIG", 
                                        '2' = "EXCN", 
                                        '3' = "GRN",
                                        '4' = "EXCN",
                                        '5' = "PYR-CA3",
                                        '6' = "PYR-CA1",
                                        '7' = "INHN",
                                        '8' = "ASTR",
                                        '9' = "MG",
                                        '10' = "BSKT",
                                        '11' = "MOSS",
                                        '12' = "OPC",
                                        '13' = "CPC",
                                        '14' = "VLMC",
                                        '15' = "MG-A",
                                        '16' = "PERI",
                                        '18' = "ASTR-A1")

# set the order of cell types
seurat_sct_gen_clusters@active.ident <- factor(x = seurat_sct_gen_clusters@active.ident, 
                                               levels = c("INHN","BSKT","MOSS",
                                                          "EXCN","PYR-CA1","PYR-CA3","GRN",
                                                          "MG","MG-A",
                                                          "ASTR","ASTR-A1",
                                                          "OLIG","OPC",
                                                          "VLMC","CPC","PERI"))

# set cell type color
cell_cols <- c("INHN"="#00BFC4","BSKT"="#B79F00","MOSS"="#00BA38",
               "EXCN"="#F8766D","PYR-CA1"="#619CFF","PYR-CA3"="#F564E3","GRN"="#E88526",
               "MG"="#93AA00","MG-A"="#00BF74",
               "ASTR"="#00B9E3","ASTR-A1"="#F17D50",
               "OLIG"="#8E92FF","OPC"="#D39200",
               "VLMC"="#5EB300","CPC"="#00C19F","PERI"="#00ADFA")

# UMAPS with cell type labels
DimPlot(seurat_sct_gen_clusters,
        reduction = "umap",
        label = TRUE,
        label.size = 8,
        pt.size = 0.5,
        cols = cell_cols) + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","cell","types","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","cell","types","lowres.tiff", sep = "_"))
hirestiffsquare(paste(prefixPCres,"UMAP","cell","types","hires_square.tiff", sep = "_"))

DimPlot(seurat_sct_gen_clusters,
        reduction = "umap",
        label = FALSE,
        label.size = 4,
        pt.size = 0.5,
        cols = cell_cols) + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","cell","types","no","label","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","cell","types","no","label","lowres.tiff", sep = "_"))

DimPlot(seurat_sct_gen_clusters, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 3, 
        pt.size = 0.3,
        split.by = "Type",
        cols = cell_cols) + 
  NoLegend()
hirestiff(paste(prefixPCres,"UMAP","cell","types","by","type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"UMAP","cell","types","by","type","lowres.tiff", sep = "_"))

for (i in 1:length(levels(seurat_sct_gen_clusters$Type))){
  p <- (DimPlot(seurat_sct_gen_clusters, 
                reduction = "umap", 
                label = TRUE,
                label.size = 8,
                pt.size = 0.5,
                cols = cell_cols, 
                cells = c(WhichCells(seurat_sct_gen_clusters, 
                                     expression = Type == levels(seurat_sct_gen_clusters$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct_gen_clusters$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                  levels(seurat_sct_gen_clusters$Type)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                   levels(seurat_sct_gen_clusters$Type)[i],"lowres.tiff",sep = "_"))
  hirestiffsquare(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                        levels(seurat_sct_gen_clusters$Type)[i],"hires_square.tiff",sep = "_"))
}

for (i in 1:length(levels(seurat_sct_gen_clusters$Type))){
  p <- (DimPlot(seurat_sct_gen_clusters, 
                reduction = "umap", 
                label = FALSE,
                label.size = 4,
                pt.size = 0.5,
                cols = cell_cols, 
                cells = c(WhichCells(seurat_sct_gen_clusters, 
                                     expression = Type == levels(seurat_sct_gen_clusters$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct_gen_clusters$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                  levels(seurat_sct_gen_clusters$Type)[i],"no","label","hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","cell","types","split_by","type",
                   levels(seurat_sct_gen_clusters$Type)[i],"no","label","lowres.tiff",sep = "_"))
}

  #### save cell type as a metadata column ####
seurat_sct_gen_clusters[["Cell.Type"]] <- seurat_sct_gen_clusters@active.ident

  #### Save cell type labeled seurat obj ####
saveRDS(seurat_sct_gen_clusters, paste0(prefixPCres,"_seurat_cell_types_labeled.rds"))

# seurat_sct_gen_clusters <- readRDS(file = paste0(prefixPCres,"_seurat_cell_types_labeled.rds"))

#### ####
  #### sub out Mossy (MOSS) cluster ####
MOSS_cluster <- subset(seurat_sct_gen_clusters, idents = c("MOSS"))

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
MOSS_cluster <- SCTransform(MOSS_cluster, vars.to.regress = c("percent.mt"), verbose = TRUE)

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

# save scree plot
dev.copy(pdf, paste0(subprefix_MOSS,"_scree_plot.pdf"))
dev.off()

# Number of PCs to use:
# PCs_MOSS <- 13

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

saveRDS(MOSS_cluster, file = paste0(MOSS_prefixpc,"_afterUMAP.rds"))

# MOSS_cluster <- read_rds(file = paste0(MOSS_prefixpc,"_afterUMAP.rds"))

  #### Clustering and Resolution ####

# Determine the K-nearest neighbor graph
MOSS_cluster <- FindNeighbors(object = MOSS_cluster, reduction = "pca",  dims = 1:PCs_MOSS)

# Determine the clusters                              
MOSS_cluster <- FindClusters(object = MOSS_cluster,
                             resolution = c(0.2))

subres_MOSS <- "_res.0.2"
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
        label.size = 8,
        pt.size = 0.75,
        split.by = "Type",
        repel = TRUE) + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MOSS," split by Type"))

hirestiff(paste(MOSS_prefixpcres,"UMAP","by","cluster","split","by","type","hires.tiff",sep = "_"))
lowrestiff(paste(MOSS_prefixpcres,"UMAP","by","cluster","split","by","type","lowres.tiff",sep = "_"))

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

saveRDS(MOSS_cluster, file = paste0(MOSS_prefixpcres,"_after_reclustering.rds"))

# MOSS_cluster <- readRDS(file = paste0(MOSS_prefixpcres,"_after_reclustering.rds"))

  #### Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(MOSS_cluster, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

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

  #### save MOSS objects ####

save(subprefix_MOSS,MOSS_cluster,MOSS_prefix,MOSS_prefixpc,
     MOSS_prefixpcres,PCs_MOSS,subres_MOSS,ident_subres_MOSS,
     file = paste0(MOSS_prefixpcres,"_MOSS_objects.rdata"))

MOSS_prefixpcres <- paste0("./",date,"_",project,"_","MOSS_cluster","_","13","PCs","_res.0.2")

# load(paste0(MOSS_prefixpcres,"_MOSS_objects.rdata"))

#### ####
  #### sub out MG cluster ####

MG_cluster <- subset(seurat_sct_gen_clusters, idents = c("MG"))

DimPlot(MG_cluster, 
        reduction = "umap", 
        label = FALSE, 
        group.by = "Type",
        cols = aging_palette) + 
  ggtitle(paste0(ident_res))

hirestiff(paste(prefixPCres,"MG","cluster","subset","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"MG","cluster","subset","lowres.tiff", sep = "_"))

DimPlot(MG_cluster, 
        reduction = "umap", 
        label = FALSE, 
        split.by = "Type", 
        group.by = "Type", 
        cols = aging_palette) + 
  NoLegend()

hirestiff(paste(prefixPCres,"MG","cluster","subset","by","type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"MG","cluster","subset","by","type","lowres.tiff", sep = "_"))

for (i in 1:length(levels(MG_cluster$Type))){
  p <- (DimPlot(MG_cluster, 
                reduction = "umap", 
                label = FALSE, 
                cols = aging_palette[i],
                cells = c(WhichCells(MG_cluster, expression = Type == levels(MG_cluster$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(MG_cluster$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","MG_cluster",levels(MG_cluster$Type)[i],
                  "split_by","type","hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","MG_cluster",levels(MG_cluster$Type)[i],
                   "split_by","type","lowres.tiff",sep = "_"))
}

  #### SCtransform the subset of interest ####

subprefix_MG <- "MG_cluster"

MG_prefix <- paste0("./",date,"_",project,"_",subprefix_MG)

# Use SCTranform for the next steps
MG_cluster <- SCTransform(MG_cluster, vars.to.regress = c("percent.mt"), verbose = TRUE)

  #### PCA on the seurat object ####
MG_cluster <- RunPCA(MG_cluster, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(MG_cluster@reductions$pca@cell.embeddings)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(MG_cluster@reductions$pca@feature.loadings)

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
PCs_MG <- sum(eig > expected.contribution)

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
     main = paste0(MG_prefix,"scree plot PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = PCs_MG, col = "blue")
text(PCs_MG, cut.off, label = paste("cutoff PC",PCs_MG),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0(subprefix_MG,"_scree_plot.pdf"))
dev.off()

# Number of PCs to use:
# PCs_MG <- 9

MG_prefixpc <- paste0(MG_prefix,"_",PCs_MG,"PCs")

  #### Run UMAP ####
MG_cluster <- RunUMAP(MG_cluster, dims = 1:PCs_MG, reduction = "pca")

# Plot UMAP                             
## UMAP plot by sample name ("orig.ident")
DimPlot(MG_cluster, reduction = "umap", label = FALSE, pt.size = 0.75, group.by = "orig.ident")
hirestiff(paste0(MG_prefixpc,"_UMAP_by_sample_hires.tiff"))
lowrestiff(paste0(MG_prefixpc,"_UMAP_by_sample_lowres.tiff"))

DimPlot(MG_cluster, reduction = "umap", 
        label = FALSE, 
        pt.size = 1.25, 
        group.by = "Type", 
        split.by = "Type", 
        cols = aging_palette) + 
  NoLegend()
hirestiff(paste0(MG_prefixpc,"_UMAP_by_type_hires.tiff"))
lowrestiff(paste0(MG_prefixpc,"_UMAP_by_type_lowres.tiff"))

for (i in 1:length(levels(MG_cluster$Type))){
  p <- (DimPlot(MG_cluster, 
                reduction = "umap", 
                label = FALSE, 
                cols = aging_palette[i],
                cells = c(WhichCells(MG_cluster, expression = Type == levels(MG_cluster$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(MG_cluster$Type)[i])))
  
  print(p)
  
  hirestiff(paste(MG_prefixpc,"UMAP","MGA_cluster",levels(MG_cluster$Type)[i],
                  "split_by","type","hires.tiff",sep = "_"))
  lowrestiff(paste(MG_prefixpc,"UMAP","MGA_cluster",levels(MG_cluster$Type)[i],
                   "split_by","type","lowres.tiff",sep = "_"))
}

saveRDS(MG_cluster, file = paste0(MG_prefixpc,"_afterUMAP.rds"))

# MG_cluster <- read_rds(file = paste0(MG_prefixpc,"_afterUMAP.rds"))

  #### Clustering and Resolution ####

# Determine the K-nearest neighbor graph
MG_cluster <- FindNeighbors(object = MG_cluster, reduction = "pca",  dims = 1:PCs_MG)

# Determine the clusters
MG_cluster <- FindClusters(object = MG_cluster,
                           resolution = c(0.8))

subres_MG <- "_res.0.8"
ident_subres_MG <- paste0("SCT_snn",subres_MG)

Idents(MG_cluster) <- ident_subres_MG
DimPlot(MG_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        pt.size = 2) + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MG))

MG_prefixpcres <- paste0(MG_prefixpc,subres_MG)

# set cluster colors
mg_clust_cols <- scales::hue_pal()(length(levels(MG_cluster)))
names(mg_clust_cols) <- levels(MG_cluster)

DimPlot(MG_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        pt.size = 2) + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MG))

hirestiff(paste(MG_prefixpcres,"UMAP","by","cluster","hires.tiff",sep = "_"))
lowrestiff(paste(MG_prefixpcres,"UMAP","by","cluster","lowres.tiff",sep = "_"))

DimPlot(MG_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 4,
        pt.size = 2,
        split.by = "Type") + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MG," split by Type"))

hirestiff(paste(MG_prefixpcres,"UMAP","by","cluster","split","by","type","hires.tiff",sep = "_"))
lowrestiff(paste(MG_prefixpcres,"UMAP","by","cluster","split","by","type","lowres.tiff",sep = "_"))

for (i in 1:length(levels(MG_cluster$Type))){
  p <- (DimPlot(MG_cluster, 
                reduction = "umap", 
                label = TRUE,
                label.size = 5, 
                pt.size = 2,
                cols = mg_clust_cols, 
                cells = c(WhichCells(MG_cluster, expression = Type == levels(MG_cluster$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(MG_cluster$Type)[i]," in MG cluster")))
  
  print(p)
  
  hirestiff(paste(MG_prefixpcres,"UMAP","split_by","type",
                  levels(MG_cluster$Type)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(MG_prefixpcres,"UMAP","split_by","type",
                   levels(MG_cluster$Type)[i],"lowres.tiff",sep = "_"))
}

saveRDS(MG_cluster, file = paste0(MG_prefixpcres,"_after_reclustering.rds"))

# MG_cluster <- readRDS(file = paste0(MG_prefixpcres,"_after_reclustering.rds"))

  #### Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(MG_cluster, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

write.csv(n_cells, file = paste0(MG_prefixpcres,"_cells_per_cluster.csv"))


  #### Normalize RNA ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(MG_cluster) <- "RNA"

# Normalize, find variable features, scale data 
MG_cluster <- NormalizeData(MG_cluster, verbose = FALSE)
MG_cluster <- FindVariableFeatures(MG_cluster)
all.genes <- rownames(MG_cluster)
MG_cluster <- ScaleData(MG_cluster, features = all.genes)

saveRDS(MG_cluster, file = paste0(MG_prefixpcres,"_after_rna_norm.rds"))

# MG_cluster <- readRDS(file = paste0(MG_prefixpcres,"_after_rna_norm.rds"))

  #### save MG objects ####
save(subprefix_MG,MG_cluster,MG_prefix,MG_prefixpc,
     MG_prefixpcres,PCs_MG,subres_MG,ident_subres_MG,
     file = paste0(MG_prefixpcres,"_MG_objects.rdata"))

MG_prefixpcres <- paste0("./",date,"_",project,"_","MG_cluster","_","9","PCs","_res.0.8")

# load(paste0(MG_prefixpcres,"_MG_objects.rdata"))

#### ####
  #### sub out Activated MG (MGA)cluster ####

MGA_cluster <- subset(seurat_sct_gen_clusters, idents = c("MG-A"))

DimPlot(MGA_cluster, 
        reduction = "umap", 
        label = FALSE, 
        group.by = "Type",
        cols = aging_palette) + 
  ggtitle(paste0(ident_res))

hirestiff(paste(prefixPCres,"MG-A","cluster","subset","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"MG-A","cluster","subset","lowres.tiff", sep = "_"))

# as there is one cell way outside of the others, filter to remove it
MGA_cluster <- subset(MGA_cluster, subset = UMAP_1 < -2)

DimPlot(MGA_cluster, 
        reduction = "umap", 
        label = FALSE, 
        group.by = "Type",
        cols = aging_palette) + 
  ggtitle(paste0(ident_res))

hirestiff(paste(prefixPCres,"MG-A","cluster","subset","filtered","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"MG-A","cluster","subset","filtered","lowres.tiff", sep = "_"))

DimPlot(MGA_cluster, 
        reduction = "umap", 
        label = FALSE, 
        split.by = "Type", 
        group.by = "Type", 
        cols = aging_palette) + 
  NoLegend()

hirestiff(paste(prefixPCres,"MG-A","cluster","subset","by","type","hires.tiff", sep = "_"))
lowrestiff(paste(prefixPCres,"MG-A","cluster","subset","by","type","lowres.tiff", sep = "_"))

for (i in 1:length(levels(MGA_cluster$Type))){
  p <- (DimPlot(MGA_cluster, 
                reduction = "umap", 
                label = FALSE, 
                cols = aging_palette[i],
                cells = c(WhichCells(MGA_cluster, expression = Type == levels(MGA_cluster$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(MGA_cluster$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","MGA_cluster",levels(MGA_cluster$Type)[i],
                  "split_by","type","hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","MGA_cluster",levels(MGA_cluster$Type)[i],
                   "split_by","type","lowres.tiff",sep = "_"))
}

  #### SCtransform the subset of interest ####
subprefix_MGA <- "MGA_cluster"

MGA_prefix <- paste0("./",date,"_",project,"_",subprefix_MGA)

# Use SCTranform for the next steps
MGA_cluster <- SCTransform(MGA_cluster, vars.to.regress = c("percent.mt"), verbose = TRUE)

  #### PCA on the seurat object ####
MGA_cluster <- RunPCA(MGA_cluster, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(MGA_cluster@reductions$pca@cell.embeddings)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(MGA_cluster@reductions$pca@feature.loadings)

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
PCs_MGA <- sum(eig > expected.contribution)

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
     main = paste0(MGA_prefix,"scree plot PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = PCs_MGA, col = "blue")
text(PCs_MGA, cut.off, label = paste("cutoff PC",PCs_MGA),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0(subprefix_MGA,"_scree_plot.pdf"))
dev.off()

# Number of PCs to use:
# PCs_MGA <- 12

MGA_prefixpc <- paste0(MGA_prefix,"_",PCs_MGA,"PCs")

  #### Run UMAP ####
MGA_cluster <- RunUMAP(MGA_cluster, dims = 1:PCs_MGA, reduction = "pca")

# Plot UMAP                             
## UMAP plot by sample name ("orig.ident")
DimPlot(MGA_cluster, reduction = "umap", label = FALSE, pt.size = 0.75, group.by = "orig.ident")
hirestiff(paste0(MGA_prefixpc,"_UMAP_by_sample_hires.tiff"))
lowrestiff(paste0(MGA_prefixpc,"_UMAP_by_sample_lowres.tiff"))

DimPlot(MGA_cluster, reduction = "umap", 
        label = FALSE, 
        pt.size = 1.25, 
        group.by = "Type", 
        split.by = "Type", 
        cols = aging_palette) + 
  NoLegend()
hirestiff(paste0(MGA_prefixpc,"_UMAP_by_type_hires.tiff"))
lowrestiff(paste0(MGA_prefixpc,"_UMAP_by_type_lowres.tiff"))

for (i in 1:length(levels(MGA_cluster$Type))){
  p <- (DimPlot(MGA_cluster, 
                reduction = "umap", 
                label = FALSE, 
                cols = aging_palette[i],
                cells = c(WhichCells(MGA_cluster, expression = Type == levels(MGA_cluster$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(MGA_cluster$Type)[i])))
  
  print(p)
  
  hirestiff(paste(MGA_prefixpc,"UMAP","MGA_cluster",levels(MGA_cluster$Type)[i],
                  "split_by","type","hires.tiff",sep = "_"))
  lowrestiff(paste(MGA_prefixpc,"UMAP","MGA_cluster",levels(MGA_cluster$Type)[i],
                   "split_by","type","lowres.tiff",sep = "_"))
}

saveRDS(MGA_cluster, file = paste0(MGA_prefixpc,"_afterUMAP.rds"))

# MGA_cluster <- read_rds(file = paste0(MGA_prefixpc,"_afterUMAP.rds"))

  #### Clustering and Resolution ####

# Determine the K-nearest neighbor graph
MGA_cluster <- FindNeighbors(object = MGA_cluster, reduction = "pca",  dims = 1:PCs_MGA)

# Determine the clusters                              
MGA_cluster <- FindClusters(object = MGA_cluster,
                            resolution = c(0.2))

subres_MGA <- "_res.0.2"
ident_subres_MGA <- paste0("SCT_snn",subres_MGA)

Idents(MGA_cluster) <- ident_subres_MGA
DimPlot(MGA_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        pt.size = 2) + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MGA))

MGA_prefixpcres <- paste0(MGA_prefixpc,subres_MGA)

DimPlot(MGA_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        pt.size = 2) + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MGA))

hirestiff(paste(MGA_prefixpcres,"UMAP","by","cluster","hires.tiff",sep = "_"))
lowrestiff(paste(MGA_prefixpcres,"UMAP","by","cluster","lowres.tiff",sep = "_"))

DimPlot(MGA_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        pt.size = 1.5,
        split.by = "Type") + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MGA," split by Type"))

hirestiff(paste(MGA_prefixpcres,"UMAP","by","cluster","split","by","type","hires.tiff",sep = "_"))
lowrestiff(paste(MGA_prefixpcres,"UMAP","by","cluster","split","by","type","lowres.tiff",sep = "_"))

for (i in 1:length(levels(MGA_cluster$Type))){
  p <- (DimPlot(MGA_cluster, 
                reduction = "umap", 
                label = TRUE,
                label.size = 5, 
                pt.size = 2, 
                cells = c(WhichCells(MGA_cluster, expression = Type == levels(MGA_cluster$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(MGA_cluster$Type)[i])))
  
  print(p)
  
  hirestiff(paste(MGA_prefixpcres,"UMAP","split_by","type",
                  levels(MGA_cluster$Type)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(MGA_prefixpcres,"UMAP","split_by","type",
                   levels(MGA_cluster$Type)[i],"lowres.tiff",sep = "_"))
}

saveRDS(MGA_cluster, file = paste0(MGA_prefixpcres,"_after_reclustering.rds"))

# MGA_cluster <- readRDS(file = paste0(MGA_prefixpcres,"_after_reclustering.rds"))

  #### Extract number of cells per cluster per orig.ident ####
n_cells <- FetchData(MGA_cluster, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

write.csv(n_cells, file = paste0(MGA_prefixpcres,"_cells_per_cluster.csv"))

  #### Normalize RNA ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(MGA_cluster) <- "RNA"

# Normalize, find variable features, scale data 
MGA_cluster <- NormalizeData(MGA_cluster, verbose = FALSE)
MGA_cluster <- FindVariableFeatures(MGA_cluster)
all.genes <- rownames(MGA_cluster)
MGA_cluster <- ScaleData(MGA_cluster, features = all.genes)

saveRDS(MGA_cluster, file = paste0(MGA_prefixpcres,"_after_rna_norm.rds"))

MGA_cluster <- readRDS(file = paste0(MGA_prefixpcres,"_after_rna_norm.rds"))

  #### save MGA objects ####
save(subprefix_MGA,MGA_cluster,MGA_prefix,MGA_prefixpc,
     MGA_prefixpcres,PCs_MGA,subres_MGA,ident_subres_MGA,
     file = paste0(MGA_prefixpcres,"_MGA_objects.rdata"))

MGA_prefixpcres <- paste0("./",date,"_",project,"_","MGA_cluster","_","12","PCs","_res.0.2")

# load(paste0(MGA_prefixpcres,"_MGA_objects.rdata"))

#### ####
  #### Load proteomics data ####

protdata <- read_csv("moser_et_al_2022_iMP_Plasma_Proteomics_Raw.csv")

# Differentially expressed proteins ("DEPs") for AgedVeh vs Young were determned by:
# 1) filtered by column "imputed FDR" <=0.05
# 2) filtered FDR-filtered DEPs for absolute value of Log2FC >= 0.58 (~1.5 FC)
# 3) this left 75 DEPs with pos L2FC and 34 DEPs with neg L2FC

deps <- read_csv("moser_et_al_2022_iMP_Plasma_Proteomics_deps.csv")

# filter the raw data by the DEPs
filteredprots <- protdata %>% filter(protdata$Gene %in% deps$Gene)

write.csv(filteredprots, file = "moser_et_al_2022_iMP_Plasma_Proteomics_Raw_Filtered_DEPs.csv")

# Excel was then used to calculate ZScores
# ZScores were copied to new csv file
zscores <- read.csv("moser_et_al_2022_iMP_Plasma_Proteomics_ZScore.csv")

# ZScore data was then prepped for heatmap use using "melt"
zdata <- melt(zscores, id.vars = "Gene")

# Data for vocano plots
volc <- read_csv("moser_et_al_2022_iMP_Plasma_Proteomics_For_Volcano.csv")

volc_colors <- c("#707070", "#00ff00", "#00ffff", "#ff00ff")

#### ####
   #### Figure 2a ####
# volcano plot of DEPs Aging Veh vs Young
EnhancedVolcano(volc,
                lab = volc$Gene,
                x = 'log2FC_AVvsY',
                y = 'FDR_AVvsY',
                xlim = c(min(volc$log2FC_AVvsY, na.rm = TRUE) - 0.1, 
                         max(volc$log2FC_AVvsY, na.rm = TRUE) + 0.1),
                ylim = c(0, max(-log10(volc$FDR_AVvsY), na.rm = TRUE) + 2),
                title = "Aging Veh vs Young",
                titleLabSize = 30,
                caption = NULL,
                subtitle = NULL,
                xlab = bquote(bold(~Log["2"] ~ "fold change")),
                ylab = bquote(bold(~"-" ~Log["10"] ~ FDR)),
                axisLabSize = 25,
                col = volc_colors,
                colAlpha = 1/1,
                legendLabels = c("NS", 
                                 expression(Log[2] ~ FC), 
                                 "FDR", 
                                 expression(FDR ~ and ~ log[2] ~ FC)),
                legendPosition = 'bottom',
                legendLabSize = 25,
                legendIconSize = 7,
                drawConnectors = TRUE,
                min.segment.length = 1,
                labSize = 7,
                max.overlaps = Inf,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 7
) 

   #### Figure 2b ####
# heatmap of DEPs by ZScore
ggplot(zdata, aes(x = variable,
                  y = Gene,
                  fill = value)) + 
  geom_tile() +
  theme(axis.text.y = element_text(size = 5),
        axis.text.x = element_text(size = 0)) +
  guides(fill=guide_legend(title="Z-Score")) +
  scale_fill_viridis_c(option = "A", direction = -1) 

# NOTE: Image was modified as follows for better representation of the data
#       1) image was inverted to put the genes in A to Z order instead of default Z to A
#       2) only every 4th gene name was listed to allow for readability

   #### Figure 2f ####
# volcano plot of DEPs Aging iMPs vs Aging Veh
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
                axisLabSize = 25,
                col = volc_colors,
                colAlpha = 1/1,
                legendLabels = c("NS", 
                                 expression(Log[2] ~ FC), 
                                 "FDR", 
                                 expression(FDR ~ and ~ log[2] ~ FC)),
                legendPosition = 'bottom',
                legendLabSize = 25,
                legendIconSize = 7,
                drawConnectors = TRUE,
                min.segment.length = 1,
                labSize = 7,
                max.overlaps = Inf,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 7
) 

   #### Figure 4a ####

DimPlot(seurat_sct_gen_clusters,
        reduction = "umap",
        label = TRUE,
        label.size = 8,
        pt.size = 0.5,
        cols = cell_cols) + 
  NoLegend()
hirestiffsquare(paste(prefixPCres,"Fig4a","UMAP","cell","types","hires_square.tiff", sep = "_"))

   #### Figure 4b ####
VlnPlot(seurat_sct_gen_clusters,features = c("Gad1","Gad2",            ## INHN
                                             "Dlx1","Col19a1",         ## BSKT
                                             "Pcp4","Calb2",           ## MOSS
                                             "Slc17a7","Neurod6",      ## EXCN
                                             "Satb2","Grm3",           ## PYR-C1
                                             "Chrm1","Actn2",          ## PYR-C3
                                             "Calb1","Prox1",          ## GRN
                                             "Tmem119","Itgam",        ## MG
                                             "Mrc1","Jun",             ## MG-A 
                                             "Aqp4","Gfap",            ## ASTR
                                             "S100b","C3",             ## ASTR-A1
                                             "Mbp","Mog",              ## OLIG
                                             "Olig1","Sox10",          ## OPC
                                             "Pdgfra","Slc6a13",       ## VLMC
                                             "Prlr","Kl",              ## CPC
                                             "Kcnj8","Rgs5"            ## PERI
                                             ),
        stack = TRUE,
        flip = TRUE,
        sort = FALSE) + 
  NoLegend() + 
  ggtitle("All Cell Types") +
  theme(axis.text = element_text(size = 18), 
        axis.title = element_text(size = 22),
        axis.text.y = element_text(size = 12))

hirestiffsquare(paste(prefixPCres,"Fig4b","VLN","cell","type","markers","hires_square.tiff", sep = "_"))

   #### Figure 4c ####
DimPlot(MGA_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5,
        pt.size = 1.5,
        split.by = "Type") + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MGA," split by Type"))

hirestiff(paste(MGA_prefixpcres,"Fig4c","UMAP","by","cluster",
                "split","by","type","hires.tiff",sep = "_"))

   #### Figure 4e ####
# data extracted as follows and graphs generated in Prism
# Expression of Select Genes

expression_data <- FetchData(object = seurat_sct_gen_clusters, 
                             cells = WhichCells(seurat_sct_gen_clusters,
                                                idents = "MG-A"),
                             vars = c("orig.ident","Type","Cell.Type",
                                      "Saa3", "Il6","Tnf","Il10"))

write.csv(expression_data, 
          file = paste(prefixPCres,
                       "Fig4e","data","MG-A",
                       "select","gene","expression",
                       "by","type.csv", sep = "_"))


   #### Figure 4f ####
# data extracted as follows and graphs generated in Prism
# Expression of Select Genes

expression_data <- FetchData(object = seurat_sct_gen_clusters, 
                             cells = WhichCells(seurat_sct_gen_clusters,
                                                idents = "MG-A"),
                             vars = c("orig.ident","Type","Cell.Type",
                                      "C3","Cfp","Cfh","C4b",
                                      "C2","C1qc","C1qb","C1qa"))

write.csv(expression_data, 
          file = paste(prefixPCres,
                       "Fig4f","data","MG-A",
                       "select","gene","expression",
                       "by","type.csv", sep = "_"))
   #### Figure 5a ####
for (i in 1:length(levels(seurat_sct_gen_clusters$Type))){
  p <- (DimPlot(seurat_sct_gen_clusters, 
                reduction = "umap", 
                label = TRUE,
                label.size = 8,
                pt.size = 0.5,
                cols = cell_cols, 
                cells = c(WhichCells(seurat_sct_gen_clusters, 
                                     expression = Type == levels(seurat_sct_gen_clusters$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(seurat_sct_gen_clusters$Type)[i])))
  
  print(p)
  
  hirestiffsquare(paste(prefixPCres,"Fig5a","UMAP","cell","types","split_by","type",
                        levels(seurat_sct_gen_clusters$Type)[i],"hires_square.tiff",sep = "_"))
}

   #### Figure 5c ####
DimPlot(MOSS_cluster, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 8,
        pt.size = 0.75,
        split.by = "Type",
        repel = TRUE) + 
  NoLegend() + 
  ggtitle(paste0(ident_subres_MOSS," split by Type"))

hirestiff(paste(MOSS_prefixpcres,"Fig5c","UMAP","by","cluster",
                "split","by","type","hires.tiff",sep = "_"))
   #### Figure 5e ####
# data extracted as follows and graphs generated in Prism
# Expression of Select Genes

expression_data <- FetchData(object = seurat_sct_gen_clusters, 
                             cells = WhichCells(seurat_sct_gen_clusters,
                                                idents = "MOSS"),
                             vars = c("orig.ident","Type","Cell.Type",
                                      "Kcnma1","Syt9","Vav3","Kcnd2",
                                      "Rgs6","Prkd1","Thsd7a","Trpm3"))

write.csv(expression_data, 
          file = paste(prefixPCres,
                       "Fig5e","data","MOSS","select","gene","expression",
                       "by","type.csv", sep = "_"))

   #### Figure 5f ####
# data extracted as follows and graphs generated in Prism
# Expression of Select Genes
expression_data <- FetchData(object = seurat_sct_gen_clusters, 
                             cells = WhichCells(seurat_sct_gen_clusters,
                                                idents = "MOSS"),
                             vars = c("orig.ident","Type","Cell.Type",
                                      "Dcx"))

write.csv(expression_data, 
          file = paste(prefixPCres,
                       "Fig5f","data","MOSS","Dcx","expression",
                       "by","type.csv", sep = "_"))

   #### Supplementary Figure 3a ####
# volcano plot of DEPs Aging iMPs vs Young
EnhancedVolcano(volc,
                lab = volc$Gene,
                x = 'log2FC_AMvsY',
                y = 'FDR_AMvsY',
                xlim = c(min(volc$log2FC_AMvsY, na.rm = TRUE) - 0.1, 
                         max(volc$log2FC_AMvsY, na.rm = TRUE) + 0.1),
                ylim = c(0, max(-log10(volc$FDR_AMvsY), na.rm = TRUE) + 2),
                title = "Aging iMPs vs Young",
                titleLabSize = 30,
                caption = NULL,
                subtitle = NULL,
                xlab = bquote(bold(~Log["2"] ~ "fold change")),
                ylab = bquote(bold(~"-" ~Log["10"] ~ FDR)),
                axisLabSize = 25,
                col = volc_colors,
                colAlpha = 1/1,
                legendLabels = c("NS", 
                                 expression(Log[2] ~ FC), 
                                 "FDR", 
                                 expression(FDR ~ and ~ log[2] ~ FC)),
                legendPosition = 'bottom',
                legendLabSize = 25,
                legendIconSize = 7,
                drawConnectors = TRUE,
                min.segment.length = 1,
                labSize = 7,
                max.overlaps = Inf,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                arrowheads = FALSE,
                pointSize = 7
                ) 

   #### Supplementary Figure 5a ####
# "Microgilia" column (clusters split by Type)
for (i in 1:length(levels(MG_cluster$Type))){
  p <- (DimPlot(MG_cluster, 
                reduction = "umap", 
                label = FALSE,
                label.size = 5, 
                pt.size = 2,
                cols = mg_clust_cols, 
                cells = c(WhichCells(MG_cluster, 
                                     expression = Type == levels(MG_cluster$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(MG_cluster$Type)[i]," in MG cluster")))
  
  print(p)
  
  hirestiff(paste(MG_prefixpcres,"SuppFig5a","Microglia","UMAP","split_by","type",
                  levels(MG_cluster$Type)[i],"hires.tiff",sep = "_"))
}

# "Tmem119" and "P2ry12" columns (expression split by Type)

genelist <- c("Tmem119","P2ry12")

for (k in 1:length(genelist)){ 
  for (i in 1:length(levels(MG_cluster$Type))){
    p <- (FeaturePlot(MG_cluster, 
                      reduction = "umap", 
                      features = genelist[k], 
                      label = FALSE,
                      label.size = 4,
                      pt.size = 0.75,
                      cells = c(WhichCells(MG_cluster, 
                                           expression = Type == levels(MG_cluster$Type)[i])))) +
      ggtitle(paste0(genelist[k]," in ",levels(MG_cluster$Type)[i], " in MG Cluster"))
    
    print(p)
    
    hirestiff(paste(MG_prefixpcres,"SuppFig5a","Columns2_3","UMAP",
                    genelist[k],"in",levels(MG_cluster$Type)[i],
                    "hires.tiff",sep = "_"))
  }
}
   #### Supplementary Figure 5b ####
# ARMS
VlnPlot(MGA_cluster, 
        features = c("Apoe","Cd74","H2-Ab1",
                     "H2-Aa","Ctsb","Ctsd"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("ARMS")

hirestiff(paste(MGA_prefixpcres,"SuppFig5b","vln","ARMS","genes","hires.tiff", sep = "_"))

VlnPlot(MGA_cluster, 
        features = c("Apoe","Cd74","H2-Ab1",
                     "H2-Aa","Ctsb","Ctsd"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "Type",
        cols = aging_palette) +
  ggtitle("ARMS by Type")

hirestiff(paste(MGA_prefixpcres,"SuppFig5b","vln","ARMS","genes","by","type","hires.tiff", sep = "_"))

   #### Supplementary Figure 5c ####
#IRMS
VlnPlot(MGA_cluster, 
        features = c("Ifi204","Ifi27l2a","Irf7",
                     "Isg15","Oasl2","Rnf213",
                     "Rtp4","Usp18","Xaf1"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("IRMS")

hirestiff(paste(MGA_prefixpcres,"SuppFig5c","vln","IRMS","genes","hires.tiff", sep = "_"))

VlnPlot(MGA_cluster, 
        features = c("Ifi204","Ifi27l2a","Irf7",
                     "Isg15","Oasl2","Rnf213",
                     "Rtp4","Usp18","Xaf1"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "Type",
        cols = aging_palette) +
  ggtitle("IRMS by Type")

hirestiff(paste(MGA_prefixpcres,"SuppFig5c","vln","IRMS","genes","by","type","hires.tiff", sep = "_"))

   #### Supplementary Figure 5d ####
#Aging MG
VlnPlot(MGA_cluster, 
        features = c("Ccl4","Il1b","Ifitm3","Rtp4","Irf7",
                     "Isg15","Oasl2","F13a1","Mrc1","Pf4",
                     "Clec12a","Ms4a7"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("Aging MG Genes")

hirestiff(paste(MGA_prefixpcres,"SuppFig5d","vln","Aging_MG","genes","hires.tiff", sep = "_"))

VlnPlot(MGA_cluster, 
        features = c("Ccl4","Il1b","Ifitm3","Rtp4","Irf7",
                     "Isg15","Oasl2","F13a1","Mrc1","Pf4",
                     "Clec12a","Ms4a7"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "Type",
        cols = aging_palette) + 
  ggtitle("Aging MG Genes by Type")

hirestiff(paste(MGA_prefixpcres,"SuppFig5d","vln","Aging_MG","genes","by","type","hires.tiff", sep = "_"))

   #### Supplementary Figure 5e ####
# stage 1 DAM upreg
VlnPlot(MGA_cluster, 
        features = c("Tyrobp","Apoe","B2m","Trem2"), 
        flip = TRUE, 
        stack = TRUE) +
  NoLegend() +
  ggtitle("DAM Stage 1, Upregulated Genes")

hirestiff(paste(MGA_prefixpcres,"SuppFig5e","vln","DAM_stg1_up","genes","hires.tiff", sep = "_"))

VlnPlot(MGA_cluster, 
        features = c("Tyrobp","Apoe","B2m","Trem2"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "Type",
        cols = aging_palette) +
  ggtitle("DAM Stage 1, Upregulated Genes by Type")

hirestiff(paste(MGA_prefixpcres,"SuppFig5e","vln","DAM_stg1_up","genes","by","type","hires.tiff", sep = "_"))

   #### Supplementary Figure 5f ####
# Serum Amyloid Genes and downstream

VlnPlot(MGA_cluster, 
        features = c("Saa1","Saa2","Saa3",
                     "Il6","Tnf","Il10"), 
        flip = TRUE, 
        stack = TRUE,
        split.by = "Type",
        cols = aging_palette) +
  ggtitle("Serum Amyloid Genes by Type")

hirestiff(paste(MGA_prefixpcres,"SuppFig5f","vln","saa","genes","by","type","hires.tiff", sep = "_"))


   #### Supplementary Figure 5g ####
#  Comp pathway genes 

VlnPlot(MGA_cluster,
        features = c("C3","Cfp","Cfh","C4b","C2",
                     "C1qc","C1qb","C1qa"),
        flip = TRUE,
        stack = TRUE,
        split.by = "Type",
        cols = aging_palette) 
hirestiff(paste(MGA_prefixpcres,"SuppFig5g","vln","Complement","pathway","genes","hires.tiff", sep = "_"))


   #### Supplementary Figure 5h ####
#  Comp pathway genes by Type 

VlnPlot(seurat_sct_gen_clusters,
        features = c("C3","Cfp","Cfh","C4b","C2",
                     "C1qc","C1qb","C1qa"),
        idents = c("MG","MG-A","ASTR","ASTR-A1","VLMC"),
        flip = TRUE,
        stack = TRUE,
        split.by = "Type",
        cols = aging_palette) 
hirestiff(paste(prefixPCres,"SuppFig5h","VLN","Complement","pathway","genes","hires.tiff", sep = "_"))

   #### Supplementary Figure 5i ####

expression_data <- FetchData(object = seurat_sct_gen_clusters, 
                             cells = WhichCells(seurat_sct_gen_clusters,
                                                idents = c("ASTR","ASTR-A1","VLMC")),
                             vars = c("orig.ident","Type","Cell.Type",
                                      "C4b","Cfh","C3"))

write.csv(expression_data, 
          file = paste(prefixPCres,
                       "SuppFig5i","data","A_A1_V","select","gene","expression",
                       "by","type.csv", sep = "_"))

   #### Supplementary Figure 6b ####
# data extracted as follows and graphs generated in Prism
# Expression of Select Genes

expression_data <- FetchData(object = seurat_sct_gen_clusters, 
                             cells = WhichCells(seurat_sct_gen_clusters,
                                                idents = "MOSS"),
                             vars = c("orig.ident","Type","Cell.Type",
                                      "Nwd2","Scube1","Kctd8","Tac2",
                                      "D130079A08Rik","Cadps2","Cpne4","Gfra1",
                                      "Etv1","Slit2","D130009I18Rik","Necab3")
                             )

write.csv(expression_data, 
          file = paste(prefixPCres,
                       "SuppFig6b","data","MOSS","select","gene","expression",
                       "by","type.csv", sep = "_")
          )

  
