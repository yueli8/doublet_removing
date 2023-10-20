
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(devtools)
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(cowplot)
library(reshape2)
library(tidyverse)
library(garnett)
library(data.table)
library(DoubletFinder)
library(harmony)

setwd("D:\\MM\\scRNA.MM\\GSE223972_RAW\\Rdata\\")

load("G_CSF_H1.rdata")
load("G_CSF_H2.rdata")
load("G_CSF_M1.rdata")
load("G_CSF_M2.rdata")
load("G_CSF_M3.rdata")
load("G_CSF_M4.rdata")
load("M-H1.rdata")
load("M-H2.rdata")
load("M_G_CSF_M1.rdata")
load("M_G_CSF_M2.rdata")
load("M_G_CSF_M3.rdata")
load("M_G_CSF_M4.rdata")
load("P_H1.rdata")
load("P_H2.rdata")
load("P_G_CSF_M1.rdata")
load("P_G_CSF_M2.rdata")
load("P_G_CSF_M3.rdata")
load("P_G_CSF_M4.rdata")

##计算质控指标
#计算细胞中线粒体基因比例
G_CSF_H1[["percent.mt"]] <- PercentageFeatureSet(G_CSF_H1, pattern = "^MT-")
G_CSF_H2[["percent.mt"]] <- PercentageFeatureSet(G_CSF_H2, pattern = "^MT-")
G_CSF_M1[["percent.mt"]] <- PercentageFeatureSet(G_CSF_M1, pattern = "^MT-")
G_CSF_M2[["percent.mt"]] <- PercentageFeatureSet(G_CSF_M2, pattern = "^MT-")
G_CSF_M3[["percent.mt"]] <- PercentageFeatureSet(G_CSF_M3, pattern = "^MT-")
G_CSF_M4[["percent.mt"]] <- PercentageFeatureSet(G_CSF_M4, pattern = "^MT-")
M_H1[["percent.mt"]] <- PercentageFeatureSet(M_H1, pattern = "^MT-")
M_H2[["percent.mt"]] <- PercentageFeatureSet(M_H2, pattern = "^MT-")
M_G_CSF_M1[["percent.mt"]] <- PercentageFeatureSet(M_G_CSF_M1, pattern = "^MT-")
M_G_CSF_M2[["percent.mt"]] <- PercentageFeatureSet(M_G_CSF_M2, pattern = "^MT-")
M_G_CSF_M3[["percent.mt"]] <- PercentageFeatureSet(M_G_CSF_M3, pattern = "^MT-")
M_G_CSF_M4[["percent.mt"]] <- PercentageFeatureSet(M_G_CSF_M4, pattern = "^MT-")
P_H1[["percent.mt"]] <- PercentageFeatureSet(P_H1, pattern = "^MT-")
P_H2[["percent.mt"]] <- PercentageFeatureSet(P_H2, pattern = "^MT-")
P_G_CSF_M1[["percent.mt"]] <- PercentageFeatureSet(P_G_CSF_M1, pattern = "^MT-")
P_G_CSF_M2[["percent.mt"]] <- PercentageFeatureSet(P_G_CSF_M2, pattern = "^MT-")
P_G_CSF_M3[["percent.mt"]] <- PercentageFeatureSet(P_G_CSF_M3, pattern = "^MT-")
P_G_CSF_M4[["percent.mt"]] <- PercentageFeatureSet(P_G_CSF_M4, pattern = "^MT-")

#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m_G_CSF_H1 <- match(HB.genes, rownames(G_CSF_H1@assays$RNA)) 
HB.genes_G_CSF_H1 <- rownames(G_CSF_H1@assays$RNA)[HB_m_G_CSF_H1] 
HB.genes_G_CSF_H1 <- HB.genes_G_CSF_H1[!is.na(HB.genes_G_CSF_H1)] 
G_CSF_H1[["percent.HB"]]<-PercentageFeatureSet(G_CSF_H1, features=HB.genes_G_CSF_H1)

HB_m_G_CSF_H2 <- match(HB.genes, rownames(G_CSF_H2@assays$RNA)) 
HB.genes_G_CSF_H2 <- rownames(G_CSF_H2@assays$RNA)[HB_m_G_CSF_H2] 
HB.genes_G_CSF_H2 <- HB.genes_G_CSF_H2[!is.na(HB.genes_G_CSF_H2)] 
G_CSF_H2[["percent.HB"]]<-PercentageFeatureSet(G_CSF_H2, features=HB.genes_G_CSF_H2)

HB_m_G_CSF_M1 <- match(HB.genes, rownames(G_CSF_M1@assays$RNA)) 
HB.genes_G_CSF_M1 <- rownames(G_CSF_M1@assays$RNA)[HB_m_G_CSF_M1] 
HB.genes_G_CSF_M1 <- HB.genes_G_CSF_M1[!is.na(HB.genes_G_CSF_M1)] 
G_CSF_M1[["percent.HB"]]<-PercentageFeatureSet(G_CSF_M1, features=HB.genes_G_CSF_M1)

HB_m_G_CSF_M2 <- match(HB.genes, rownames(G_CSF_M2@assays$RNA)) 
HB.genes_G_CSF_M2 <- rownames(G_CSF_M2@assays$RNA)[HB_m_G_CSF_M2] 
HB.genes_G_CSF_M2 <- HB.genes_G_CSF_M2[!is.na(HB.genes_G_CSF_M2)] 
G_CSF_M2[["percent.HB"]]<-PercentageFeatureSet(G_CSF_M2, features=HB.genes_G_CSF_M2)

HB_m_G_CSF_M3 <- match(HB.genes, rownames(G_CSF_M3@assays$RNA)) 
HB.genes_G_CSF_M3 <- rownames(G_CSF_M3@assays$RNA)[HB_m_G_CSF_M3] 
HB.genes_G_CSF_M3 <- HB.genes_G_CSF_M3[!is.na(HB.genes_G_CSF_M3)] 
G_CSF_M3[["percent.HB"]]<-PercentageFeatureSet(G_CSF_M3, features=HB.genes_G_CSF_M3)

HB_m_G_CSF_M4 <- match(HB.genes, rownames(G_CSF_M4@assays$RNA)) 
HB.genes_G_CSF_M4 <- rownames(G_CSF_M4@assays$RNA)[HB_m_G_CSF_M4] 
HB.genes_G_CSF_M4 <- HB.genes_G_CSF_M4[!is.na(HB.genes_G_CSF_M4)] 
G_CSF_M4[["percent.HB"]]<-PercentageFeatureSet(G_CSF_M4, features=HB.genes_G_CSF_M4)

HB_m_M_H1 <- match(HB.genes, rownames(M_H1@assays$RNA)) 
HB.genes_M_H1 <- rownames(M_H1@assays$RNA)[HB_m_M_H1] 
HB.genes_M_H1 <- HB.genes_M_H1[!is.na(HB.genes_M_H1)] 
M_H1[["percent.HB"]]<-PercentageFeatureSet(M_H1, features=HB.genes_M_H1)

HB_m_M_H2 <- match(HB.genes, rownames(M_H2@assays$RNA)) 
HB.genes_M_H2 <- rownames(M_H2@assays$RNA)[HB_m_M_H2] 
HB.genes_M_H2 <- HB.genes_M_H2[!is.na(HB.genes_M_H2)] 
M_H2[["percent.HB"]]<-PercentageFeatureSet(M_H2, features=HB.genes_M_H2)

HB_m_M_G_CSF_M1 <- match(HB.genes, rownames(M_G_CSF_M1@assays$RNA)) 
HB.genes_M_G_CSF_M1 <- rownames(M_G_CSF_M1@assays$RNA)[HB_m_M_G_CSF_M1] 
HB.genes_M_G_CSF_M1 <- HB.genes_M_G_CSF_M1[!is.na(HB.genes_M_G_CSF_M1)] 
M_G_CSF_M1[["percent.HB"]]<-PercentageFeatureSet(M_G_CSF_M1, features=HB.genes_M_G_CSF_M1)

HB_m_M_G_CSF_M2 <- match(HB.genes, rownames(M_G_CSF_M2@assays$RNA)) 
HB.genes_M_G_CSF_M2 <- rownames(M_G_CSF_M2@assays$RNA)[HB_m_M_G_CSF_M2] 
HB.genes_M_G_CSF_M2 <- HB.genes_M_G_CSF_M2[!is.na(HB.genes_M_G_CSF_M2)] 
M_G_CSF_M2[["percent.HB"]]<-PercentageFeatureSet(M_G_CSF_M2, features=HB.genes_M_G_CSF_M2)

HB_m_M_G_CSF_M3 <- match(HB.genes, rownames(M_G_CSF_M3@assays$RNA)) 
HB.genes_M_G_CSF_M3 <- rownames(M_G_CSF_M3@assays$RNA)[HB_m_M_G_CSF_M3] 
HB.genes_M_G_CSF_M3 <- HB.genes_M_G_CSF_M3[!is.na(HB.genes_M_G_CSF_M3)] 
M_G_CSF_M3[["percent.HB"]]<-PercentageFeatureSet(M_G_CSF_M3, features=HB.genes_M_G_CSF_M3)

HB_m_M_G_CSF_M4 <- match(HB.genes, rownames(M_G_CSF_M4@assays$RNA)) 
HB.genes_M_G_CSF_M4 <- rownames(M_G_CSF_M4@assays$RNA)[HB_m_M_G_CSF_M4] 
HB.genes_M_G_CSF_M4 <- HB.genes_M_G_CSF_M4[!is.na(HB.genes_M_G_CSF_M4)] 
M_G_CSF_M4[["percent.HB"]]<-PercentageFeatureSet(M_G_CSF_M4, features=HB.genes_M_G_CSF_M4)

HB_m_P_H1 <- match(HB.genes, rownames(P_H1@assays$RNA)) 
HB.genes_P_H1 <- rownames(P_H1@assays$RNA)[HB_m_P_H1] 
HB.genes_P_H1 <- HB.genes_P_H1[!is.na(HB.genes_P_H1)] 
P_H1[["percent.HB"]]<-PercentageFeatureSet(P_H1, features=HB.genes_P_H1)

HB_m_P_H2 <- match(HB.genes, rownames(P_H2@assays$RNA)) 
HB.genes_P_H2 <- rownames(P_H2@assays$RNA)[HB_m_P_H2] 
HB.genes_P_H2 <- HB.genes_P_H2[!is.na(HB.genes_P_H2)] 
P_H2[["percent.HB"]]<-PercentageFeatureSet(P_H2, features=HB.genes_P_H2)

HB_m_P_G_CSF_M1 <- match(HB.genes, rownames(P_G_CSF_M1@assays$RNA)) 
HB.genes_P_G_CSF_M1 <- rownames(P_G_CSF_M1@assays$RNA)[HB_m_P_G_CSF_M1] 
HB.genes_P_G_CSF_M1 <- HB.genes_P_G_CSF_M1[!is.na(HB.genes_P_G_CSF_M1)] 
P_G_CSF_M1[["percent.HB"]]<-PercentageFeatureSet(P_G_CSF_M1, features=HB.genes_P_G_CSF_M1)

HB_m_P_G_CSF_M2 <- match(HB.genes, rownames(P_G_CSF_M2@assays$RNA)) 
HB.genes_P_G_CSF_M2 <- rownames(P_G_CSF_M2@assays$RNA)[HB_m_P_G_CSF_M2] 
HB.genes_P_G_CSF_M2 <- HB.genes_P_G_CSF_M2[!is.na(HB.genes_P_G_CSF_M2)] 
P_G_CSF_M2[["percent.HB"]]<-PercentageFeatureSet(P_G_CSF_M2, features=HB.genes_P_G_CSF_M2)

HB_m_P_G_CSF_M3 <- match(HB.genes, rownames(P_G_CSF_M3@assays$RNA)) 
HB.genes_P_G_CSF_M3 <- rownames(P_G_CSF_M3@assays$RNA)[HB_m_P_G_CSF_M3] 
HB.genes_P_G_CSF_M3 <- HB.genes_P_G_CSF_M3[!is.na(HB.genes_P_G_CSF_M3)] 
P_G_CSF_M3[["percent.HB"]]<-PercentageFeatureSet(P_G_CSF_M3, features=HB.genes_P_G_CSF_M3)

HB_m_P_G_CSF_M4 <- match(HB.genes, rownames(P_G_CSF_M4@assays$RNA)) 
HB.genes_P_G_CSF_M4 <- rownames(P_G_CSF_M4@assays$RNA)[HB_m_P_G_CSF_M4] 
HB.genes_P_G_CSF_M4 <- HB.genes_P_G_CSF_M4[!is.na(HB.genes_P_G_CSF_M4)] 
P_G_CSF_M4[["percent.HB"]]<-PercentageFeatureSet(P_G_CSF_M4, features=HB.genes_P_G_CSF_M4)

####Feature、count、线粒体基因、红细胞基因占比可视化。
violin.G_CSF_H1 <- VlnPlot(G_CSF_H1,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.G_CSF_H2 <- VlnPlot(G_CSF_H2,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.G_CSF_M1 <- VlnPlot(G_CSF_M1,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.G_CSF_M2 <- VlnPlot(G_CSF_M2,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.G_CSF_M3 <- VlnPlot(G_CSF_M3,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.G_CSF_M4 <- VlnPlot(G_CSF_M4,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 

violin.M_H1 <- VlnPlot(M_H1,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.M_H2 <- VlnPlot(M_H2,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.M_G_CSF_M1 <- VlnPlot(M_G_CSF_M1,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.M_G_CSF_M2 <- VlnPlot(M_G_CSF_M2,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.M_G_CSF_M3 <- VlnPlot(M_G_CSF_M3,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.M_G_CSF_M4 <- VlnPlot(M_G_CSF_M4,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4)

violin.P_H1 <- VlnPlot(P_H1,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.P_H2 <- VlnPlot(P_H2,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.P_G_CSF_M1 <- VlnPlot(P_G_CSF_M1,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.P_G_CSF_M2 <- VlnPlot(P_G_CSF_M2,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.P_G_CSF_M3 <- VlnPlot(P_G_CSF_M3,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 
violin.P_G_CSF_M4 <- VlnPlot(P_G_CSF_M4,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"),pt.size = 0.01, ncol = 4) 

#####也可手动保存图片
setwd("D:\\MM\\scRNA.MM\\GSE223972_RAW\\QC\\")
ggsave("vlnplot_before_qc_G_CSF_H1.tiff",dpi = 300, plot = violin.G_CSF_H1, width = 12, height = 6) 
ggsave("vlnplot_before_qc_G_CSF_H2.tiff",dpi = 300, plot = violin.G_CSF_H2, width = 12, height = 6)
ggsave("vlnplot_before_qc_G_CSF_M1.tiff",dpi = 300, plot = violin.G_CSF_M1, width = 12, height = 6)
ggsave("vlnplot_before_qc_G_CSF_M2.tiff",dpi = 300, plot = violin.G_CSF_M2, width = 12, height = 6) 
ggsave("vlnplot_before_qc_G_CSF_M3.tiff",dpi = 300, plot = violin.G_CSF_M3, width = 12, height = 6) 
ggsave("vlnplot_before_qc_G_CSF_M4.tiff",dpi = 300, plot = violin.G_CSF_M4, width = 12, height = 6)

ggsave("vlnplot_before_qc_M_H1.tiff",dpi = 300, plot = violin.M_H1, width = 12, height = 6)
ggsave("vlnplot_before_qc_M_H2.tiff",dpi = 300, plot = violin.M_H2, width = 12, height = 6) 
ggsave("vlnplot_before_qc_M_G_CSF_M1.tiff",dpi = 300, plot = violin.M_G_CSF_M1, width = 12, height = 6) 
ggsave("vlnplot_before_qc_M_G_CSF_M2.tiff",dpi = 300, plot = violin.M_G_CSF_M2, width = 12, height = 6)
ggsave("vlnplot_before_qc_M_G_CSF_M3.tiff",dpi = 300, plot = violin.M_G_CSF_M3, width = 12, height = 6)
ggsave("vlnplot_before_qc_M_G_CSF_M4.tiff",dpi = 300, plot = violin.M_G_CSF_M4, width = 12, height = 6) 

ggsave("vlnplot_before_qc_P_H1.tiff",dpi = 300, plot = violin.P_H1, width = 12, height = 6) 
ggsave("vlnplot_before_qc_P_H2.tiff",dpi = 300, plot = violin.P_H2, width = 12, height = 6)
ggsave("vlnplot_before_qc_P_G_CSF_M1.tiff",dpi = 300, plot = violin.P_G_CSF_M1, width = 12, height = 6)
ggsave("vlnplot_before_qc_P_G_CSF_M2.tiff",dpi = 300, plot = violin.P_G_CSF_M2, width = 12, height = 6) 
ggsave("vlnplot_before_qc_P_G_CSF_M3.tiff",dpi = 300, plot = violin.P_G_CSF_M3, width = 12, height = 6) 
ggsave("vlnplot_before_qc_P_G_CSF_M4.tiff",dpi = 300, plot = violin.P_G_CSF_M4, width = 12, height = 6)

###这几个指标之间的相关性。 把图画到画板上，然后手动保存
plot.G_CSF_H1_1=FeatureScatter(G_CSF_H1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.G_CSF_H1_2=FeatureScatter(G_CSF_H1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.G_CSF_H1_3=FeatureScatter(G_CSF_H1, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.G_CSF_H1 <- CombinePlots(plots = list(plot.G_CSF_H1_1, plot.G_CSF_H1_2, plot.G_CSF_H1_3), nrow=1, legend="none") 
ggsave("pearplot.G_CSF_H1.tiff",dpi = 300, plot = pearplot.G_CSF_H1, width = 12, height = 6)

plot.G_CSF_H2_1=FeatureScatter(G_CSF_H2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.G_CSF_H2_2=FeatureScatter(G_CSF_H2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.G_CSF_H2_3=FeatureScatter(G_CSF_H2, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.G_CSF_H2 <- CombinePlots(plots = list(plot.G_CSF_H2_1, plot.G_CSF_H2_2, plot.G_CSF_H2_3), nrow=1, legend="none") 
ggsave("pearplot.G_CSF_H2.tiff",dpi = 300, plot = pearplot.G_CSF_H2, width = 12, height = 6)

plot.G_CSF_M1_1=FeatureScatter(G_CSF_M1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.G_CSF_M1_2=FeatureScatter(G_CSF_M1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.G_CSF_M1_3=FeatureScatter(G_CSF_M1, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.G_CSF_M1 <- CombinePlots(plots = list(plot.G_CSF_M1_1, plot.G_CSF_M1_2, plot.G_CSF_M1_3), nrow=1, legend="none") 
ggsave("pearplot.G_CSF_M1.tiff",dpi = 300, plot = pearplot.G_CSF_M1, width = 12, height = 6)

plot.G_CSF_M2_1=FeatureScatter(G_CSF_M2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.G_CSF_M2_2=FeatureScatter(G_CSF_M2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.G_CSF_M2_3=FeatureScatter(G_CSF_M2, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.G_CSF_M2 <- CombinePlots(plots = list(plot.G_CSF_M2_1, plot.G_CSF_M2_2, plot.G_CSF_M2_3), nrow=1, legend="none") 
ggsave("pearplot.G_CSF_M2.tiff",dpi = 300, plot = pearplot.G_CSF_M2, width = 12, height = 6)

plot.G_CSF_M3_1=FeatureScatter(G_CSF_M3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.G_CSF_M3_2=FeatureScatter(G_CSF_M3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.G_CSF_M3_3=FeatureScatter(G_CSF_M3, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.G_CSF_M3 <- CombinePlots(plots = list(plot.G_CSF_M3_1, plot.G_CSF_M3_2, plot.G_CSF_M3_3), nrow=1, legend="none") 
ggsave("pearplot.G_CSF_M3.tiff",dpi = 300, plot = pearplot.G_CSF_M3, width = 12, height = 6)

plot.G_CSF_M4_1=FeatureScatter(G_CSF_M4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.G_CSF_M4_2=FeatureScatter(G_CSF_M4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.G_CSF_M4_3=FeatureScatter(G_CSF_M4, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.G_CSF_M4 <- CombinePlots(plots = list(plot.G_CSF_M4_1, plot.G_CSF_M4_2, plot.G_CSF_M4_3), nrow=1, legend="none") 
ggsave("pearplot.G_CSF_M4.tiff",dpi = 300, plot = pearplot.G_CSF_M4, width = 12, height = 6)

plot.M_H1_1=FeatureScatter(M_H1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.M_H1_2=FeatureScatter(M_H1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.M_H1_3=FeatureScatter(M_H1, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.M_H1 <- CombinePlots(plots = list(plot.M_H1_1, plot.M_H1_2, plot.M_H1_3), nrow=1, legend="none") 
ggsave("pearplot.M_H1.tiff",dpi = 300, plot = pearplot.M_H1, width = 12, height = 6)

plot.M_H2_1=FeatureScatter(M_H2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.M_H2_2=FeatureScatter(M_H2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.M_H2_3=FeatureScatter(M_H2, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.M_H2 <- CombinePlots(plots = list(plot.M_H2_1, plot.M_H2_2, plot.M_H2_3), nrow=1, legend="none") 
ggsave("pearplot.M_H2.tiff",dpi = 300, plot = pearplot.M_H2, width = 12, height = 6)

plot.M_G_CSF_M1_1=FeatureScatter(M_G_CSF_M1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.M_G_CSF_M1_2=FeatureScatter(M_G_CSF_M1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.M_G_CSF_M1_3=FeatureScatter(M_G_CSF_M1, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.M_G_CSF_M1 <- CombinePlots(plots = list(plot.M_G_CSF_M1_1, plot.M_G_CSF_M1_2, plot.M_G_CSF_M1_3), nrow=1, legend="none") 
ggsave("pearplot.M_G_CSF_M1.tiff",dpi = 300, plot = pearplot.M_G_CSF_M1, width = 12, height = 6)

plot.M_G_CSF_M2_1=FeatureScatter(M_G_CSF_M2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.M_G_CSF_M2_2=FeatureScatter(M_G_CSF_M2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.M_G_CSF_M2_3=FeatureScatter(M_G_CSF_M2, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.M_G_CSF_M2 <- CombinePlots(plots = list(plot.M_G_CSF_M2_1, plot.M_G_CSF_M2_2, plot.M_G_CSF_M2_3), nrow=1, legend="none") 
ggsave("pearplot.M_G_CSF_M2.tiff",dpi = 300, plot = pearplot.M_G_CSF_M2, width = 12, height = 6)

plot.M_G_CSF_M3_1=FeatureScatter(M_G_CSF_M3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.M_G_CSF_M3_2=FeatureScatter(M_G_CSF_M3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.M_G_CSF_M3_3=FeatureScatter(M_G_CSF_M3, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.M_G_CSF_M3 <- CombinePlots(plots = list(plot.M_G_CSF_M3_1, plot.M_G_CSF_M3_2, plot.M_G_CSF_M3_3), nrow=1, legend="none") 
ggsave("pearplot.M_G_CSF_M3.tiff",dpi = 300, plot = pearplot.M_G_CSF_M3, width = 12, height = 6)

plot.M_G_CSF_M4_1=FeatureScatter(M_G_CSF_M4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.M_G_CSF_M4_2=FeatureScatter(M_G_CSF_M4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.M_G_CSF_M4_3=FeatureScatter(M_G_CSF_M4, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.M_G_CSF_M4 <- CombinePlots(plots = list(plot.M_G_CSF_M4_1, plot.M_G_CSF_M4_2, plot.M_G_CSF_M4_3), nrow=1, legend="none") 
ggsave("pearplot.M_G_CSF_M4.tiff",dpi = 300, plot = pearplot.M_G_CSF_M4, width = 12, height = 6)

plot.P_H1_1=FeatureScatter(P_H1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.P_H1_2=FeatureScatter(P_H1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.P_H1_3=FeatureScatter(P_H1, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.P_H1 <- CombinePlots(plots = list(plot.P_H1_1, plot.P_H1_2, plot.P_H1_3), nrow=1, legend="none") 
ggsave("pearplot.P_H1.tiff",dpi = 300, plot = pearplot.P_H1, width = 12, height = 6)

plot.P_H2_1=FeatureScatter(P_H2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.P_H2_2=FeatureScatter(P_H2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.P_H2_3=FeatureScatter(P_H2, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.P_H2 <- CombinePlots(plots = list(plot.P_H2_1, plot.P_H2_2, plot.P_H2_3), nrow=1, legend="none") 
ggsave("pearplot.P_H2.tiff",dpi = 300, plot = pearplot.P_H2, width = 12, height = 6)

plot.P_G_CSF_M1_1=FeatureScatter(P_G_CSF_M1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.P_G_CSF_M1_2=FeatureScatter(P_G_CSF_M1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.P_G_CSF_M1_3=FeatureScatter(P_G_CSF_M1, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.P_G_CSF_M1 <- CombinePlots(plots = list(plot.P_G_CSF_M1_1, plot.P_G_CSF_M1_2, plot.P_G_CSF_M1_3), nrow=1, legend="none") 
ggsave("pearplot.P_G_CSF_M1.tiff",dpi = 300, plot = pearplot.P_G_CSF_M1, width = 12, height = 6)

plot.P_G_CSF_M2_1=FeatureScatter(P_G_CSF_M2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.P_G_CSF_M2_2=FeatureScatter(P_G_CSF_M2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.P_G_CSF_M2_3=FeatureScatter(P_G_CSF_M2, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.P_G_CSF_M2 <- CombinePlots(plots = list(plot.P_G_CSF_M2_1, plot.P_G_CSF_M2_2, plot.P_G_CSF_M2_3), nrow=1, legend="none") 
ggsave("pearplot.P_G_CSF_M2.tiff",dpi = 300, plot = pearplot.P_G_CSF_M2, width = 12, height = 6)

plot.P_G_CSF_M3_1=FeatureScatter(P_G_CSF_M3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.P_G_CSF_M3_2=FeatureScatter(P_G_CSF_M3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.P_G_CSF_M3_3=FeatureScatter(P_G_CSF_M3, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.P_G_CSF_M3 <- CombinePlots(plots = list(plot.P_G_CSF_M3_1, plot.P_G_CSF_M3_2, plot.P_G_CSF_M3_3), nrow=1, legend="none") 
ggsave("pearplot.P_G_CSF_M3.tiff",dpi = 300, plot = pearplot.P_G_CSF_M3, width = 12, height = 6)

plot.P_G_CSF_M4_1=FeatureScatter(P_G_CSF_M4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot.P_G_CSF_M4_2=FeatureScatter(P_G_CSF_M4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot.P_G_CSF_M4_3=FeatureScatter(P_G_CSF_M4, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot.P_G_CSF_M4 <- CombinePlots(plots = list(plot.P_G_CSF_M4_1, plot.P_G_CSF_M4_2, plot.P_G_CSF_M4_3), nrow=1, legend="none") 
ggsave("pearplot.P_G_CSF_M4.tiff",dpi = 300, plot = pearplot.P_G_CSF_M4, width = 12, height = 6)

# filter everything to 500 unique genes/cell
G_CSF_H1_1 <- subset(G_CSF_H1, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
G_CSF_H2_1 <- subset(G_CSF_H2, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
G_CSF_M1_1 <- subset(G_CSF_M1, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
G_CSF_M2_1 <- subset(G_CSF_M2, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
G_CSF_M3_1 <- subset(G_CSF_M3, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
G_CSF_M4_1 <- subset(G_CSF_M4, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)

M_H1_1 <- subset(M_H1, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
M_H2_1 <- subset(M_H2, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
M_G_CSF_M1_1 <- subset(M_G_CSF_M1, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
M_G_CSF_M2_1 <- subset(M_G_CSF_M2, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
M_G_CSF_M3_1 <- subset(M_G_CSF_M3, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
M_G_CSF_M4_1 <- subset(M_G_CSF_M4, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)

P_H1_1 <- subset(P_H1, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
P_H2_1 <- subset(P_H2, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
P_G_CSF_M1_1 <- subset(P_G_CSF_M1, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
P_G_CSF_M2_1 <- subset(P_G_CSF_M2, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
P_G_CSF_M3_1 <- subset(P_G_CSF_M3, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)
P_G_CSF_M4_1 <- subset(P_G_CSF_M4, subset = nFeature_RNA > 500 & percent.mt < 5 & percent.HB < 1 & nCount_RNA > 1000)

# Normalize and make UMAP
G_CSF_H1_1 <- NormalizeData(G_CSF_H1_1, normalization.method = "LogNormalize", 
                            scale.factor = 10000)
G_CSF_H1_1 <- FindVariableFeatures(G_CSF_H1_1, selection.method = "vst", nfeatures = 2000)
all.genes_G_CSF_H1_1 <- rownames(G_CSF_H1_1)
G_CSF_H1_1 <- ScaleData(G_CSF_H1_1, features = all.genes_G_CSF_H1_1)
G_CSF_H1_1 <- RunPCA(G_CSF_H1_1, features = VariableFeatures(object = G_CSF_H1_1))
G_CSF_H1_1 <- FindNeighbors(G_CSF_H1_1, dims = 1:20)
G_CSF_H1_1 <- FindClusters(G_CSF_H1_1, resolution = 0.5)
G_CSF_H1_1 <- RunUMAP(G_CSF_H1_1, dims = 1:20)

G_CSF_H2_1 <- NormalizeData(G_CSF_H2_1, normalization.method = "LogNormalize", 
                            scale.factor = 10000)
G_CSF_H2_1 <- FindVariableFeatures(G_CSF_H2_1, selection.method = "vst", nfeatures = 2000)
all.genes_G_CSF_H2_1 <- rownames(G_CSF_H2_1)
G_CSF_H2_1 <- ScaleData(G_CSF_H2_1, features = all.genes_G_CSF_H2_1)
G_CSF_H2_1 <- RunPCA(G_CSF_H2_1, features = VariableFeatures(object = G_CSF_H2_1))
G_CSF_H2_1 <- FindNeighbors(G_CSF_H2_1, dims = 1:20)
G_CSF_H2_1 <- FindClusters(G_CSF_H2_1, resolution = 0.5)
G_CSF_H2_1 <- RunUMAP(G_CSF_H2_1, dims = 1:20)

G_CSF_M1_1 <- NormalizeData(G_CSF_M1_1, normalization.method = "LogNormalize", 
                            scale.factor = 10000)
G_CSF_M1_1 <- FindVariableFeatures(G_CSF_M1_1, selection.method = "vst", nfeatures = 2000)
all.genes_G_CSF_M1_1 <- rownames(G_CSF_M1_1)
G_CSF_M1_1 <- ScaleData(G_CSF_M1_1, features = all.genes_G_CSF_M1_1)
G_CSF_M1_1 <- RunPCA(G_CSF_M1_1, features = VariableFeatures(object = G_CSF_M1_1))
G_CSF_M1_1 <- FindNeighbors(G_CSF_M1_1, dims = 1:20)
G_CSF_M1_1 <- FindClusters(G_CSF_M1_1, resolution = 0.5)
G_CSF_M1_1 <- RunUMAP(G_CSF_M1_1, dims = 1:20)

G_CSF_M2_1 <- NormalizeData(G_CSF_M2_1, normalization.method = "LogNormalize", 
                            scale.factor = 10000)
G_CSF_M2_1 <- FindVariableFeatures(G_CSF_M2_1, selection.method = "vst", nfeatures = 2000)
all.genes_G_CSF_M2_1 <- rownames(G_CSF_M2_1)
G_CSF_M2_1 <- ScaleData(G_CSF_M2_1, features = all.genes_G_CSF_M2_1)
G_CSF_M2_1 <- RunPCA(G_CSF_M2_1, features = VariableFeatures(object = G_CSF_M2_1))
G_CSF_M2_1 <- FindNeighbors(G_CSF_M2_1, dims = 1:20)
G_CSF_M2_1 <- FindClusters(G_CSF_M2_1, resolution = 0.5)
G_CSF_M2_1 <- RunUMAP(G_CSF_M2_1, dims = 1:20)

G_CSF_M3_1 <- NormalizeData(G_CSF_M3_1, normalization.method = "LogNormalize", 
                            scale.factor = 10000)
G_CSF_M3_1 <- FindVariableFeatures(G_CSF_M3_1, selection.method = "vst", nfeatures = 2000)
all.genes_G_CSF_M3_1 <- rownames(G_CSF_M3_1)
G_CSF_M3_1 <- ScaleData(G_CSF_M3_1, features = all.genes_G_CSF_M3_1)
G_CSF_M3_1 <- RunPCA(G_CSF_M3_1, features = VariableFeatures(object = G_CSF_M3_1))
G_CSF_M3_1 <- FindNeighbors(G_CSF_M3_1, dims = 1:20)
G_CSF_M3_1 <- FindClusters(G_CSF_M3_1, resolution = 0.5)
G_CSF_M3_1 <- RunUMAP(G_CSF_M3_1, dims = 1:20)

G_CSF_M4_1 <- NormalizeData(G_CSF_M4_1, normalization.method = "LogNormalize", 
                            scale.factor = 10000)
G_CSF_M4_1 <- FindVariableFeatures(G_CSF_M4_1, selection.method = "vst", nfeatures = 2000)
all.genes_G_CSF_M4_1 <- rownames(G_CSF_M4_1)
G_CSF_M4_1 <- ScaleData(G_CSF_M4_1, features = all.genes_G_CSF_M4_1)
G_CSF_M4_1 <- RunPCA(G_CSF_M4_1, features = VariableFeatures(object = G_CSF_M4_1))
G_CSF_M4_1 <- FindNeighbors(G_CSF_M4_1, dims = 1:20)
G_CSF_M4_1 <- FindClusters(G_CSF_M4_1, resolution = 0.5)
G_CSF_M4_1 <- RunUMAP(G_CSF_M4_1, dims = 1:20)

M_H1_1 <- NormalizeData(M_H1_1, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
M_H1_1 <- FindVariableFeatures(M_H1_1, selection.method = "vst", nfeatures = 2000)
all.genes_M_H1_1 <- rownames(M_H1_1)
M_H1_1 <- ScaleData(M_H1_1, features = all.genes_M_H1_1)
M_H1_1 <- RunPCA(M_H1_1, features = VariableFeatures(object = M_H1_1))
M_H1_1 <- FindNeighbors(M_H1_1, dims = 1:20)
M_H1_1 <- FindClusters(M_H1_1, resolution = 0.5)
M_H1_1 <- RunUMAP(M_H1_1, dims = 1:20)

M_H2_1 <- NormalizeData(M_H2_1, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
M_H2_1 <- FindVariableFeatures(M_H2_1, selection.method = "vst", nfeatures = 2000)
all.genes_M_H2_1 <- rownames(M_H2_1)
M_H2_1 <- ScaleData(M_H2_1, features = all.genes_M_H2_1)
M_H2_1 <- RunPCA(M_H2_1, features = VariableFeatures(object = M_H2_1))
M_H2_1 <- FindNeighbors(M_H2_1, dims = 1:20)
M_H2_1 <- FindClusters(M_H2_1, resolution = 0.5)
M_H2_1 <- RunUMAP(M_H2_1, dims = 1:20)

M_G_CSF_M1_1 <- NormalizeData(M_G_CSF_M1_1, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
M_G_CSF_M1_1 <- FindVariableFeatures(M_G_CSF_M1_1, selection.method = "vst", nfeatures = 2000)
all.genes_M_G_CSF_M1_1 <- rownames(M_G_CSF_M1_1)
M_G_CSF_M1_1 <- ScaleData(M_G_CSF_M1_1, features = all.genes_M_G_CSF_M1_1)
M_G_CSF_M1_1 <- RunPCA(M_G_CSF_M1_1, features = VariableFeatures(object = M_G_CSF_M1_1))
M_G_CSF_M1_1 <- FindNeighbors(M_G_CSF_M1_1, dims = 1:20)
M_G_CSF_M1_1 <- FindClusters(M_G_CSF_M1_1, resolution = 0.5)
M_G_CSF_M1_1 <- RunUMAP(M_G_CSF_M1_1, dims = 1:20)

M_G_CSF_M2_1 <- NormalizeData(M_G_CSF_M2_1, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
M_G_CSF_M2_1 <- FindVariableFeatures(M_G_CSF_M2_1, selection.method = "vst", nfeatures = 2000)
all.genes_M_G_CSF_M2_1 <- rownames(M_G_CSF_M2_1)
M_G_CSF_M2_1 <- ScaleData(M_G_CSF_M2_1, features = all.genes_M_G_CSF_M2_1)
M_G_CSF_M2_1 <- RunPCA(M_G_CSF_M2_1, features = VariableFeatures(object = M_G_CSF_M2_1))
M_G_CSF_M2_1 <- FindNeighbors(M_G_CSF_M2_1, dims = 1:20)
M_G_CSF_M2_1 <- FindClusters(M_G_CSF_M2_1, resolution = 0.5)
M_G_CSF_M2_1 <- RunUMAP(M_G_CSF_M2_1, dims = 1:20)

M_G_CSF_M3_1 <- NormalizeData(M_G_CSF_M3_1, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
M_G_CSF_M3_1 <- FindVariableFeatures(M_G_CSF_M3_1, selection.method = "vst", nfeatures = 2000)
all.genes_M_G_CSF_M3_1 <- rownames(M_G_CSF_M3_1)
M_G_CSF_M3_1 <- ScaleData(M_G_CSF_M3_1, features = all.genes_M_G_CSF_M3_1)
M_G_CSF_M3_1 <- RunPCA(M_G_CSF_M3_1, features = VariableFeatures(object = M_G_CSF_M3_1))
M_G_CSF_M3_1 <- FindNeighbors(M_G_CSF_M3_1, dims = 1:20)
M_G_CSF_M3_1 <- FindClusters(M_G_CSF_M3_1, resolution = 0.5)
M_G_CSF_M3_1 <- RunUMAP(M_G_CSF_M3_1, dims = 1:20)

M_G_CSF_M4_1 <- NormalizeData(M_G_CSF_M4_1, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
M_G_CSF_M4_1 <- FindVariableFeatures(M_G_CSF_M4_1, selection.method = "vst", nfeatures = 2000)
all.genes_M_G_CSF_M4_1 <- rownames(M_G_CSF_M4_1)
M_G_CSF_M4_1 <- ScaleData(M_G_CSF_M4_1, features = all.genes_M_G_CSF_M4_1)
M_G_CSF_M4_1 <- RunPCA(M_G_CSF_M4_1, features = VariableFeatures(object = M_G_CSF_M4_1))
M_G_CSF_M4_1 <- FindNeighbors(M_G_CSF_M4_1, dims = 1:20)
M_G_CSF_M4_1 <- FindClusters(M_G_CSF_M4_1, resolution = 0.5)
M_G_CSF_M4_1 <- RunUMAP(M_G_CSF_M4_1, dims = 1:20)

P_H1_1 <- NormalizeData(P_H1_1, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
P_H1_1 <- FindVariableFeatures(P_H1_1, selection.method = "vst", nfeatures = 2000)
all.genes_P_H1_1 <- rownames(P_H1_1)
P_H1_1 <- ScaleData(P_H1_1, features = all.genes_P_H1_1)
P_H1_1 <- RunPCA(P_H1_1, features = VariableFeatures(object = P_H1_1))
P_H1_1 <- FindNeighbors(P_H1_1, dims = 1:20)
P_H1_1 <- FindClusters(P_H1_1, resolution = 0.5)
P_H1_1 <- RunUMAP(P_H1_1, dims = 1:20)

P_H2_1 <- NormalizeData(P_H2_1, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
P_H2_1 <- FindVariableFeatures(P_H2_1, selection.method = "vst", nfeatures = 2000)
all.genes_P_H2_1 <- rownames(P_H2_1)
P_H2_1 <- ScaleData(P_H2_1, features = all.genes_P_H2_1)
P_H2_1 <- RunPCA(P_H2_1, features = VariableFeatures(object = P_H2_1))
P_H2_1 <- FindNeighbors(P_H2_1, dims = 1:20)
P_H2_1 <- FindClusters(P_H2_1, resolution = 0.5)
P_H2_1 <- RunUMAP(P_H2_1, dims = 1:20)

P_G_CSF_M1_1 <- NormalizeData(P_G_CSF_M1_1, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
P_G_CSF_M1_1 <- FindVariableFeatures(P_G_CSF_M1_1, selection.method = "vst", nfeatures = 2000)
all.genes_P_G_CSF_M1_1 <- rownames(P_G_CSF_M1_1)
P_G_CSF_M1_1 <- ScaleData(P_G_CSF_M1_1, features = all.genes_P_G_CSF_M1_1)
P_G_CSF_M1_1 <- RunPCA(P_G_CSF_M1_1, features = VariableFeatures(object = P_G_CSF_M1_1))
P_G_CSF_M1_1 <- FindNeighbors(P_G_CSF_M1_1, dims = 1:20)
P_G_CSF_M1_1 <- FindClusters(P_G_CSF_M1_1, resolution = 0.5)
P_G_CSF_M1_1 <- RunUMAP(P_G_CSF_M1_1, dims = 1:20)

P_G_CSF_M2_1 <- NormalizeData(P_G_CSF_M2_1, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
P_G_CSF_M2_1 <- FindVariableFeatures(P_G_CSF_M2_1, selection.method = "vst", nfeatures = 2000)
all.genes_P_G_CSF_M2_1 <- rownames(P_G_CSF_M2_1)
P_G_CSF_M2_1 <- ScaleData(P_G_CSF_M2_1, features = all.genes_P_G_CSF_M2_1)
P_G_CSF_M2_1 <- RunPCA(P_G_CSF_M2_1, features = VariableFeatures(object = P_G_CSF_M2_1))
P_G_CSF_M2_1 <- FindNeighbors(P_G_CSF_M2_1, dims = 1:20)
P_G_CSF_M2_1 <- FindClusters(P_G_CSF_M2_1, resolution = 0.5)
P_G_CSF_M2_1 <- RunUMAP(P_G_CSF_M2_1, dims = 1:20)

P_G_CSF_M3_1 <- NormalizeData(P_G_CSF_M3_1, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
P_G_CSF_M3_1 <- FindVariableFeatures(P_G_CSF_M3_1, selection.method = "vst", nfeatures = 2000)
all.genes_P_G_CSF_M3_1 <- rownames(P_G_CSF_M3_1)
P_G_CSF_M3_1 <- ScaleData(P_G_CSF_M3_1, features = all.genes_P_G_CSF_M3_1)
P_G_CSF_M3_1 <- RunPCA(P_G_CSF_M3_1, features = VariableFeatures(object = P_G_CSF_M3_1))
P_G_CSF_M3_1 <- FindNeighbors(P_G_CSF_M3_1, dims = 1:20)
P_G_CSF_M3_1 <- FindClusters(P_G_CSF_M3_1, resolution = 0.5)
P_G_CSF_M3_1 <- RunUMAP(P_G_CSF_M3_1, dims = 1:20)

P_G_CSF_M4_1 <- NormalizeData(P_G_CSF_M4_1, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
P_G_CSF_M4_1 <- FindVariableFeatures(P_G_CSF_M4_1, selection.method = "vst", nfeatures = 2000)
all.genes_P_G_CSF_M4_1 <- rownames(P_G_CSF_M4_1)
P_G_CSF_M4_1 <- ScaleData(P_G_CSF_M4_1, features = all.genes_P_G_CSF_M4_1)
P_G_CSF_M4_1 <- RunPCA(P_G_CSF_M4_1, features = VariableFeatures(object = P_G_CSF_M4_1))
P_G_CSF_M4_1 <- FindNeighbors(P_G_CSF_M4_1, dims = 1:20)
P_G_CSF_M4_1 <- FindClusters(P_G_CSF_M4_1, resolution = 0.5)
P_G_CSF_M4_1 <- RunUMAP(P_G_CSF_M4_1, dims = 1:20)

## Run doublet finder , plot pre and post doublet finder results,and save data.
nExp_poi_G_CSF_H1_1 <- round(0.08*length(G_CSF_H1_1@meta.data$orig.ident)*length(G_CSF_H1_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
G_CSF_H1_1_DF <- doubletFinder_v3(G_CSF_H1_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_G_CSF_H1_1, reuse.pANN = FALSE, sct = FALSE)
print(head(G_CSF_H1_1_DF@meta.data))
table(G_CSF_H1_1_DF$DF.classifications_0.25_0.09_423)

DimPlotG_CSF_H1 <- DimPlot(G_CSF_H1_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_G_CSF_H1.tiff",dpi = 300, plot = DimPlotG_CSF_H1, width = 6, height = 6)
G_CSF_H1_1s=subset(G_CSF_H1_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.G_CSF_H1.1 <- DimPlot(G_CSF_H1_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_G_CSF_H1.tiff",dpi = 300, plot = DimPlot.G_CSF_H1.1, width = 6, height = 6)
saveRDS(G_CSF_H1_1s, file = "G_CSF_H1_1s.rds")

nExp_poi_G_CSF_H2_1 <- round(0.08*length(G_CSF_H2_1@meta.data$orig.ident)*length(G_CSF_H2_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
G_CSF_H2_1_DF <- doubletFinder_v3(G_CSF_H2_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_G_CSF_H2_1, reuse.pANN = FALSE, sct = FALSE)
print(head(G_CSF_H2_1_DF@meta.data))
table(G_CSF_H2_1_DF$DF.classifications_0.25_0.09_0)

DimPlotG_CSF_H2 <- DimPlot(G_CSF_H2_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_G_CSF_H2.tiff",dpi = 300, plot = DimPlotG_CSF_H2, width = 6, height = 6)
G_CSF_H2_1s=subset(G_CSF_H2_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.G_CSF_H2.1 <- DimPlot(G_CSF_H2_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_G_CSF_H2.tiff",dpi = 300, plot = DimPlot.G_CSF_H2.1, width = 6, height = 6)
saveRDS(G_CSF_H2_1s, file = "G_CSF_H2_1s.rds")

nExp_poi_G_CSF_M1_1 <- round(0.08*length(G_CSF_M1_1@meta.data$orig.ident)*length(G_CSF_M1_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
G_CSF_M1_1_DF <- doubletFinder_v3(G_CSF_M1_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_G_CSF_M1_1, reuse.pANN = FALSE, sct = FALSE)
print(head(G_CSF_M1_1_DF@meta.data))
table(G_CSF_M1_1_DF$DF.classifications_0.25_0.09_0)

DimPlotG_CSF_M1 <- DimPlot(G_CSF_M1_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_G_CSF_M1.tiff",dpi = 300, plot = DimPlotG_CSF_M1, width = 6, height = 6)
G_CSF_M1_1s=subset(G_CSF_M1_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.G_CSF_M1.1 <- DimPlot(G_CSF_M1_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_G_CSF_M1.tiff",dpi = 300, plot = DimPlot.G_CSF_M1.1, width = 6, height = 6)
saveRDS(G_CSF_M1_1s, file = "G_CSF_M1_1s.rds")

nExp_poi_G_CSF_M2_1 <- round(0.08*length(G_CSF_M2_1@meta.data$orig.ident)*length(G_CSF_M2_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
G_CSF_M2_1_DF <- doubletFinder_v3(G_CSF_M2_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_G_CSF_M2_1, reuse.pANN = FALSE, sct = FALSE)
print(head(G_CSF_M2_1_DF@meta.data))
table(G_CSF_M2_1_DF$DF.classifications_0.25_0.09_0)

DimPlotG_CSF_M2 <- DimPlot(G_CSF_M2_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_G_CSF_M2.tiff",dpi = 300, plot = DimPlotG_CSF_M2, width = 6, height = 6)
G_CSF_M2_1s=subset(G_CSF_M2_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.G_CSF_M2.1 <- DimPlot(G_CSF_M2_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_G_CSF_M2.tiff",dpi = 300, plot = DimPlot.G_CSF_M2.1, width = 6, height = 6)
saveRDS(G_CSF_M2_1s, file = "G_CSF_M2_1s.rds")

nExp_poi_G_CSF_M3_1 <- round(0.08*length(G_CSF_M3_1@meta.data$orig.ident)*length(G_CSF_M3_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
G_CSF_M3_1_DF <- doubletFinder_v3(G_CSF_M3_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_G_CSF_M3_1, reuse.pANN = FALSE, sct = FALSE)
print(head(G_CSF_M3_1_DF@meta.data))
table(G_CSF_M3_1_DF$DF.classifications_0.25_0.09_0)

DimPlotG_CSF_M3 <- DimPlot(G_CSF_M3_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_G_CSF_M3.tiff",dpi = 300, plot = DimPlotG_CSF_M3, width = 6, height = 6)
G_CSF_M3_1s=subset(G_CSF_M3_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.G_CSF_M3.1 <- DimPlot(G_CSF_M3_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_G_CSF_M3.tiff",dpi = 300, plot = DimPlot.G_CSF_M3.1, width = 6, height = 6)
saveRDS(G_CSF_M3_1s, file = "G_CSF_M3_1s.rds")

nExp_poi_G_CSF_M4_1 <- round(0.08*length(G_CSF_M4_1@meta.data$orig.ident)*length(G_CSF_M4_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
G_CSF_M4_1_DF <- doubletFinder_v3(G_CSF_M4_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_G_CSF_M4_1, reuse.pANN = FALSE, sct = FALSE)
print(head(G_CSF_M4_1_DF@meta.data))
table(G_CSF_M4_1_DF$DF.classifications_0.25_0.09_0)

DimPlotG_CSF_M4 <- DimPlot(G_CSF_M4_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_G_CSF_M4.tiff",dpi = 300, plot = DimPlotG_CSF_M4, width = 6, height = 6)
G_CSF_M4_1s=subset(G_CSF_M4_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.G_CSF_M4.1 <- DimPlot(G_CSF_M4_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_G_CSF_M4.tiff",dpi = 300, plot = DimPlot.G_CSF_M4.1, width = 6, height = 6)
saveRDS(G_CSF_M4_1s, file = "G_CSF_M4_1s.rds")

nExp_poi_M_H1_1 <- round(0.08*length(M_H1_1@meta.data$orig.ident)*length(M_H1_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
M_H1_1_DF <- doubletFinder_v3(M_H1_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_M_H1_1, reuse.pANN = FALSE, sct = FALSE)
print(head(M_H1_1_DF@meta.data))
table(M_H1_1_DF$DF.classifications_0.25_0.09_0)

DimPlotM_H1 <- DimPlot(M_H1_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_M_H1.tiff",dpi = 300, plot = DimPlotM_H1, width = 6, height = 6)
M_H1_1s=subset(M_H1_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.M_H1.1 <- DimPlot(M_H1_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_M_H1.tiff",dpi = 300, plot = DimPlot.M_H1.1, width = 6, height = 6)
saveRDS(M_H1_1s, file = "M_H1_1s.rds")

nExp_poi_M_H2_1 <- round(0.08*length(M_H2_1@meta.data$orig.ident)*length(M_H2_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
M_H2_1_DF <- doubletFinder_v3(M_H2_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_M_H2_1, reuse.pANN = FALSE, sct = FALSE)
print(head(M_H2_1_DF@meta.data))
table(M_H2_1_DF$DF.classifications_0.25_0.09_0)

DimPlotM_H2 <- DimPlot(M_H2_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_M_H2.tiff",dpi = 300, plot = DimPlotM_H2, width = 6, height = 6)
M_H2_1s=subset(M_H2_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.M_H2.1 <- DimPlot(M_H2_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_M_H2.tiff",dpi = 300, plot = DimPlot.M_H2.1, width = 6, height = 6)
saveRDS(M_H2_1s, file = "M_H2_1s.rds")

nExp_poi_M_G_CSF_M1_1 <- round(0.08*length(M_G_CSF_M1_1@meta.data$orig.ident)*length(M_G_CSF_M1_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
M_G_CSF_M1_1_DF <- doubletFinder_v3(M_G_CSF_M1_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_M_G_CSF_M1_1, reuse.pANN = FALSE, sct = FALSE)
print(head(M_G_CSF_M1_1_DF@meta.data))
table(M_G_CSF_M1_1_DF$DF.classifications_0.25_0.09_0)

DimPlotM_G_CSF_M1 <- DimPlot(M_G_CSF_M1_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_M_G_CSF_M1.tiff",dpi = 300, plot = DimPlotM_G_CSF_M1, width = 6, height = 6)
M_G_CSF_M1_1s=subset(M_G_CSF_M1_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.M_G_CSF_M1.1 <- DimPlot(M_G_CSF_M1_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_M_G_CSF_M1.tiff",dpi = 300, plot = DimPlot.M_G_CSF_M1.1, width = 6, height = 6)
saveRDS(M_G_CSF_M1_1s, file = "M_G_CSF_M1_1s.rds")

nExp_poi_M_G_CSF_M2_1 <- round(0.08*length(M_G_CSF_M2_1@meta.data$orig.ident)*length(M_G_CSF_M2_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
M_G_CSF_M2_1_DF <- doubletFinder_v3(M_G_CSF_M2_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_M_G_CSF_M2_1, reuse.pANN = FALSE, sct = FALSE)
print(head(M_G_CSF_M2_1_DF@meta.data))
table(M_G_CSF_M2_1_DF$DF.classifications_0.25_0.09_0)

DimPlotM_G_CSF_M2 <- DimPlot(M_G_CSF_M2_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_M_G_CSF_M2.tiff",dpi = 300, plot = DimPlotM_G_CSF_M2, width = 6, height = 6)
M_G_CSF_M2_1s=subset(M_G_CSF_M2_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.M_G_CSF_M2.1 <- DimPlot(M_G_CSF_M2_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_M_G_CSF_M2.tiff",dpi = 300, plot = DimPlot.M_G_CSF_M2.1, width = 6, height = 6)
saveRDS(M_G_CSF_M2_1s, file = "M_G_CSF_M2_1s.rds")

nExp_poi_M_G_CSF_M3_1 <- round(0.08*length(M_G_CSF_M3_1@meta.data$orig.ident)*length(M_G_CSF_M3_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
M_G_CSF_M3_1_DF <- doubletFinder_v3(M_G_CSF_M3_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_M_G_CSF_M3_1, reuse.pANN = FALSE, sct = FALSE)
print(head(M_G_CSF_M3_1_DF@meta.data))
table(M_G_CSF_M3_1_DF$DF.classifications_0.25_0.09_0)

DimPlotM_G_CSF_M3 <- DimPlot(M_G_CSF_M3_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_M_G_CSF_M3.tiff",dpi = 300, plot = DimPlotM_G_CSF_M3, width = 6, height = 6)
M_G_CSF_M3_1s=subset(M_G_CSF_M3_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.M_G_CSF_M3.1 <- DimPlot(M_G_CSF_M3_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_M_G_CSF_M3.tiff",dpi = 300, plot = DimPlot.M_G_CSF_M3.1, width = 6, height = 6)
saveRDS(M_G_CSF_M3_1s, file = "M_G_CSF_M3_1s.rds")

nExp_poi_M_G_CSF_M4_1 <- round(0.08*length(M_G_CSF_M4_1@meta.data$orig.ident)*length(M_G_CSF_M4_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
M_G_CSF_M4_1_DF <- doubletFinder_v3(M_G_CSF_M4_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_M_G_CSF_M4_1, reuse.pANN = FALSE, sct = FALSE)
print(head(M_G_CSF_M4_1_DF@meta.data))
table(M_G_CSF_M4_1_DF$DF.classifications_0.25_0.09_0)

DimPlotM_G_CSF_M4 <- DimPlot(M_G_CSF_M4_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_M_G_CSF_M4.tiff",dpi = 300, plot = DimPlotM_G_CSF_M4, width = 6, height = 6)
M_G_CSF_M4_1s=subset(M_G_CSF_M4_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.M_G_CSF_M4.1 <- DimPlot(M_G_CSF_M4_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_M_G_CSF_M4.tiff",dpi = 300, plot = DimPlot.M_G_CSF_M4.1, width = 6, height = 6)
saveRDS(M_G_CSF_M4_1s, file = "M_G_CSF_M4_1s.rds")

nExp_poi_P_H1_1 <- round(0.08*length(P_H1_1@meta.data$orig.ident)*length(P_H1_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
P_H1_1_DF <- doubletFinder_v3(P_H1_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_P_H1_1, reuse.pANN = FALSE, sct = FALSE)
print(head(P_H1_1_DF@meta.data))
table(P_H1_1_DF$DF.classifications_0.25_0.09_0)

DimPlotP_H1 <- DimPlot(P_H1_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_P_H1.tiff",dpi = 300, plot = DimPlotP_H1, width = 6, height = 6)
P_H1_1s=subset(P_H1_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.P_H1.1 <- DimPlot(P_H1_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_P_H1.tiff",dpi = 300, plot = DimPlot.P_H1.1, width = 6, height = 6)
saveRDS(P_H1_1s, file = "P_H1_1s.rds")

nExp_poi_P_H2_1 <- round(0.08*length(P_H2_1@meta.data$orig.ident)*length(P_H2_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
P_H2_1_DF <- doubletFinder_v3(P_H2_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_P_H2_1, reuse.pANN = FALSE, sct = FALSE)
print(head(P_H2_1_DF@meta.data))
table(P_H2_1_DF$DF.classifications_0.25_0.09_0)

DimPlotP_H2 <- DimPlot(P_H2_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_P_H2.tiff",dpi = 300, plot = DimPlotP_H2, width = 6, height = 6)
P_H2_1s=subset(P_H2_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.P_H2.1 <- DimPlot(P_H2_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_P_H2.tiff",dpi = 300, plot = DimPlot.P_H2.1, width = 6, height = 6)
saveRDS(P_H2_1s, file = "P_H2_1s.rds")

nExp_poi_P_G_CSF_M1_1 <- round(0.08*length(P_G_CSF_M1_1@meta.data$orig.ident)*length(P_G_CSF_M1_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
P_G_CSF_M1_1_DF <- doubletFinder_v3(P_G_CSF_M1_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_P_G_CSF_M1_1, reuse.pANN = FALSE, sct = FALSE)
print(head(P_G_CSF_M1_1_DF@meta.data))
table(P_G_CSF_M1_1_DF$DF.classifications_0.25_0.09_0)

DimPlotP_G_CSF_M1 <- DimPlot(P_G_CSF_M1_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_P_G_CSF_M1.tiff",dpi = 300, plot = DimPlotP_G_CSF_M1, width = 6, height = 6)
P_G_CSF_M1_1s=subset(P_G_CSF_M1_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.P_G_CSF_M1.1 <- DimPlot(P_G_CSF_M1_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_P_G_CSF_M1.tiff",dpi = 300, plot = DimPlot.P_G_CSF_M1.1, width = 6, height = 6)
saveRDS(P_G_CSF_M1_1s, file = "P_G_CSF_M1_1s.rds")

nExp_poi_P_G_CSF_M2_1 <- round(0.08*length(P_G_CSF_M2_1@meta.data$orig.ident)*length(P_G_CSF_M2_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
P_G_CSF_M2_1_DF <- doubletFinder_v3(P_G_CSF_M2_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_P_G_CSF_M2_1, reuse.pANN = FALSE, sct = FALSE)
print(head(P_G_CSF_M2_1_DF@meta.data))
table(P_G_CSF_M2_1_DF$DF.classifications_0.25_0.09_0)

DimPlotP_G_CSF_M2 <- DimPlot(P_G_CSF_M2_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_P_G_CSF_M2.tiff",dpi = 300, plot = DimPlotP_G_CSF_M2, width = 6, height = 6)
P_G_CSF_M2_1s=subset(P_G_CSF_M2_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.P_G_CSF_M2.1 <- DimPlot(P_G_CSF_M2_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_P_G_CSF_M2.tiff",dpi = 300, plot = DimPlot.P_G_CSF_M2.1, width = 6, height = 6)
saveRDS(P_G_CSF_M2_1s, file = "P_G_CSF_M2_1s.rds")

nExp_poi_P_G_CSF_M3_1 <- round(0.08*length(P_G_CSF_M3_1@meta.data$orig.ident)*length(P_G_CSF_M3_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
P_G_CSF_M3_1_DF <- doubletFinder_v3(P_G_CSF_M3_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_P_G_CSF_M3_1, reuse.pANN = FALSE, sct = FALSE)
print(head(P_G_CSF_M3_1_DF@meta.data))
table(P_G_CSF_M3_1_DF$DF.classifications_0.25_0.09_0)

DimPlotP_G_CSF_M3 <- DimPlot(P_G_CSF_M3_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_P_G_CSF_M3.tiff",dpi = 300, plot = DimPlotP_G_CSF_M3, width = 6, height = 6)
P_G_CSF_M3_1s=subset(P_G_CSF_M3_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.P_G_CSF_M3.1 <- DimPlot(P_G_CSF_M3_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_P_G_CSF_M3.tiff",dpi = 300, plot = DimPlot.P_G_CSF_M3.1, width = 6, height = 6)
saveRDS(P_G_CSF_M3_1s, file = "P_G_CSF_M3_1s.rds")

nExp_poi_P_G_CSF_M4_1 <- round(0.08*length(P_G_CSF_M4_1@meta.data$orig.ident)*length(P_G_CSF_M4_1@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
P_G_CSF_M4_1_DF <- doubletFinder_v3(P_G_CSF_M4_1, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi_P_G_CSF_M4_1, reuse.pANN = FALSE, sct = FALSE)
print(head(P_G_CSF_M4_1_DF@meta.data))
table(P_G_CSF_M4_1_DF$DF.classifications_0.25_0.09_0)

DimPlotP_G_CSF_M4 <- DimPlot(P_G_CSF_M4_1_DF,reduction = "umap",group.by = "DF.classifications_0.25_0.09_0")
ggsave("UMAP_pre_double_removal_P_G_CSF_M4.tiff",dpi = 300, plot = DimPlotP_G_CSF_M4, width = 6, height = 6)
P_G_CSF_M4_1s=subset(P_G_CSF_M4_1_DF, subset = DF.classifications_0.25_0.09_0 != "Doublet")
DimPlot.P_G_CSF_M4.1 <- DimPlot(P_G_CSF_M4_1s,reduction = "umap")
ggsave("UMAP_post_double_removal_P_G_CSF_M4.tiff",dpi = 300, plot = DimPlot.P_G_CSF_M4.1, width = 6, height = 6)
saveRDS(P_G_CSF_M4_1s, file = "P_G_CSF_M4_1s.rds")