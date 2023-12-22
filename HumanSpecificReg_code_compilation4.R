#R scripts for Human-specific epigenomic states in spermatogenesis
#Data Source: Lesch, B.J., Silber, S.J., McCarrey, J.R. & Page, D.C. 
#Parallel evolution of male germline epigenetic poising and somatic 
#development in animals. Nat Genet 48, 888-94 (2016).

#Code written by: Caiyun Liao MD MPH, Bluma Lesch MD PhD
#2022-2023

#Part 4: expression bar plots

#setwd ("") #specify work directory here

#This program uses the expss package
library(expss)
library(dplyr)
library(pipeR)
library(dbplyr)
library(ggplot2)
library(RColorBrewer)

rm(list=ls())

raw <- read.table("bigTable_160126_clustering.txt",
                  header= TRUE,
                  sep="",
                  row.names = NULL,
                  na.strings = "NA")

average <- function (x1,x2) {
        
        (x1+x2)/2
        
}

#RHESUS        
raw %>>%
        select(contains(c("rh","hsapiens_homolog_ensembl_gene","gene_name"))) %>>%
        (~rhesus)

rhesus <- within(rhesus,
                 {pachy_k4_mean=average(rh1_pachy_k4,rh2_pachy_k4)
                 pachy_k27_mean=average(rh1_pachy_k27,rh2_pachy_k27)
                 rs_k4_mean=average(rh1_rs_k4,rh2_rs_k4)
                 rs_k27_mean=average(rh1_rs_k27,rh2_rs_k27)
                 pachyFPKM_mean=average(rh1_pachyFPKM,rh2_pachyFPKM)
                 rsFPKM_mean=average(rh1_rsFPKM,rh2_rsFPKM) })

###MOUSE
raw %>>%
        select(contains(c("mm","hsapiens_homolog_ensembl_gene","gene_name"))) %>>%
        (~mouse)

mouse <- within(mouse,
                {pachy_k4_mean=average(mm1_pachy_k4,mm2_pachy_k4)
                pachy_k27_mean=average(mm1_pachy_k27,mm2_pachy_k27)
                rs_k4_mean=average(mm1_rs_k4,mm2_rs_k4)
                rs_k27_mean=average(mm1_rs_k27,mm2_rs_k27)
                pachyFPKM_mean=average(mm1_pachyFPKM,mm2_pachyFPKM)
                rsFPKM_mean=average(mm1_rsFPKM,mm2_rsFPKM) })

mouse %>>%
        select(contains(c("_mean","hsapiens_homolog_ensembl_gene","gene_name"))) %>>%
        (~mouse)

###OPOSSUM
raw %>>%
        select(contains(c("md","hsapiens_homolog_ensembl_gene","gene_name"))) %>>%
        (~opossum)

opossum <- within(opossum,
                  {pachy_k4_mean=average(md1_pachy_k4,md2_pachy_k4)
                  pachy_k27_mean=average(md1_pachy_k27,md2_pachy_k27)
                  rs_k4_mean=average(md1_rs_k4,md2_rs_k4)
                  rs_k27_mean=average(md1_rs_k27,md2_rs_k27)
                  pachyFPKM_mean=average(md1_pachyFPKM,md2_pachyFPKM)
                  rsFPKM_mean=average(md1_rsFPKM,md2_rsFPKM) })

opossum %>>%
        select(contains(c("_mean","hsapiens_homolog_ensembl_gene","gene_name"))) %>>%
        (~opossum)

#Human- excludes replicate 1, only use data from replicates 2 and 3
raw %>>%
        select(contains(c("hs2","hs3","hsapiens_homolog_ensembl_gene","gene_name"))) %>>%
        (~human_data_23)

human_data_23<- within(human_data_23,
                       {pachy_k4_mean=average(hs2_pachy_k4,hs3_pachy_k4)
                       pachy_k27_mean=average(hs2_pachy_k27,hs3_pachy_k27)
                       rs_k4_mean=average(hs2_rs_k4,hs3_rs_k4)
                       rs_k27_mean=average(hs2_rs_k27,hs3_rs_k27)
                       pachyFPKM_mean=average(hs2_pachyFPKM,hs3_pachyFPKM)
                       rsFPKM_mean=average(hs2_rsFPKM,hs3_rsFPKM) })

human_data_23 %>>%
        select(contains(c("_mean","hsapiens_homolog_ensembl_gene","gene_name"))) %>>%
        (~human_data_23)

library(tidyverse)

human_data_23 %>>%
        select(contains(c("gene_name", "hsapiens_homolog_ensembl_gene","pachyFPKM","rsFPKM"))) %>>%
        (~hs_exp)

colnames(hs_exp)[colnames(hs_exp)=="pachyFPKM"] <- "ps_hs"
colnames(hs_exp)[colnames(hs_exp)=="rsFPKM"] <-"rs_hs"

#Rhesus
rhesus %>>%
        select(contains(c("gene_name","hsapiens_homolog_ensembl_gene","pachyFPKM","rsFPKM"))) %>>%
        (~rh_exp)

colnames(rh_exp)[colnames(rh_exp)=="pachyFPKM"] <- "ps_rh"
colnames(rh_exp)[colnames(rh_exp)=="rsFPKM"] <- "rs_rh"

#Mouse
mouse %>>%
        select(contains(c("gene_name","hsapiens_homolog_ensembl_gene","pachyFPKM","rsFPKM"))) %>>%
        (~mm_exp)

colnames(mm_exp)[colnames(mm_exp)=="pachyFPKM"]<- "ps_mm"
colnames(mm_exp)[colnames(mm_exp)=="rsFPKM"] <- "rs_mm"

#Opossum
opossum %>>%
        select(contains(c("gene_name","hsapiens_homolog_ensembl_gene","pachyFPKM","rsFPKM"))) %>>%
        (~md_exp)

colnames(md_exp)[colnames(md_exp)=="pachyFPKM"]<- "ps_md"
colnames(md_exp)[colnames(md_exp)=="rsFPKM"] <-"rs_md"

#Combine these into one single file for the expression data
list (hs_exp,rh_exp,mm_exp,md_exp) %>>%
        reduce(full_join,by=c("hsapiens_homolog_ensembl_gene","gene_name")) %>>%
        (~all_exp)

saveRDS(all_exp,file="all_exp.rds")

#Use the expression figure of BVES as an example
bves <- all_exp[all_exp$gene_name=="BVES",]

bves %>>% select(starts_with(c("ps","rs")))

long <- bves %>>%
        pivot_longer(col=contains(c("ps","rs")),
                     names_to= c("cell_type","species"),
                     names_sep = "_")

long$species<- factor(long$species,
                      level=c("hs",
                              "rh",
                              "mm",
                              "md"))

#Draw bar graph for expression data
color <- brewer.pal(n=1, name="Greys")[8]

bves<-ggplot(data = long) +
        geom_bar(aes(x = cell_type, y = value), stat = "identity", width = 0.3) +
        facet_grid(.~species)+
        labs(title="BVES",
             x="cell types and species",y="expression level")+
        scale_x_discrete(labels=c("PS","RS"))+
        theme_classic()+ 
        theme(plot.title=element_text(hjust = 0.5,size=10)) + 
        scale_fill_manual(values = color) 

bves

#Code for making the collage of genes with human-specific active state

library(tidyr)
library(gridExtra)

raw %>>%
        select(contains(c("gene_name","hs2_pachyFPKM","hs2_rsFPKM",
                          "hs3_pachyFPKM","hs3_rsFPKM",
                          "rh1_pachyFPKM","rh1_rsFPKM",
                          "rh2_pachyFPKM","rh2_rsFPKM",
                          "mm1_pachyFPKM","mm1_rsFPKM",
                          "mm2_pachyFPKM","mm2_rsFPKM",
                          "md1_pachyFPKM","md1_rsFPKM",
                          "md2_pachyFPKM","md2_rsFPKM"))) %>>%
        (~raw_exp)

bar_graph_exp <- merge(raw_exp,all_exp,by="gene_name",all= FALSE)

#Filter out the genes that are uniquely "active" in human
#NEK5 was manually added on 12/21/23 as this was included in the human-specific
#active list based on recent literature reports
bar_graph_exp_active <- bar_graph_exp[bar_graph_exp$gene_name %in%
                                              c("AGMO","ANKFN1","ATRNL1","BVES","FER1L6","FRS2","IGLL1",
                                                "ITGB6","KRT23","LRFN4","MYBPC1","MYO3A","NSG2","OPN5",
                                                "OR2H1","POPDC3","SCN3A","SLC24A1","SNAP91","TNFAIP6",
                                                "XIRP2","ZNF385B","NEK5"),]

#Switch to long format for graphing
bar_graph_exp_active %>>% select(starts_with(c("gene_name","ps","rs","hs2","hs3",
                                               "rh1","rh2","mm1","mm2",
                                               "md1","md2"))) %>>%
        ~bar_graph_exp_active2

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="ps_hs"]<-"ps_hs_mean"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="ps_rh"]<-"ps_rh_mean"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="ps_mm"]<-"ps_mm_mean"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="ps_md"]<-"ps_md_mean"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="rs_hs"]<-"rs_hs_mean"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="rs_rh"]<-"rs_rh_mean"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="rs_mm"]<-"rs_mm_mean"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="rs_md"]<-"rs_md_mean"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="hs2_pachyFPKM"]<-"ps_hs_rep1"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="hs2_rsFPKM"]<-"rs_hs_rep1"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="hs3_pachyFPKM"]<-"ps_hs_rep2"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="hs3_rsFPKM"]<-"rs_hs_rep2"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="rh1_pachyFPKM"]<-"ps_rh_rep1"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="rh1_rsFPKM"]<-"rs_rh_rep1"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="rh2_pachyFPKM"]<-"ps_rh_rep2"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="rh2_rsFPKM"]<-"rs_rh_rep2"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="mm1_pachyFPKM"]<-"ps_mm_rep1"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="mm1_rsFPKM"]<-"rs_mm_rep1"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="mm2_pachyFPKM"]<-"ps_mm_rep2"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="mm2_rsFPKM"]<-"rs_mm_rep2"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="md1_pachyFPKM"]<-"ps_md_rep1"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="md1_rsFPKM"]<-"rs_md_rep1"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="md2_pachyFPKM"]<-"ps_md_rep2"

names(bar_graph_exp_active2)[names(bar_graph_exp_active2)=="md2_rsFPKM"]<-"rs_md_rep2"

long <- bar_graph_exp_active2 %>>%
        pivot_longer(col=contains(c("ps","rs")),
                     names_to= c("cell_type","species","mean_rep"),
                     names_sep = "_")

long$species<- factor(long$species,
                      level=c("hs","rh","mm","md"))

long$mean_rep <- factor(long$mean_rep,
                        level=c("mean","rep1","rep2"))

#Create a function to draw the bar/ scattered plot for each gene
plot_gene <- function(df, gene_id) {
        
        ggplot(data=df[df$gene_name==gene_id,]) +
                geom_bar(data=filter(df,mean_rep=="mean"),stat= "identity",aes(x=cell_type,y=value),width=0.3) +
                geom_point(data=filter(df,mean_rep %in% c("rep1","rep2")),
                           aes(x=cell_type,y=value),
                           color= "black",size=2) +
                facet_grid(.~species)+        
                labs(title=gene_id,x="Cell type/Species",y="Expression Value")+
                scale_x_discrete(labels=c("PS","RS"))+
                theme_classic()+ 
                theme(plot.title=element_text(hjust = 0.5,size=10)) 
}

#Create a list of plots for each gene
gene_plots <- lapply(unique(long$gene_name),function(gene) {
        
        plot_gene(long[long$gene_name==gene,],gene)
})

#Arrange and display the plots in a grid
active_grid<-grid.arrange(grobs=gene_plots,ncol=3)
