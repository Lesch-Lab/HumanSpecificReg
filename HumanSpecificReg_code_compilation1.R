#R scripts for Human-specific epigenomic states in spermatogenesis
#Data Source: Lesch, B.J., Silber, S.J., McCarrey, J.R. & Page, D.C. 
#Parallel evolution of male germline epigenetic poising and somatic 
#development in animals. Nat Genet 48, 888-94 (2016).

#Code written by: Caiyun Liao MD MPH, Bluma Lesch MD PhD
#2022-2023

#Part 1: Derive mean of data from the biological replicates, define gene states
#Sort genes into groups based on pattern of consistency/divergence across
#the species/cell type combinations, and plot these distributions

#setwd ("") #specify work directory here

#This program uses the expss package
library(expss)
library(dplyr)
library(pipeR)
library(dbplyr)
library(tidyverse)

rm(list=ls())

#Import original ChIP & RNA Seq data
raw <- read.table("bigTable_160126_clustering.txt",
                  header= TRUE,
                  sep="",
                  row.names = NULL,
                  na.strings = "NA")

anyNA(raw) #No missing data

anyDuplicated(raw[,1]) #No duplicates in the gene names

#Rhesus/ mouse/ opossum each have two replicates
#Get means to use for analysis

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


#Cutoffs used to define chromatin status: 
#According to Lesch et al Nature Genetics 2016
#K4 & K27: 0.5 input-subtracted read/million
#RNA Seq: 5 FPKM

rhesus %>>%
        select(contains(c("_mean","hsapiens_homolog_ensembl_gene","gene_name"))) %>>%
        (~rhesus)

pachy_rh <- rep(0,nrow(rhesus))
rs_rh <- rep(0,nrow(rhesus))                

rhesus_data <- cbind(rhesus,pachy_rh,rs_rh)

#Define status of pachy & rs, each has 8 permutations

#default is 0 (repressed)

colnames(rhesus_data) <- gsub("_mean","",colnames(rhesus_data))

rhesus_data$pachy_rh <- with(rhesus_data,
                             ifelse(pachy_k4 >= 0.5 &
                             pachy_k27 < 0.5 &
                             pachyFPKM >= 5,
                             pachy_rh <- 1,
                             ifelse(pachy_k4 >= 0.5 &
                             pachy_k27 >= 0.5 &
                             pachyFPKM < 5,
                             pachy_rh <- 2,
                             ifelse(pachy_k4 < 0.5 &
                             pachy_k27 < 0.5 &
                             pachyFPKM < 5,
                             pachy_rh <- 3,
                             ifelse(pachy_k4 >= 0.5 &
                             pachy_k27 < 0.5 &
                             pachyFPKM < 5,
                             pachy_rh <- 4,
                             ifelse(pachy_k4 >= 0.5 &
                             pachy_k27 >= 0.5 &
                             pachyFPKM >= 5,
                             pachy_rh <- 5,
                             ifelse(pachy_k4 < 0.5 &
                             pachy_k27 >= 0.5 &
                             pachyFPKM >= 5,
                             pachy_rh <- 6,
                             ifelse(pachy_k4 < 0.5 &
                             pachy_k27<0.5 &                                                
                             pachyFPKM>=5, 
                             pachy_rh<-7, 
                             pachy_rh <- 0))))))))

#round spermatids; default is 0 (repressed)

rhesus_data$rs_rh <- with(rhesus_data,
                          ifelse(rs_k4>=0.5 &
                                 rs_k27<0.5 &
                                 rsFPKM>=5, 
                                 rs_rh<- 1,
                                 ifelse(rs_k4>=0.5 & 
                                 rs_k27>=0.5 & 
                                 rsFPKM< 5,
                                 rs_rh<- 2,
                                 ifelse(rs_k4<0.5 &
                                 rs_k27<0.5 &
                                 rsFPKM<5, 
                                 rs_rh <-3,
                                 ifelse(rs_k4>=0.5 &
                                 rs_k27<0.5 &
                                 rsFPKM<5, 
                                 rs_rh <-4,
                                 ifelse(rs_k4>=0.5 &
                                 rs_k27>=0.5 &
                                 rsFPKM>=5,
                                 rs_rh<-5,
                                 ifelse(rs_k4<0.5 &
                                 rs_k27>=0.5 &
                                 rsFPKM>=5, 
                                 rs_rh<-6,
                                 ifelse(rs_k4<0.5 &
                                 rs_k27<0.5 &
                                 rsFPKM>=5, 
                                 rs_rh <- 7, 
                                 rs_rh <- 0 ))))))))

#give value labels to mark histone modification status
val_lab <- "0 Repressed               
           1 Activated
           2 Poised
           3 All neg
           4 K4 only
           5 All pos
           6 K27 & exp
           7 Exp only"


rhesus_data = apply_labels(rhesus_data,
                           pachy_rh= num_lab(val_lab),
                           rs_rh= num_lab(val_lab)
)

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

pachy_mm <- rep(0,nrow(mouse))
rs_mm <- rep(0,nrow(mouse))                

mouse_data <- cbind(mouse, pachy_mm,rs_mm)

#Define status of pachy & rs, each has 8 permutations

#default is 0 (repressed)

colnames(mouse_data) <- gsub("_mean","",colnames(mouse_data))

mouse_data$pachy_mm <- with(mouse_data,
                            ifelse(pachy_k4 >= 0.5 &
                                   pachy_k27 < 0.5 &
                                   pachyFPKM >= 5,
                                   pachy_mm <- 1,
                                   ifelse(pachy_k4 >= 0.5 &
                                   pachy_k27 >= 0.5 &
                                   pachyFPKM < 5,
                                   pachy_mm <- 2,
                                   ifelse(pachy_k4 < 0.5 &
                                   pachy_k27 < 0.5 &
                                   pachyFPKM < 5,
                                   pachy_mm <- 3,
                                   ifelse(pachy_k4 >= 0.5 &
                                   pachy_k27 < 0.5 &
                                   pachyFPKM < 5,
                                   pachy_mm <- 4,
                                   ifelse(pachy_k4 >= 0.5 &
                                   pachy_k27 >= 0.5 &
                                   pachyFPKM >= 5,
                                   pachy_mm <- 5,
                                   ifelse(pachy_k4 < 0.5 &
                                   pachy_k27 >= 0.5 &
                                   pachyFPKM >= 5,
                                   pachy_mm <- 6,
                                   ifelse(pachy_k4 < 0.5 &
                                   pachy_k27<0.5 &                                                
                                   pachyFPKM>=5, 
                                   pachy_mm<-7, 
                                   pachy_mm <- 0))))))))

#round spermatids; default is 0 (repressed)

mouse_data$rs_mm <- with(mouse_data,
                         ifelse(rs_k4>=0.5 &
                                rs_k27<0.5 &
                                rsFPKM>=5, 
                                rs_mm<- 1,
                                ifelse(rs_k4>=0.5 & 
                                rs_k27>=0.5 & 
                                rsFPKM< 5,
                                rs_mm<- 2,
                                ifelse(rs_k4<0.5 &
                                rs_k27<0.5 &
                                rsFPKM<5, 
                                rs_mm <-3,
                                ifelse(rs_k4>=0.5 &
                                rs_k27<0.5 &
                                rsFPKM<5, 
                                rs_mm <-4,
                                ifelse(rs_k4>=0.5 &
                                rs_k27>=0.5 &
                                rsFPKM>=5,
                                rs_mm<-5,
                                ifelse(rs_k4<0.5 &
                                rs_k27>=0.5 &
                                rsFPKM>=5, 
                                rs_mm<-6,
                                ifelse(rs_k4<0.5 &
                                rs_k27<0.5 &
                                rsFPKM>=5, 
                                rs_mm <- 7, 
                                rs_mm <- 0 ))))))))

#give value labels to mark histone modification status
mouse_data = apply_labels(mouse_data,
                          pachy_mm= num_lab(val_lab),
                          rs_mm= num_lab(val_lab)
)

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

pachy_md <- rep(0,nrow(opossum))
rs_md <- rep(0,nrow(opossum))                

opossum_data <- cbind(opossum,pachy_md,rs_md)

#Define status of pachy & rs, each has 8 permutations

#default is 0 (repressed)

colnames(opossum_data) <- gsub("_mean","",colnames(opossum_data))

opossum_data$pachy_md <- with(opossum_data,
                              ifelse(pachy_k4 >= 0.5 &
                              pachy_k27 < 0.5 &
                              pachyFPKM >= 5,
                              pachy_md <- 1,
                              ifelse(pachy_k4 >= 0.5 &
                              pachy_k27 >= 0.5 &
                              pachyFPKM < 5,
                              pachy_md <- 2,
                              ifelse(pachy_k4 < 0.5 &
                              pachy_k27 < 0.5 &
                              pachyFPKM < 5,
                              pachy_md <- 3,
                              ifelse(pachy_k4 >= 0.5 &
                              pachy_k27 < 0.5 &
                              pachyFPKM < 5,
                              pachy_md <- 4,
                              ifelse(pachy_k4 >= 0.5 &
                              pachy_k27 >= 0.5 &
                              pachyFPKM >= 5,
                              pachy_md <- 5,
                              ifelse(pachy_k4 < 0.5 &
                              pachy_k27 >= 0.5 &
                              pachyFPKM >= 5,
                              pachy_md <- 6,
                              ifelse(pachy_k4 < 0.5 &
                              pachy_k27<0.5 &                                                
                              pachyFPKM>=5, 
                              pachy_md<-7, 
                              pachy_md <- 0))))))))

#round spermatids; default is 0 (repressed)

opossum_data$rs_md <- with(opossum_data,
                           ifelse(rs_k4>=0.5 &
                           rs_k27<0.5 &
                           rsFPKM>=5, 
                           rs_md<- 1,
                           ifelse(rs_k4>=0.5 & 
                           rs_k27>=0.5 & 
                           rsFPKM< 5,
                           rs_md<- 2,
                           ifelse(rs_k4<0.5 &
                           rs_k27<0.5 &
                           rsFPKM<5, 
                           rs_md <-3,
                           ifelse(rs_k4>=0.5 &
                           rs_k27<0.5 &
                           rsFPKM<5, 
                           rs_md <-4,
                           ifelse(rs_k4>=0.5 &
                           rs_k27>=0.5 &
                           rsFPKM>=5,
                           rs_md<-5,
                           ifelse(rs_k4<0.5 &
                           rs_k27>=0.5 &
                           rsFPKM>=5, 
                           rs_md<-6,
                           ifelse(rs_k4<0.5 &
                           rs_k27<0.5 &
                           rsFPKM>=5, 
                           rs_md <- 7, 
                           rs_md <- 0 ))))))))

#give value labels to mark histone modification status
opossum_data = apply_labels(opossum_data,
                            pachy_md= num_lab(val_lab),
                            rs_md= num_lab(val_lab)
)


###HUMANS- 3 REPLICATES, replicates 2&3 were used for analysis. See manuscript for 
#detailed explanations
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

pachy_hs_23 <- rep(0,nrow(human_data_23))
rs_hs_23 <- rep(0,nrow(human_data_23))                

human_data_23 <- cbind(human_data_23,pachy_hs_23,rs_hs_23)

#Define status of pachy & rs, each has 8 permutations

#default is 0 (repressed)

colnames(human_data_23) <- gsub("_mean","",colnames(human_data_23))

human_data_23$pachy_hs_23 <- with(human_data_23,
                                  ifelse(pachy_k4 >= 0.5 &
                                         pachy_k27 < 0.5 &
                                         pachyFPKM >= 5,
                                         pachy_hs_23 <- 1,
                                         ifelse(pachy_k4 >= 0.5 &
                                         pachy_k27 >= 0.5 &
                                         pachyFPKM < 5,
                                         pachy_hs_23 <- 2,
                                         ifelse(pachy_k4 < 0.5 &
                                         pachy_k27 < 0.5 &
                                         pachyFPKM < 5,
                                         pachy_hs_23 <- 3,
                                         ifelse(pachy_k4 >= 0.5 &
                                         pachy_k27 < 0.5 &
                                         pachyFPKM < 5,
                                         pachy_hs_23 <- 4,
                                         ifelse(pachy_k4 >= 0.5 &
                                         pachy_k27 >= 0.5 &
                                         pachyFPKM >= 5,
                                         pachy_hs_23 <- 5,
                                         ifelse(pachy_k4 < 0.5 &
                                         pachy_k27 >= 0.5 &
                                         pachyFPKM >= 5,
                                         pachy_hs_23 <- 6,
                                         ifelse(pachy_k4 < 0.5 &
                                         pachy_k27<0.5 &                                                
                                         pachyFPKM>=5, 
                                         pachy_hs_23<-7, 
                                         pachy_hs_23 <- 0))))))))

#round spermatids; default is 0 (repressed)

human_data_23$rs_hs_23 <- with(human_data_23,
                               ifelse(rs_k4>=0.5 &
                               rs_k27<0.5 &
                               rsFPKM>=5, 
                               rs_hs_23<- 1,
                               ifelse(rs_k4>=0.5 & 
                               rs_k27>=0.5 & 
                               rsFPKM< 5,
                               rs_hs_23<- 2,
                               ifelse(rs_k4<0.5 &
                               rs_k27<0.5 &
                               rsFPKM<5, 
                               rs_hs_23 <-3,
                               ifelse(rs_k4>=0.5 &
                               rs_k27<0.5 &
                               rsFPKM<5, 
                               rs_hs_23 <-4,
                               ifelse(rs_k4>=0.5 &
                               rs_k27>=0.5 &
                               rsFPKM>=5,
                               rs_hs_23<-5,
                               ifelse(rs_k4<0.5 &
                               rs_k27>=0.5 &
                               rsFPKM>=5, 
                               rs_hs_23<-6,
                               ifelse(rs_k4<0.5 &
                               rs_k27<0.5 &
                               rsFPKM>=5, 
                               rs_hs_23 <- 7, 
                               rs_hs_23 <- 0 ))))))))

#give value labels to mark histone modification status

human_data_23 = apply_labels(human_data_23,
                             pachy_hs_23= num_lab(val_lab),
                             rs_hs_23= num_lab(val_lab))

#Below are processing data to classify gene states for each species

#Rhesus
rhesus_data %>>%
        select(contains(c("pachy_rh","rs_rh","gene"))) %>>%
        (~rhesus_pachy_rs)

#Mouse
mouse_data %>>%
        select (contains(c("pachy_mm","rs_mm","gene"))) %>>%
        (~mouse_pachy_rs)

#Opossum
opossum_data %>>%
        select(contains(c("pachy_md","rs_md","gene"))) %>>%
        (~opossum_pachy_rs)

#Human
human_data_23 %>>%
        select(contains(c("pachy_hs_23","rs_hs_23","gene"))) %>>%
        (~human_pachy_rs)

#All four species combined 
list (rhesus_pachy_rs,mouse_pachy_rs,opossum_pachy_rs,
      human_pachy_rs) %>>%
        reduce(full_join,by=c("hsapiens_homolog_ensembl_gene")) %>>%
        (~all_pachy_rs)

saveRDS(all_pachy_rs,file="all_pachy_rs.rds")

all_pachy_rs$all_same <- 0

all_pachy_rs$id <- seq_along(all_pachy_rs[,1])

#Identify genes in the same status in pachy & RS across all four species
#Must convert each row to a list in order for the code to work
for (i in 1:nrow(all_pachy_rs)) {
        
        all_pachy_rs %>>%
                filter(id==i) %>>%
                select(contains(c("pachy_","rs_"))) %>>%
                as.list() %>>%
                (~temp_list)
        
        ifelse( length(unique(temp_list))==1,
                all_pachy_rs[i,"all_same"]<-1,
                all_pachy_rs[i,"all_same"]<-0) 
}

#Sort the data frame based on status and whether status remains the same
all_pachy_rs %>>% 
        filter(all_same==1) %>>%
        tab_sort_asc((c("pachy_rh","rs_rh",
                        "pachy_mm","rs_mm",
                        "pachy_md","rs_md",
                        "pachy_hs_23","rs_hs_23"))) %>>% 
        select(contains(c("pachy_rh","gene","all_same"))) %>>%
        write_csv(file="all_same.csv") 

#Sanity check: compare the statistics from this set to those in the 2016 Lesch paper
#We know there are 405 core poised genes
all_pachy_rs %>>%
        filter(all_same==1) %>>%
        filter(pachy_rh==2) %>>%
        select(contains(c("pachy_rh","gene"))) %>>%
        (~stably_poised)

#In the 2016 Lesch paper, there were 331 genes that were poised in hs/rh/mm/md
print(nrow(stably_poised)) 

#In light of challenges from different reconciliation of human replicates, removing 
#human replicates and repeat sanity checks

#rh-mm-md all poised
all_pachy_rs %>>%
        filter(all_same==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~all_same_3species)

sanity_check<- data.frame(count=1:21,
                          row.names=c("all_same_3species",
                                      "two_states_3species",
                                      "two_states_rh_mm",
                                      "two_states_rh_md",
                                      "two_states_mm_md",
                                      "three_states_3species",
                                      "three_states_rh_mm",
                                      "three_states_rh_md",
                                      "three_states_mm_md",
                                      "four_states_rh_mm",
                                      "four_states_rh_md",
                                      "four_states_mm_md",
                                      "five_states_rh_mm",
                                      "five_states_rh_md",
                                      "five_states_mm_md",
                                      "six_states_rh_mm",
                                      "six_states_rh_md",
                                      "six_states_mm_md",
                                      "seven_states_rh_mm",
                                      "seven_states_rh_md",
                                      "seven_states_mm_md"))

sanity_check ["all_same_3species",] <- nrow(all_same_3species)

#Generate a gene set that excludes states 4-7 and are active only in PS
#or RS of human

subset_not_same <- filter(all_pachy_rs,all_same==0)

hs_activated1 <- subset_not_same %>>%
        subset((pachy_hs_23==1 | rs_hs_23==1) &
                       pachy_rh!=1 & rs_rh!=1 &
                       pachy_mm!=1 & rs_mm!=1 &
                       pachy_md!=1 & rs_md!=1) %>>%
        select(!contains(c("id","all_same")))

hs_activated2 <- hs_activated1 %>>%
        subset(pachy_hs_23 %in% c(0,1,2,3) &
                       rs_hs_23 %in% c(0,1,2,3) &
                       pachy_rh %in% c(0,1,2,3) & 
                       rs_rh %in% c(0,1,2,3) &
                       pachy_mm %in% c(0,1,2,3) &
                       rs_mm %in% c(0,1,2,3) &
                       pachy_md %in% c(0,1,2,3) &
                       rs_md %in% c(0,1,2,3)) %>>%
        select(!contains(c("id","all_same")))

saveRDS(hs_activated2,file= "hs_activated_0to3_only.rds")

#Poised in human PS and/or RS only
hs_poised1 <- subset_not_same %>>%
        subset((pachy_hs_23==2 | rs_hs_23==2) &
                       pachy_rh!=2 & rs_rh!=2 &
                       pachy_mm!=2 & rs_mm!=2 &
                       pachy_md!=2 & rs_md!=2) %>>%
        select(!contains(c("id","all_same")))

hs_poised2 <- hs_poised1 %>>%
        subset(pachy_hs_23 %in% c(0,1,2,3) &
                       rs_hs_23 %in% c(0,1,2,3) &
                       pachy_rh %in% c(0,1,2,3) & 
                       rs_rh %in% c(0,1,2,3) &
                       pachy_mm %in% c(0,1,2,3) &
                       rs_mm %in% c(0,1,2,3) &
                       pachy_md %in% c(0,1,2,3) &
                       rs_md %in% c(0,1,2,3)) %>>%
        select(!contains(c("id","all_same")))

saveRDS(hs_poised2,file= "hs_poised_0to3_only.rds")


#Now, isolate those with two distinct states across 8 columns

subset_not_same <- filter(all_pachy_rs,all_same==0)
subset_not_same$id <- seq_along(subset_not_same[,1])
subset_not_same <- select(subset_not_same,!contains("all_same"))
subset_not_same$two_states <- 0

for (i in 1:nrow(subset_not_same)) {
        
        subset_not_same %>>%
                filter(id==i) %>>%
                select(contains(c("pachy_","rs_"))) %>>%
                as.list() %>>%
                (~temp_list)
        
        ifelse( length(unique(temp_list))==2,
                subset_not_same[i,"two_states"]<-1,
                subset_not_same[i,"two_states"]<-0) 
}

#Sort the data frame based on status and whether there were 2 states across 8 columns
subset_not_same %>>% 
        filter(two_states==1) %>>%
        tab_sort_asc((c("pachy_rh","pachy_mm",
                        "pachy_md","pachy_hs_23",
                        "rs_rh","rs_mm",
                        "rs_md","rs_hs_23"))) %>>% 
        select(!contains("id")) %>>%
        write_csv(file="two_states.csv")

#Sanity check: how many genes were poised in both cell types in rhesus, mouse, and opossum

#poised in rh-mm-md
subset_not_same %>>%
        filter(two_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~two_states_3species)

sanity_check ["two_states_3species",] <- nrow(two_states_3species)

#poised in rh-mm only
subset_not_same %>>%
        filter(two_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        select(contains("gene")) %>>%
        (~two_states_rh_mm)

sanity_check ["two_states_rh_mm",] <- nrow(two_states_rh_mm)

#poised in rh-md only
subset_not_same %>>%
        filter(two_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~two_states_rh_md)

sanity_check ["two_states_rh_md",] <- nrow(two_states_rh_md)

#poised in mm-md only
subset_not_same %>>%
        filter(two_states==1) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~two_states_mm_md)

sanity_check ["two_states_mm_md",] <- nrow(two_states_mm_md)

#Isolate those with 3 distinct statuses across the 8 columns
subset_for_three_states <- filter(subset_not_same,two_states==0)
subset_for_three_states$id <- seq_along(subset_for_three_states[,1])
subset_for_three_states <- select(subset_for_three_states,!contains("two_states"))
subset_for_three_states$three_states <- 0

for (i in 1:nrow(subset_for_three_states)) {
        
        subset_for_three_states %>>%
                filter(id==i) %>>%
                select(contains(c("pachy_","rs_"))) %>>%
                as.list() %>>%
                (~temp_list)
        
        ifelse( length(unique(temp_list))==3,
                subset_for_three_states[i,"three_states"]<-1,
                subset_for_three_states[i,"three_states"]<-0) 
}

#Sort the data frame based on status and whether there were 3 states across 8 columns
subset_for_three_states %>>% 
        filter(three_states==1) %>>%
        tab_sort_asc((c("pachy_rh","pachy_mm",
                        "pachy_md","pachy_hs_23",
                        "rs_rh","rs_mm",
                        "rs_md","rs_hs_23"))) %>>% 
        select(!contains("id")) %>>%
        write_csv(file="three_states.csv")

#Sanity check

#poised in rh-mm-md
subset_for_three_states %>>%
        filter(three_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~three_states_3species)

sanity_check ["three_states_3species",] <- nrow(three_states_3species)

#poised in rh-mm only
subset_for_three_states %>>%
        filter(three_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        select(contains("gene")) %>>%
        (~three_states_rh_mm)

sanity_check ["three_states_rh_mm",] <- nrow(three_states_rh_mm)

#poised in rh-md only
subset_for_three_states %>>%
        filter(three_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~three_states_rh_md)

sanity_check ["three_states_rh_md",] <- nrow(three_states_rh_md)

#poised in mm-md only
subset_for_three_states %>>%
        filter(three_states==1) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~three_states_mm_md)

sanity_check ["three_states_mm_md",] <- nrow(three_states_mm_md)

#Isolate those with 4 distinct statuses across the 8 columns
subset_for_four_states <- filter(subset_for_three_states,three_states==0)

subset_for_four_states$id <- seq_along(subset_for_four_states[,1])
subset_for_four_states <- select(subset_for_four_states,!contains("three_states"))
subset_for_four_states$four_states <- 0

for (i in 1:nrow(subset_for_four_states)) {
        
        subset_for_four_states %>>%
                filter(id==i) %>>%
                select(contains(c("pachy_","rs_"))) %>>%
                as.list() %>>%
                (~temp_list)
        
        ifelse( length(unique(temp_list))==4,
                subset_for_four_states[i,"four_states"]<-1,
                subset_for_four_states[i,"four_states"]<-0) 
}

#Sort the data frame based on status and whether there were 4 states across 8 columns
subset_for_four_states %>>% 
        filter(four_states==1) %>>%
        tab_sort_asc((c("pachy_rh","pachy_mm",
                        "pachy_md","pachy_hs_23",
                        "rs_rh","rs_mm",
                        "rs_md","rs_hs_23"))) %>>% 
        select(!contains("id")) %>>%
        write_csv(file="four_states.csv")

#poised in rh-mm only
subset_for_four_states %>>%
        filter(four_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        select(contains("gene")) %>>%
        (~four_states_rh_mm)

sanity_check ["four_states_rh_mm",] <- nrow(four_states_rh_mm)

#poised in rh-md only
subset_for_four_states %>>%
        filter(four_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~four_states_rh_md)

sanity_check ["four_states_rh_md",] <- nrow(four_states_rh_md)

#poised in mm-md only
subset_for_four_states %>>%
        filter(four_states==1) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~four_states_mm_md)

sanity_check ["four_states_mm_md",] <- nrow(four_states_mm_md)

#Isolate those with 5 distinct statuses across the 8 columns
subset_for_five_states <- filter(subset_for_four_states,four_states==0)
subset_for_five_states$id <- seq_along(subset_for_five_states[,1])
subset_for_five_states <- select(subset_for_five_states,!contains("four_states"))
subset_for_five_states$five_states <- 0

for (i in 1:nrow(subset_for_five_states)) {
        
        subset_for_five_states %>>%
                filter(id==i) %>>%
                select(contains(c("pachy_","rs_"))) %>>%
                as.list() %>>%
                (~temp_list)
        
        ifelse( length(unique(temp_list))==5,
                subset_for_five_states[i,"five_states"]<-1,
                subset_for_five_states[i,"five_states"]<-0) 
}

#Sort the data frame based on status and whether there were 5 states across 8 columns
subset_for_five_states %>>% 
        filter(five_states==1) %>>%
        tab_sort_asc((c("pachy_rh","pachy_mm",
                        "pachy_md","pachy_hs_23",
                        "rs_rh","rs_mm",
                        "rs_md","rs_hs_23"))) %>>% 
        select(!contains("id")) %>>%
        write_csv(file="five_states.csv")

#Sanity check
#poised in rh-mm only
subset_for_five_states %>>%
        filter(five_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        select(contains("gene")) %>>%
        (~five_states_rh_mm)

sanity_check ["five_states_rh_mm",] <- nrow(five_states_rh_mm)

#poised in rh-md only
subset_for_five_states %>>%
        filter(five_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~five_states_rh_md)

sanity_check ["five_states_rh_md",] <- nrow(five_states_rh_md)

#poised in mm-md only
subset_for_five_states %>>%
        filter(five_states==1) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~five_states_mm_md)

sanity_check ["five_states_mm_md",] <- nrow(five_states_mm_md)


#Isolate those with 6 distinct statuses across the 8 columns
subset_for_six_states <- filter(subset_for_five_states,five_states==0)
subset_for_six_states$id <- seq_along(subset_for_six_states[,1])
subset_for_six_states <- select(subset_for_six_states,!contains("five_states"))
subset_for_six_states$six_states <- 0

for (i in 1:nrow(subset_for_six_states)) {
        
        subset_for_six_states %>>%
                filter(id==i) %>>%
                select(contains(c("pachy_","rs_"))) %>>%
                as.list() %>>%
                (~temp_list)
        
        ifelse( length(unique(temp_list))==6,
                subset_for_six_states[i,"six_states"]<-1,
                subset_for_six_states[i,"six_states"]<-0) 
}


#Sort the data frame based on status and whether there were 6 states across 8 columns
subset_for_six_states %>>% 
        filter(six_states==1) %>>%
        tab_sort_asc((c("pachy_rh","pachy_mm",
                        "pachy_md","pachy_hs_23",
                        "rs_rh","rs_mm",
                        "rs_md","rs_hs_23"))) %>>% 
        select(!contains("id")) %>>%
        write_csv(file="six_states.csv")

#In theory there should not be overlap in the six-state group but still running a 
#sanity check here
subset_for_six_states %>>%
        filter(six_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        select(contains("gene")) %>>%
        (~six_states_rh_mm)

sanity_check["six_states_rh_mm",] <- nrow(six_states_rh_mm)

subset_for_six_states %>>%
        filter(six_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~six_states_rh_md)

sanity_check["six_states_rh_md",] <- nrow(six_states_rh_md)

subset_for_six_states %>>%
        filter(six_states==1) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~six_states_mm_md)

sanity_check["six_states_mm_md",] <- nrow(six_states_mm_md)

#Isolate those with 7 distinct statuses across the 8 columns
subset_for_seven_states <- filter(subset_for_six_states,six_states==0)
subset_for_seven_states$id <- seq_along(subset_for_seven_states[,1])
subset_for_seven_states <- select(subset_for_seven_states,!contains("six_states"))
subset_for_seven_states$seven_states <- 0

for (i in 1:nrow(subset_for_seven_states)) {
        
        subset_for_seven_states %>>%
                filter(id==i) %>>%
                select(contains(c("pachy_","rs_"))) %>>%
                as.list() %>>%
                (~temp_list)
        
        ifelse( length(unique(temp_list))==7,
                subset_for_seven_states[i,"seven_states"]<-1,
                subset_for_seven_states[i,"seven_states"]<-0) 
}

#Sort the data frame based on status and whether there were 7 states across 8 columns
subset_for_seven_states %>>% 
        filter(seven_states==1) %>>%
        tab_sort_asc((c("pachy_rh","pachy_mm",
                        "pachy_md","pachy_hs_23",
                        "rs_rh","rs_mm",
                        "rs_md","rs_hs_23"))) %>>% 
        select(!contains("id")) %>>%
        write_csv(file="seven_states.csv")

#Same, still doing a sanity check here
subset_for_seven_states %>>%
        filter(seven_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        select(contains("gene")) %>>%
        (~seven_states_rh_mm)

sanity_check["seven_states_rh_mm",] <- nrow(seven_states_rh_mm)

subset_for_seven_states %>>%
        filter(seven_states==1) %>>%
        filter(pachy_rh==2) %>>%
        filter(rs_rh==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~seven_states_rh_md)

sanity_check["seven_states_rh_md",] <- nrow(seven_states_rh_md)

subset_for_seven_states %>>%
        filter(seven_states==1) %>>%
        filter(pachy_mm==2) %>>%
        filter(rs_mm==2) %>>%
        filter(pachy_md==2) %>>%
        filter(rs_md==2) %>>%
        select(contains("gene")) %>>%
        (~seven_states_mm_md)

sanity_check["seven_states_mm_md",] <- nrow(seven_states_mm_md)

View(sanity_check)

sanity_check[,"group"] <- rownames(sanity_check)

write_excel_csv(sanity_check,"sanity_check.csv")

#Isolate those with 8 distinct statuses across the 8 columns
subset_for_eight_states <- filter(subset_for_seven_states,seven_states==0)
nrow(subset_for_eight_states) #there was no gene that had different states in each 8 columns

