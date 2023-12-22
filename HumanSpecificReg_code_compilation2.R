#R scripts for Human-specific epigenomic states in spermatogenesis
#Data Source: Lesch, B.J., Silber, S.J., McCarrey, J.R. & Page, D.C. 
#Parallel evolution of male germline epigenetic poising and somatic 
#development in animals. Nat Genet 48, 888-94 (2016).

#Code written by: Caiyun Liao MD MPH, Bluma Lesch MD PhD
#2022-2023

#Part 2: "rainbow graph" showing epigenomic states of various genes across 
#species/cell type combinations

#setwd ("") #specify work directory here

#This program uses the expss package
library(expss)
library(dplyr)
library(pipeR)
library(dbplyr)
library(tidyverse)

rm(list=ls())

#Plot the number of states in different cell type/species combinations
#All genes are included in the figure below; corresponds to Figure 2a

pooled_data<- readRDS("all_pachy_rs.rds") 

View(pooled_data)

long <- pooled_data %>>%
        pivot_longer(cols=contains(c("rh","mm","md","hs_23")),
                     names_to= c("cell_type","species"),
                     names_sep="_")

View(long)

data_count <- long %>>%
        group_by(cell_type, species, value) %>>%
        summarise(n=n())

View(data_count)

data_count$species[data_count$species=="hs"]<-"Human"
data_count$species[data_count$species=="rh"]<-"Rhesus"
data_count$species[data_count$species=="mm"]<-"Mouse"
data_count$species[data_count$species=="md"]<-"Opossum"

data_count$species<- factor(data_count$species,
                            level=c("Human",
                                    "Rhesus",
                                    "Mouse",
                                    "Opossum"))

#Plotting
colors <- brewer.pal(n=8, name="Set2")[c(7,2,1,8,5,4,3,6)]

data_count$value <- factor(data_count$value)

state_spread <-
        ggplot(data_count, aes(x = cell_type, y = n, fill = value)) + 
        geom_bar(stat = "identity") +
        facet_grid(.~species) +
        labs(title="Status Distribution by Species and Cell Types",
             x="Cell Types and Species",y="Counts of Specific Status")+
        scale_x_discrete(labels=c("PS","RS"))+
        theme_classic()+ 
        theme(plot.title=element_text(hjust = 0.5,size=10))+
        scale_fill_manual(values=colors)


#Distribution of states among genes that have only one state across all 8 columns
#This corresponds to Fig. 2b
all_same <- readRDS("all_same.rds")

long_all_same <- all_same %>>%
        pivot_longer(cols=contains(c("rh","mm","md","hs_23")),
                     names_to= c("cell_type","species"),
                     names_sep="_")

View(long_all_same)

data_count <- long_all_same %>>%
        group_by(cell_type, species, value) %>>%
        summarise(n=n())

View(data_count)

data_count$species[data_count$species=="hs"]<-"Human"
data_count$species[data_count$species=="rh"]<-"Rhesus"
data_count$species[data_count$species=="mm"]<-"Mouse"
data_count$species[data_count$species=="md"]<-"Opossum"

data_count$species<- factor(data_count$species,
                            level=c("Human",
                                    "Rhesus",
                                    "Mouse",
                                    "Opossum"))

#Plotting
colors <- brewer.pal(n=8, name="Set2")[c(7,2,1,8,5,4,3,6)]

data_count$value <- factor(data_count$value)

state_spread_1 <-
        ggplot(data_count, aes(x = cell_type, y = n, fill = value)) + 
        geom_bar(stat = "identity") +
        facet_grid(.~species) +
        labs(title="Uniform state in all species and cell types",
             x="Cell Types and Species",y="Counts of Specific Status")+
        scale_x_discrete(labels=c("PS","RS"))+
        theme_classic()+ 
        theme(plot.title=element_text(hjust = 0.5,size=10))+
        scale_fill_manual(values=colors)

state_spread_1

#Genes that have two states across 8 columns; corresponds to Fig. 2c
two_states <- readRDS("two_states.rds")

long_two_states <- two_states %>>%
        pivot_longer(cols=contains(c("rh","mm","md","hs_23")),
                     names_to= c("cell_type","species"),
                     names_sep="_")

data_count <- long_two_states %>>%
        group_by(cell_type, species, value) %>>%
        summarise(n=n())

data_count$species[data_count$species=="hs"]<-"Human"
data_count$species[data_count$species=="rh"]<-"Rhesus"
data_count$species[data_count$species=="mm"]<-"Mouse"
data_count$species[data_count$species=="md"]<-"Opossum"

data_count$species<- factor(data_count$species,
                            level=c("Human",
                                    "Rhesus",
                                    "Mouse",
                                    "Opossum"))

#Plotting
colors <- brewer.pal(n=8, name="Set2")[c(7,2,1,8,5,4,3,6)]

data_count$value <- factor(data_count$value)

state_spread_2 <-
        ggplot(data_count, aes(x = cell_type, y = n, fill = value)) + 
        geom_bar(stat = "identity") +
        facet_grid(.~species) +
        labs(title="Two states in all species and cell types",
             x="Cell Types and Species",y="Counts of Specific Status")+
        scale_x_discrete(labels=c("PS","RS"))+
        theme_classic()+ 
        theme(plot.title=element_text(hjust = 0.5,size=10))+
        scale_fill_manual(values=colors)

state_spread_2

#Number of genes that have 1-7 states across 8 columns; 
#Corresponds to Fig. 1d

three_states<- readRDS("three_states.rds")
four_states<- readRDS("four_states.rds")
five_states<- readRDS("five_states.rds")
six_states<- readRDS("six_states.rds")
seven_states<- readRDS("seven_states.rds")

state_count <- data_frame(n_state=ordered(c(1:8)),
                          count=0)

state_count$count[state_count$n_state==1] <- nrow(all_same)
state_count$count[state_count$n_state==2] <- nrow(two_states)
state_count$count[state_count$n_state==3] <- nrow(three_states)
state_count$count[state_count$n_state==4] <- nrow(four_states)
state_count$count[state_count$n_state==5] <- nrow(five_states)
state_count$count[state_count$n_state==6] <- nrow(six_states)
state_count$count[state_count$n_state==7] <- nrow(seven_states)

#There was 0 genes with 8 distinct states across the board
state_count$count[state_count$n_state==8] <- 0

state_count_distribution <-
        ggplot(state_count, aes(x = n_state, y = count)) + 
        geom_bar(stat = "identity", fill= "midnightblue") +
        labs(title="Number of genes with distinct states in cell types and species",
             x="Number of States",y="Number of Genes")+
        scale_x_discrete(labels=state_count$n_state)+
        scale_y_continuous(breaks = seq(0,6000,500),labels=seq(0,6000,500)) +
        theme(plot.title=element_text(hjust = 0.5,size=10),
              panel.grid = element_blank(),
              plot.background = element_rect(fill="white"),
              panel.background = element_rect(fill="white"))

state_count_distribution

#Rainbow graph showing chromatin state distribution across 4 species * 2 cell types
#Corresponds to Fig 3a

library(pipeR)
library(tidyverse)
library(data.table)
library(reshape)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

rm(list=ls())

temp <- readRDS("hs_activated_0to3_only.rds") %>>%
        select(-c("gene_name.y","gene_name.x.x","gene_name.y.y"))

temp[,c("pachy_rh","rs_rh",
        "pachy_mm","rs_mm",                        
        "pachy_md","rs_md",
        "pachy_hs_23","rs_hs_23")] <- 
        lapply(temp[,c("pachy_rh","rs_rh",
                       "pachy_mm","rs_mm",                        
                       "pachy_md","rs_md",
                       "pachy_hs_23","rs_hs_23")],as.numeric)

colnames(temp)[colnames(temp)=="pachy_rh"]<- "ps_rh"

colnames(temp)[colnames(temp)=="pachy_mm"]<- "ps_mm"

colnames(temp)[colnames(temp)=="pachy_md"]<- "ps_md"

colnames(temp)[colnames(temp)=="pachy_hs_23"]<- "ps_hs"

colnames(temp)[colnames(temp)=="rs_hs_23"]<- "rs_hs"

colnames(temp)[colnames(temp)=="hsapiens_homolog_ensembl_gene"]<- "gene" 

colnames(temp)[colnames(temp)=="gene_name.x"]<-"gene_name"

temp <- temp %>>%
        subset(select=c("gene","gene_name","ps_hs","rs_hs","ps_rh","rs_rh",
                        "ps_mm","rs_mm","ps_md","rs_md")) %>>%
        arrange()

temp<- arrange(temp,ps_hs,rs_hs,
               ps_rh,rs_rh,
               ps_mm,rs_mm,
               ps_md,rs_md)

temp$order<- seq_along(temp$gene)  

#In order to know which genes is where on the figure, need to generate 
#a separate file linking the "order" variable with gene_name

write.csv(temp,file="hs_specific_activation_rainbow_graph_gene_names.csv")

temp <- temp %>>%
        subset(select = c("order","ps_hs","rs_hs",
                          "ps_rh","rs_rh",
                          "ps_mm","rs_mm",
                          "ps_md","rs_md")) %>>%
        arrange() %>>%
        melt(id="order")

temp <- setnames(temp, c("gene","variable","value"))

#Generate the rainbow plot
colors <- brewer.pal(n=8, name="Set2")[c(7,2,1,8,5,4,3,6)] #focusing on 4 states in this subset

rainbow_plot6 <-
        ggplot(temp, 
               aes(x = variable, y = gene, fill = factor(value))) + 
        geom_tile() + 
        scale_fill_manual(values=colors) +
        labs(fill="Chromatin state") +
        theme(axis.text.y = element_text(hjust = 0, size=6),
              axis.text.x = element_text(size=6),
              axis.title = element_text(size = 8),
              legend.key.width = unit(0.2, "cm"),
              legend.title = element_text(size = 8),
              legend.margin = margin(0,0,0,0),
              plot.background = element_rect(fill="white"),
              panel.background = element_rect(fill="white"))

wrapped_title <- str_wrap("Chromatin state", width = 9)

rainbow_plot6<- rainbow_plot6 + labs(fill=wrapped_title,
                                     x="Cell type and species",
                                     y="Gene") 

#Genes that are poised in human PS and/or RS only
#Corresponds to Fig 4a

temp <- readRDS("hs_poised_0to3_only.rds") %>>%
        select(-c("gene_name.y","gene_name.x.x","gene_name.y.y"))

temp[,c("pachy_rh","rs_rh",
        "pachy_mm","rs_mm",                        
        "pachy_md","rs_md",
        "pachy_hs_23","rs_hs_23")] <- 
        lapply(temp[,c("pachy_rh","rs_rh",
                       "pachy_mm","rs_mm",                        
                       "pachy_md","rs_md",
                       "pachy_hs_23","rs_hs_23")],as.numeric)

colnames(temp)[colnames(temp)=="pachy_rh"]<- "ps_rh"

colnames(temp)[colnames(temp)=="pachy_mm"]<- "ps_mm"

colnames(temp)[colnames(temp)=="pachy_md"]<- "ps_md"

colnames(temp)[colnames(temp)=="pachy_hs_23"]<- "ps_hs"

colnames(temp)[colnames(temp)=="rs_hs_23"]<- "rs_hs"

colnames(temp)[colnames(temp)=="hsapiens_homolog_ensembl_gene"]<- "gene" 

colnames(temp)[colnames(temp)=="gene_name.x"]<-"gene_name"

temp <- temp %>>%
        subset(select=c("gene","gene_name","ps_hs","rs_hs","ps_rh","rs_rh",
                        "ps_mm","rs_mm","ps_md","rs_md")) %>>%
        arrange()

temp<- arrange(temp,ps_hs,rs_hs,
               ps_rh,rs_rh,
               ps_mm,rs_mm,
               ps_md,rs_md)

temp$order<- seq_along(temp$gene)      

temp <- temp %>>%
        subset(select = c("order","ps_hs","rs_hs",
                          "ps_rh","rs_rh",
                          "ps_mm","rs_mm",
                          "ps_md","rs_md")) %>>%
        arrange() %>>%
        melt(id="order")

temp <- setnames(temp, c("gene","variable","value"))

#Generate the rainbow plot
colors <- brewer.pal(n=8, name="Set2")[c(7,2,1,8,5,4,3,6)] #focusing on 4 states in this subset

rainbow_plot8 <-
        ggplot(temp, 
               aes(x = variable, y = gene, fill = factor(value))) + 
        geom_tile() + 
        scale_fill_manual(values=colors) +
        labs(fill="Chromatin state") +
        theme(axis.text.y = element_text(hjust = 0, size=6),
              axis.text.x = element_text(size=6),
              axis.title = element_text(size = 8),
              legend.key.width = unit(0.2, "cm"),
              legend.title = element_text(size = 8),
              legend.margin = margin(0,0,0,0),
              plot.background = element_rect(fill="white"),
              panel.background = element_rect(fill="white"))

wrapped_title <- str_wrap("Chromatin state", width = 9)

rainbow_plot8<- rainbow_plot8 + labs(fill=wrapped_title,
                                     x="Cell type and species",
                                     y="Gene") 