#R scripts for Human-specific epigenomic states in spermatogenesis
#Data Source: Lesch, B.J., Silber, S.J., McCarrey, J.R. & Page, D.C. 
#Parallel evolution of male germline epigenetic poising and somatic 
#development in animals. Nat Genet 48, 888-94 (2016).

#Code written by: Caiyun Liao MD MPH, Bluma Lesch MD PhD
#2022-2023

#Part 3: GO analysis and plotting

#setwd ("") #specify work directory here

#GO analysis for genes with human-specific active epigenetic states
#Codes for genes with human-specific poising were similar

# open required packages and databases
library("data.table")
library("biomaRt")

genes <- read.table("231211_hs_active_with_sperm.txt", header=T, sep="\t")
colnames()
genes <- genes$hsapiens_homolog_ensembl_gene

# get entrez names for all genes
# first get the ensembl Mart from biomaRt
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                   #host="grch37.ensembl.org",
                   #host="http://aug2017.archive.ensembl.org",
                   path="/biomart/martservice",
                   dataset="hsapiens_gene_ensembl",
                   #ensemblRedirect=FALSE
)

genes.entrez <- getBM(attributes = "entrezgene_id",  #updated from entrezgene on 6/9
                      filters = "ensembl_gene_id",
                      values = genes,
                      mart = ensembl)

library("GOstats")
library("org.Hs.eg.db")

# run GOstats to get enriched go terms 
# generate mapping of GO categories to genes 
frame <- toTable(org.Hs.egGO)
goframeData <- data.frame(frame$go_id, frame$Evidence, frame$gene_id)

goFrame <- GOFrame(goframeData, organism="Homo sapiens")
goAllFrame=GOAllFrame(goFrame)
# comment this out if using a custom universe
universe <- Lkeys(org.Hs.egGO)

# set gene sets to test for significant enrichment
testGenes <- genes.entrez$entrezgene
testGenes.ch <- as.character(testGenes)

# make a GOHyperGParams object for each gene set, using all species1 genes as universe, looking for over-representation (test direction)
# to correct for child dependencies, use conditional=TRUE
# set pvalueCutoff=1 to include all tested categories in final output list

params1 <- new("GOHyperGParams", geneIds=testGenes.ch, universeGeneIds=universe, annotation="org.Hs.eg.db", ontology="BP", pvalueCutoff=1, conditional=TRUE, testDirection="over")
OverCond <- hyperGTest(params1)

# Correct for multiple testing using Benjamini Hochberg approach
padj1 <- p.adjust(summary(OverCond)$Pvalue, method="BH", n=length(summary(OverCond)$GOBPID))
summaryOverCond <- summary(OverCond)
summaryOverCond$padj <- padj1

write.table(summaryOverCond, file="231218_GOstats_human_active.txt", sep="\t", quote=F, row.names=F, col.names=T)

#Generate bar plots demonstrating results of GO analysis
#Use genes that are consistently active in all species and cell types as an example

library(pipeR)
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)

#The .csv file used in this code includes the following columns:
#GOBPID, Pvalue, OddsRatio, ExpCount, Count, Size, Term, padj

GO_all_1 <- read.csv("231218_GOstats_all_active.csv", header = TRUE) %>>%
        select(contains(c("GOBPID","Term","Padj"))) 

colnames(GO_all_1)[colnames(GO_all_1) =="padj"] <- "adj_pval"

GO_all_1$adj_pval<- as.numeric(GO_all_1$adj_pval)

# Sort data by adj_pval in descending order
GO_all_1_pselect<- GO_all_1[GO_all_1$adj_pval<0.05,] #define your p value cutoff 
#for illustration on the graph

GO_all_1_pselect <- GO_all_1_pselect[order(GO_all_1_pselect$adj_pval),]

color <- brewer.pal(n=1, name="Set3")[1]

# Generate the plot
GO_all_1_bar <-
        ggplot(GO_all_1_pselect, aes(x = -log10(adj_pval), y = reorder(Term,-adj_pval))) + 
        geom_bar(stat = "identity",fill=color) +
        labs(title = "GO Enrichment Analysis, genes consistently activated",
             x = "-log10(Adjusted P-value)",
             y = "GO Term")+
        theme(axis.text.y = element_text(size = 15),
              axis.text.x = element_text(size=6),
              axis.title = element_text(size = 8))+
        theme_classic()

GO_all_1_bar

