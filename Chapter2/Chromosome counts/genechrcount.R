setwd("~/Documents/Documents/Writing/PhD/July2022/GitHub/Chapter 2 - EDA/Training gene locations")
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)

df <- fread('genes84.csv') 

gene_dist <- fread('hg19Rel92_AllgeneTypes_0kb.txt')

df_dist <- merge(df, gene_dist, by = 'Gene', all.x=TRUE)
df_dist <- df_dist[!duplicated(df_dist$Gene), ]

counts <- df_dist %>%
  group_by(Chromosome) %>%
  summarise(count= n())

png('genes84-chr-count.png', width=2000, height=1000, res=300)
p <- ggplot(data = counts, aes(x = Chromosome, y = count, fill = Chromosome)) +
  geom_bar(stat = "identity") + labs(x = "Chromosome", y = "Possible Gene count") + scale_x_continuous(breaks = seq(1, 22, by = 1)) +
  theme(legend.position = "none")
p
graphics.off()



df <- fread('gene51 chrs.csv') 

gene_dist <- fread('hg19Rel92_AllgeneTypes_0kb.txt')

df_dist <- merge(df, gene_dist, by = 'Gene', all.x=TRUE)
df_dist <- df_dist[!duplicated(df_dist$Gene), ]

counts <- df %>%
  group_by(Chromosome) %>%
  summarise(count= n())

png('genes51-chr-count.png', width=2000, height=1000, res=300)
p <- ggplot(data = counts, aes(x = Chromosome, y = count, fill = Chromosome)) +
  geom_bar(stat = "identity") + labs(x = "Chromosome", y = "Most Likely Gene count") + scale_x_continuous(breaks = seq(1, 22, by = 1)) +
  theme(legend.position = "none")
p
graphics.off()


df <- fread('genes149.csv') 

gene_dist <- fread('hg19Rel92_AllgeneTypes_0kb.txt')

df_dist <- merge(df, gene_dist, by = 'Gene', all.x=TRUE)
df_dist <- df_dist[!duplicated(df_dist$Gene), ]

filt1 <- filter(df_dist, Type == 'protein_coding')

counts <- filt1 %>%
  group_by(Chromosome) %>%
  summarise(count= n())


counts <- df_dist %>%
  group_by(Chromosome) %>%
  summarise(count= n())




png('genes149-chr-count.png', width=2000, height=1000, res=300)
p <- ggplot(data = counts, aes(x = Chromosome, y = count, fill = Chromosome)) +
  geom_bar(stat = "identity") + labs(x = "Chromosome", y = "Probable Gene count") + scale_x_continuous(breaks = seq(1, 22, by = 1)) +
  theme(legend.position = "none")
p
graphics.off()



df <- fread('genes93.csv') 

gene_dist <- fread('hg19Rel92_AllgeneTypes_0kb.txt')

df_dist <- merge(df, gene_dist, by = 'Gene', all.x=TRUE)
df_dist <- df_dist[!duplicated(df_dist$Gene), ]

filt1 <- filter(df_dist, Type == 'pseudogene')

counts <- filt1 %>%
  group_by(Chromosome) %>%
  summarise(count= n())

counts <- df_dist %>%
  group_by(Chromosome) %>%
  summarise(count= n())

png('genes93-chr-count.png', width=2000, height=1000, res=300)
p <- ggplot(data = counts, aes(x = Chromosome, y = count, fill = Chromosome)) +
  geom_bar(stat = "identity") + labs(x = "Chromosome", y = "Least Likely Gene count") + scale_x_continuous(breaks = seq(1, 22, by = 1)) +
  theme(legend.position = "none")
p
graphics.off()