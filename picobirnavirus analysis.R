---
title: "genome analysis of picobirnavirus from hecatomb output"
author: "Binita Febles"
date: '`r format(Sys.Date(), "%B %d, %Y")`
output: html_document
---

```{r global_options, include=FALSE}
knitr::opts_chunk$set(fig.width=20,
                      fig.height=15,
                      fig.path="./figures/",
                      dev='png',
                      warning=FALSE,
                      message=FALSE)
```
## Load Libraries
```{r initiate-environment}
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(data.table)

```
```{r initiate-data}
big_table <-as_tibble(fread(file = 'RC2_IBD_freeze2_viral_bigtable.tsv', header = TRUE))
```

```{r}
#get total number of RNA virus sample numbers
RNA.virus.samples<-big_table%>%
  filter(baltimoreType == 'ssRNA(+)' | baltimoreType == 'ssRNA(-)'| baltimoreType == 'dsRNA'| baltimoreType == 'ssRNA-RT')%>%
  select(sampleID,CPM,family,baltimoreType)%>%
  group_by(family,sampleID)%>%
  summarise(n = n_distinct(sampleID))%>%
  mutate(total_samples = n_distinct(sampleID))%>%
  select(family,total_samples)%>%
  distinct()

#plot total number of RNA virus
ggplot(RNA.virus.samples,aes(x = reorder(family,-total_samples),y=total_samples))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  #geom_text(aes(label=total_samples), hjust=-0.5)+
  scale_y_continuous(expand = c(0,0))+
  labs(y = 'Number of samples',x = 'RNA viral family')+
  theme(axis.text=element_text(size=10))+
  theme(axis.title = element_text(size = 10, face = 'bold'))
```

```{r}
#get filtered total number with CPM > 5

RNA.virus.samples.filtered<-big_table%>%
  filter(baltimoreType == 'ssRNA(+)' | baltimoreType == 'ssRNA(-)'| baltimoreType == 'dsRNA'| baltimoreType == 'ssRNA-RT')%>%
  group_by(sampleID,CPM,family,baltimoreType)%>%
  #summarise(mean.cpm = mean(CPM))%>%
  filter(CPM >= 5)%>%
  group_by(family,sampleID)%>%
  summarise(n = n_distinct(sampleID))%>%
  mutate(total_samples = n_distinct(sampleID))%>%
  select(family,total_samples)%>%
  distinct()

ggplot(RNA.virus.samples.filtered,aes(x = reorder(family,-total_samples),y=total_samples))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  #geom_text(aes(label=total_samples), hjust=-0.5)+
  scale_y_continuous(expand = c(0,0))+
  labs(y = 'Number of samples',x = 'RNA viral family')+
  theme(axis.text=element_text(size=10))+
  theme(axis.title = element_text(size = 10, face = 'bold'))
```

```{r}
#plot alignment of RNA viruses against reference database
RNA.viro.plot<-big_table%>%
  filter(baltimoreType == 'ssRNA(+)' | baltimoreType == 'ssRNA(-)'| baltimoreType == 'dsRNA'| baltimoreType == 'ssRNA-RT', CPM > 5)%>%
  ggplot(aes(x=alnlen,y=pident))+
  geom_jitter(aes(color = alnType),position = position_jitter(h=0.1,w=0.1),alpha=0.6)+
  facet_wrap(~family)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(strip.text = element_text(face = 'bold'))+
  geom_vline(xintercept=150,colour='black',linetype='longdash') +
  geom_hline(yintercept=70,colour='black',linetype='longdash')+
  theme(axis.text=element_text(size=10))+
  theme(axis.title = element_text(size = 10, face = 'bold'))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))+
  labs(x= 'alignment length', y = 'percentage identity')
  
RNA.viro.plot
```
```{r}
# filter bigtable with only 'picobirnaviridae' family
pico_bigtable<-big_table%>%
  filter(family == 'Picobirnaviridae')%>%
  droplevels()

#plot pico raw counts
pico_bigtable%>%
  ggplot(aes(x = sampleID, y = log2(mean(CPM))))+
  geom_bar(stat = 'identity')+
  theme(axis.text.x = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = 'bold'))+
  theme(axis.title = element_text(size = 10, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'))+
  labs(title='Raw Count per sample', x = 'sample', y = 'log2(CPM)')

#filter with count greater than 5
filtered_pico<-pico_bigtable%>%
  filter(count>=5)%>%
  group_by(sampleID,CPM)%>%
  summarise(mean.cpm = mean(CPM))

#plot filtered counts
ggplot(filtered_pico,aes(x=reorder(sampleID,-log2(mean.cpm),sum), y= log2(mean.cpm)))+
  geom_bar(stat = 'identity')+
  theme(axis.text.x=element_blank())+
  scale_y_continuous(expand = c(0,0))+
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = 'bold'))+
  theme(axis.title = element_text(size = 10, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'))+
  labs(title=' Filtered Count per sample', x = 'sample', y = 'log2(CPM)')
```
# contig analysis

```{r}
#load contig annotation table
tax_df <-as_tibble(fread(file = 'contigAnnotations.tsv',stringsAsFactors = TRUE, header = TRUE))


#load contig abundance table
contig_df<-as_tibble(fread(file='contig_count_table.tsv', header = TRUE))

#rename contig id and sample id column from contig abundance table to match ids in contig annotation table
contig_df<-rename(contig_df,contigID=Contig)
contig_df<-rename(contig_df,sampleID= '#Sample')

#join contig annotationa nd contig abundance table
taxon_contig_table<-left_join(contig_df,tax_df, by="contigID")

#get contigs from picobirnaviridae family
pico_contig<-taxon_contig_table%>%
  filter(family == 'Picobirnaviridae')

# get number of total distinct contigs and number of samples in each contig
pico_samples<-pico_contig%>%
  filter(Reads >0)%>%
  distinct(sampleID,contigID,Reads)%>%
  group_by(contigID)%>%
  summarize('total samples'=n(),
            'reads' = sum(Reads))%>%
  ungroup()
```

```{r}

#get samples where contig appears. Filter by reads or CPM values > 0
pico_df<-pico_contig%>%
  filter(Reads > 0)

# add gene column based on blast results
picobirna_final<-pico_df%>%
  mutate(gene = ifelse(
    contigID %in% c('contig_28030','contig_34470','contig_36870','contig_52780','contig_57565','contig_59603'),'Segment1','Segment 2'))

#plot contig and CPM based on gene
ggplot(picobirna_final,aes(x=reorder(contigID,-CPM,sum),y=CPM, fill = gene))+
  geom_bar(stat='identity')+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(title = "Picobirnavirus Contigs",x="contig")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
  theme(axis.text=element_text(size=10))+
  theme(axis.title = element_text(size = 10, face = 'bold'))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))

#plot contig and CPM based on species

ggplot(picobirna_final, aes(x = reorder(contigID,-CPM, sum), y = CPM, fill = species))+
  geom_bar(stat = 'identity')+
  #geom_text(aes(label=pident),position = position_dodge(width = 0.9), vjust = -1)+
  scale_y_continuous(expand = c(0,0))+
  facet_grid(~factor(gene,levels=rev(unique(picobirna_final$gene))),scales = 'free',switch ='x',space = 'free_x')+
  theme(strip.background =element_rect(fill="black"))+
  theme(strip.text = element_text(colour = 'white',face = 'bold'))+
  theme(strip.placement = 'outside')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(title = "Picobirnavirus Contigs on a species level",x="contig")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
  theme(legend.position = 'bottom')

# plot contig abundance by RPKM

sample_contig_RPKM_plot<-ggplot(picobirna_final,aes(x=contigID,y=sampleID, group = gene))+
  geom_tile(aes(fill= log2(RPKM)),colour='lightgrey',size=0.2)+
  geom_vline(xintercept=seq(1.5, length(unique(picobirna_final$contigID))-0.5, 1), lwd=0.05, colour="gray")+
  geom_hline(yintercept=seq(1.5, length(unique(picobirna_final$sampleID))-0.5, 1), lwd=0.05, colour="gray")+
  theme_light() +
  scale_fill_distiller(palette = "YlOrBr", direction = 1)+
  facet_grid(~factor(gene,levels=rev(unique(picobirna_final$gene))),scales = 'free',switch ='x',space = 'free_x')+
  theme(strip.background =element_rect(fill="black"))+
  theme(strip.text = element_text(colour = 'white',face = 'bold'))+
  theme(strip.text.x = element_text(size = 10))+
  theme(strip.placement = 'outside')+
  labs(x=" ", y="sample", title = 'Contig Abundance in Samples', fill = 'Log2(RPKM)')+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
  theme(axis.text=element_text(size=10))+
  theme(axis.title = element_text(size = 10, face = 'bold'))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(position = "bottom",expand = c(0,0))+
  theme(axis.text.y=element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face='bold'))+
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line(color='white'))+
  theme(legend.position = "right", legend.title.align = 0.1,legend.justification = "center")+
  guides(fill=guide_colourbar(barheight = 15, barwidth = 1.5))

sample_contig_RPKM_plot

contig_RPKM_filtered<-picobirna_final%>%
  filter(Reads>=5)%>%
  ggplot(aes(x=contigID,y=sampleID, group = gene))+
  geom_tile(aes(fill= log2(RPKM)),colour='lightgrey',size=0.2)+
  geom_vline(xintercept=seq(1.5, length(unique(picobirna_final$contigID))-0.5, 1), lwd=0.05, colour="gray")+
  geom_hline(yintercept=seq(1.5, length(unique(picobirna_final$sampleID))-0.5, 1), lwd=0.05, colour="gray")+
  theme_light() +
  scale_fill_distiller(palette = "YlOrBr", direction = 1)+
  facet_grid(~factor(gene,levels=rev(unique(picobirna_final$gene))),scales = 'free',switch ='x',space = 'free_x')+
  theme(strip.background =element_rect(fill="black"))+
  theme(strip.text = element_text(colour = 'white',face = 'bold'))+
  theme(strip.text.x = element_text(size = 10))+
  theme(strip.placement = 'outside')+
  labs(x=" ", y="sample", title = 'Contig Abundance in Samples (reads>=5)', fill = 'Log2(RPKM)')+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
  theme(axis.text=element_text(size=10))+
  theme(axis.title = element_text(size = 10, face = 'bold'))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(position = "bottom",expand = c(0,0))+
  theme(axis.text.y=element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face='bold'))+
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_line(color='white'))+
  theme(legend.position = "right", legend.title.align = 0.1,legend.justification = "center")+
  guides(fill=guide_colourbar(barheight = 15, barwidth = 1.5))

contig_RPKM_filtered

# heatmap of samples present 

ggplot(picobirna_final,aes(x = gene, y = sampleID))+
  geom_tile(aes(fill = log2(RPKM)),colour='lightgrey',size=0.2)+
  geom_vline(xintercept=seq(1.5, length(unique(picobirna_final$gene))-0.5, 1), lwd=0.05, colour="gray")+
  theme_classic()+
  scale_fill_distiller(palette = "YlOrBr", direction = 1)+
  theme(axis.text.y=element_blank())+
  theme(axis.text=element_text(size=10))+
  theme(axis.title = element_text(size = 10, face = 'bold'))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))+
  labs(x=" ", y="sample", title = 'Gene Abundance in Samples ', fill = 'Log2(RPKM)')+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
  theme(axis.text.x = element_text(face='bold'))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))

#heatmap of samples after filtered reads
sample_gene_filtered<-picobirna_final%>%
  filter(Reads >=5)%>%
  ggplot(aes(x = gene, y = sampleID))+
  geom_tile(aes(fill = log2(RPKM)),colour='lightgrey',size=0.2)+
  geom_vline(xintercept=seq(1.5, length(unique(picobirna_final$gene))-0.5, 1), lwd=0.05, colour="gray")+
  theme_classic()+
  scale_fill_distiller(palette = "YlOrBr", direction = 1)+
  theme(axis.text.y=element_blank())+
  theme(axis.text.y=element_blank())+
  theme(axis.text=element_text(size=15))+
  theme(axis.title = element_text(size = 15, face = 'bold'))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.title = element_text(size = 15))+
  labs(x=" ", y="sample", title = 'Gene Abundance in Samples (reads>=5) ', fill = 'Log2(RPKM)')+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
  theme(axis.text.x = element_text(face='bold'))+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))

sample_gene_filtered

```


