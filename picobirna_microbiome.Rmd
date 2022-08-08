---
title: "microbiome analysis of Picobirnavirus positive samples"
author: "Binita Febles"
date: '`r format(Sys.Date(), "%B %d, %Y")`
output: html_document
---

```{r global_options, include=FALSE}
knitr::opts_chunk$set(fig.width=15,
                      fig.height=10,
                      fig.path="./figures/",
                      dev='png',
                      warning=FALSE,
                      message=FALSE)
```
## Load Libraries
```{r initiate-environment}
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(plotly)
library(data.table)

```

```{r initiate-data}
#load phyloseq object of all the samples
ps0<-readRDS('ps0.M649.M713.M749.park.only.rdp.RDS')
sample_variables(ps0)
sample_names(ps0)
rank_names(ps0)
summarize_phyloseq(ps0)
ntaxa(ps0)
```
#check phylum with unique OTUs if there is any NA. 
```{r}
ps0%>%
  psmelt%>%
  group_by(Phylum)%>%
  summarise(OTUn = (unique(OTU)%>%length))%>%
  arrange(desc(OTUn))

#remove features with 'NA' phylum annotation
ps0 = subset_taxa(ps0, !is.na(Phylum) & !Phylum %in% c('NA'))
get_taxa_unique(ps0,'Phylum')
```
#load picobirnavirus contig file with RPKM abundance and subset phyloseq object with picobirna positive samples
```{r}
p.df<-fread('picobirna_df_final.tsv',stringsAsFactors = TRUE,header = TRUE)

# rename column in pico_contigs
p.df1<-rename(p.df, seqID_virome=sampleID)%>%
  select(seqID_virome,contigID,RPKM,gene)%>%
  mutate(pico_RA=RPKM/sum(RPKM))

#subset ps0 only with pico samples
ps1<-subset_samples(ps0,seqID_virome%in%p.df1$seqID_virome)
ps1

#keep taxa only part of pico positive samples
ntaxa(ps1)
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
```
#Taxon cleaning
```{r}
#remove sequences classified as either mitochondria or chloroplast
ntaxa(ps1) # Check the number of taxa prior to removal
#remove taxa
ps1<-ps1%>%
  subset_taxa(Family  != "mitochondria" & Class   != "Chloroplast")
ntaxa(ps1) # Confirm that the taxa were removed
```
#calculate compositional version of the data
```{r}
#relative abundance
ps1.rel<-microbiome::transform(ps1,'compositional')
#Output of dada2 will most likely have seqs as rownames instead of OTU ids or taxa names
taxa_names(ps1.rel)[1:10]
#change it to ASVIDs
ps1.rel<-microbiome::add_refseq(ps1.rel)
ps1.rel
taxa_names(ps1.rel)
```
#community composition plotting
```{r}
#Make it relative abundance
p.phy<-ps1.rel%>%
  aggregate_taxa(level = 'Phylum')

p.phy<-plot_composition(ps2)+
  scale_fill_brewer("Phylum", palette = "Paired") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_continuous(expand=c(0,0))+
  labs(x='Samples',y='Relative Abundance')+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))
  
p.phy   
ggsave('microbiome_pico.png',width=15,dpi=500)

#select more colors 
cols<-60
mycolors<-colorRampPalette(brewer.pal(8, "Paired"))(cols)

# Family level plotting
p.fam<-ps1.rel%>%
  aggregate_taxa(level = 'Family')

p.fam<-plot_composition(ps.fam)+
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_manual('Family',values=mycolors)+
  scale_y_continuous(expand=c(0,0))+
  labs(x='Samples',y='Relative Abundance')+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))+
  theme(axis.text=element_text(size=10))+
  theme(axis.title = element_text(size = 10, face = 'bold'))+
  theme(legend.text = element_text(size = 8))+
  theme(legend.title = element_text(size = 8))
p.fam
ggsave('microbiome_pico_fam.png',width=20,height=15,dpi=500)  

```
# Core microbiota Analysis.We only need taxa with high prevalence for further analysis
```{r}
# this filters on a family level
# taxa without any filter on abundance and prevalence
nocore.ps1<-core_members(ps1.rel,detection = 0, prevalence = 0 )
nocore.ps1

#detection = 0.01 #filtering out low abundance taxa
#prevalence = atleast 50% prevalent in all samples
core.ps1<-core_members(ps1.rel,detection = 0.01, prevalence = 50/100 )
core.ps1 [1:10]
```
#add taxa names along with ASV ids 
```{r}
# Extract the taxonomy table

taxonomy <- as.data.frame(tax_table(ps1.rel))

# Subset this taxonomy table to include only core OTUs  
core_taxa_id <- subset(taxonomy, rownames(taxonomy) %in% core.ps1)

# also get taxonomy table to include non core OTUs  
non_core_taxa_id <- subset(taxonomy, rownames(taxonomy) %in% nocore.ps1)


DT::datatable(core_taxa_id)
DT::datatable(non_core_taxa_id)
```
#Use the format_to_besthit function from microbiomeutilities to get taxonomic identities of the ASVs.
```{r}
ps1.rel.f <- microbiomeutilities::format_to_besthit(ps1.rel)
taxa_names(ps1.rel.f)
```
#add the best taxonomic classification available to the core and non core microbiota
```{r}
nocore.ps1<-core_members(ps1.rel.f,detection = 0,prevalence = 0)
nocore.ps1

core.ps1<-core_members(ps1.rel.f,detection = 0.01, prevalence = 50/100 )
core.ps1

```
#Core heatmaps
```{r}
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-2), log10(.2), length = 10), 3)

ps1.rel.spec<-aggregate_taxa(ps1.rel,'Species')

core.plot<-plot_core(ps1.rel.spec,
                     plot.type = 'heatmap',
                     colours = rev(brewer.pal(5,'Spectral')),
                    prevalences=prevalences,detections = detections,
                    min.prevalence = 0.5)+
  #scale_y_discrete(expand=c(0,0))+
  theme(axis.text.x= element_text(size=8),
        axis.text.y= element_text(size=8,face="italic"),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

core.plot

```

#load metadata and join with picobirna contig data
```{r}
meta<-fread('Freeze3_metadata.csv',header = TRUE)
#join metadata and pico contig data
merged.data<-left_join(meta,p.df1,by='seqID_virome')
#str(merged.data)

# convert na's to zeros
merged.data[is.na(merged.data)]=0

#reassign 0 to absent in Present column
merged.data$Present[merged.data$Present=='0']<-'absent'
#str(merged.data)

final.data<-merged.data%>%
  select(database_ID,contigID,RPKM,gene,pico_RA,Present)
str(final.data)
```
#Create dataframe from phyloseq object with high prevalent taxa
```{r}
taxa_names(ps1.rel) 

core.p<-core(ps1.rel,detection = 0.01,prevalence=50/100)
core.p
taxa_names(core.p)

#create dataframe from phyloseq object
core.df<-as_tibble(psmelt(core.p))

#merge df from phyloseq object to final metadat with picobirnavirus
ps.df<-left_join(core.df,final.data,by='database_ID')
str(ps.df)
```
#Create dataframe from phyloseq object without high prevalent taxa
```{r}
taxa_names(ps1.rel) 

non.core.p<-core(ps1.rel,detection = 0,prevalence=0)
non.core.p
taxa_names(non.core.p)

#create dataframe from phyloseq object
non.core.df<-as_tibble(psmelt(non.core.p))

#merge df from phyloseq object to final metadat with picobirnavirus
non.ps.df<-left_join(non.core.df,final.data,by='database_ID')

```
#relative abundance plot of  contigs
```{r}
str(non.ps.df)
#convert all string into factors
non.ps.df[sapply(non.ps.df, is.character)] <- lapply(non.ps.df[sapply(non.ps.df, is.character)], as.factor)
str(non.ps.df)
```
#contig plot with rel abundance of taxa
```{r}
non.ps.df<-non.ps.df[!is.na(non.ps.df$contigID), ]
#non.ps.df<-non.ps.df[!is.na(non.ps.df$Species), ]


contig.ps.plot<-non.ps.df%>%
  #filter(Abundance>0)%>%
  ggplot(aes(x=contigID,y=Abundance,fill=Family))+
  geom_bar(stat = 'identity',position='fill',width = 1)+
  facet_wrap(~Phylum)+
  scale_fill_brewer(palette = 'Paired')+
  theme(axis.text.x=element_text(angle = 90))

contig.ps.plot

ps.df<-ps.df[!is.na(ps.df$contigID), ]
contig.core.plot<-ps.df%>%
  ggplot(aes(x=contigID,y=Abundance,fill=Species))+
  geom_bar(stat = 'identity',position='fill',width = 1)+
  scale_fill_brewer(palette = 'Paired')+
  theme(axis.text.x=element_text(angle = 90))

contig.core.plot

```
#linear model family
```{r}
non.ps.df$RPKM.log<-log(non.ps.df$RPKM)

ggscatter(non.ps.df,x='RPKM.log',y='Abundance', size = 1,
                        add='reg.line',conf.int = TRUE,
          facet.by = 'Species',
          add.params = list(color = "blue",
                            fill = "lightgray"))+
  stat_cor(method = "spearman",)


ggsave('non_corepicobirna_Species.png',width = 20,height = 15)


```
