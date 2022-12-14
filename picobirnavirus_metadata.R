---
title: "Picobirnavirus association with IBD metadata"
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
meta<-as_tibble(fread(file = 'Freeze3_metadata.csv', stringsAsFactors = TRUE, header = TRUE))
```

```{r}
#only 777 samples were used in the picobirnavirus genome analysis.
#get sample ids for all the 777 samples used for picobirna genome analysis
#load viral big table
df.virus<-as_tibble(fread(file = 'RC2_IBD_freeze2_viral_bigtable.tsv', stringsAsFactors = TRUE, header = TRUE))
sampleIDs<-summarise(df.virus,sampleID)%>%
  distinct()

#extract only samples present in analysis data from viral_bigtable
new_meta<-filter(meta,seqID_virome %in% sampleIDs$sampleID)

# rename column in new metadata
meta_df<-rename(new_meta, sampleID=seqID_virome)

#count total number of samples in each disease in metadata
sample.diag.count<-meta_df%>%
  select(database_ID,diagnosis)%>%
  group_by(diagnosis)%>%
  summarise(total_count=n())
```

#count picobirnavirus samples in metadata
```{r}
# count hits for picobirnaviridae for each sample
df.picobirna<-df.virus%>%
  filter(family == 'Picobirnaviridae')%>%
  group_by(sampleID)%>%
  summarise(n = sum(CPM))%>%
  mutate(Pico_PA=ifelse(n>0,'present','absent'))
  
# merge in the metadata
# Using all=T will do an outer join and bring in the absent samples
meta_pico<-merge(df.picobirna,meta_df,by='sampleID',all = T)
str(meta_pico)

# convert na's to zeros
meta_pico[is.na(meta_pico)]=0
str(meta_pico)

#reassign 0 to absent in Pico_PA column
meta_pico$Pico_PA[meta_pico$Pico_PA=='0']<-'absent'

#count picobirna present and absent samples
pico.diag.count<-meta_pico%>%
  select(sampleID,diagnosis,Pico_PA)%>%
  distinct()%>%
  group_by(diagnosis,Pico_PA)%>%
  summarise(total_count=n())
```
# picobirnavirus samples count in IBD vs control
```{r}
#samples association with IBD only
df.picobirna.count<-meta_pico%>%
  select(sampleID,diagnosis,Pico_PA)%>%
  filter(diagnosis == "HHC" | diagnosis == "CD" | diagnosis == "UC") %>%
  distinct()%>%
  group_by(diagnosis,Pico_PA)%>%
  summarise(count=n())%>%
  group_by(diagnosis)%>%
  mutate(Percentage = round(count/sum(count)*100,2))
  
#create matrix  for present and absent for each diagnosis group
pico.HHC = c(
  df.picobirna.count%>%
    filter(diagnosis =='HHC', Pico_PA == 'present')%>%
    pull(count),
  df.picobirna.count%>%
    filter(diagnosis =='HHC',Pico_PA == 'absent')%>%
    pull(count))

pico.UC = c(
  df.picobirna.count%>%
    filter(diagnosis =='UC', Pico_PA == 'present')%>%
    pull(count),
  df.picobirna.count%>%
    filter(diagnosis =='UC',Pico_PA == 'absent')%>%
    pull(count))

pico.CD = c(
  df.picobirna.count%>%
    filter(diagnosis =='CD', Pico_PA == 'present')%>%
    pull(count),
  df.picobirna.count%>%
    filter(diagnosis =='CD',Pico_PA == 'absent')%>%
    pull(count))

# create the 2x2 matrix for each group
pico.fish.mtx1 = matrix(c(pico.HHC,pico.UC),nrow=2)
colnames(pico.fish.mtx1) = c('HHC','UC')
row.names(pico.fish.mtx1) = c('present','absent')

pico.fish.mtx2 = matrix(c(pico.HHC,pico.CD),nrow=2)
colnames(pico.fish.mtx2) = c('HHC','CD')
row.names(pico.fish.mtx2) = c('present','absent')

pico.fish.mtx3 = matrix(c(pico.UC,pico.CD),nrow=2)
colnames(pico.fish.mtx3) = c('UC','CD')
row.names(pico.fish.mtx3) = c('present','absent')

#run Fisher's exact test for HHC & UC
pico.fish.mtx1.fisher<-fisher.test(pico.fish.mtx1)
pico.fish.mtx1.fisher$p.value

#run Fisher's exact test for HHC & CD
pico.fish.mtx2.fisher<-fisher.test(pico.fish.mtx2)
pico.fish.mtx2.fisher$p.value

#run Fisher's exact test for UC & CD
pico.fish.mtx3.fisher<-fisher.test(pico.fish.mtx3)
pico.fish.mtx3.fisher$p.value

ggplot(df.picobirna.count,aes(x=reorder(diagnosis, count),y=count,fill = Pico_PA))+
  geom_bar(stat = 'identity',position='stack')+
  #geom_text(aes(label=Percentage),position=position_dodge(width=0.9),vjust=3)+
  labs(x=" ", y="Number of samples", title = 'Picobirnavirus samples in disease status', fill = 'Picobirnavirus')+
  scale_y_continuous(breaks=c(0,50,100,150,200,250,300))+
  scale_fill_brewer(palette="Dark2")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
  theme(axis.title.y = element_text(face="bold",size=12),
        axis.text.x = element_text(face="bold",size=12))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_text(size = 12))
```
#picorbirnavirus viral loads in disease vs control

```{r}
#picobirnavirus viral loads in each group
df.viral.count<-meta_pico%>%
  select(sampleID,diagnosis,Pico_PA,n,flare,inflammation)%>%
  filter(diagnosis == 'HHC' |diagnosis == "CD" | diagnosis == "UC",Pico_PA == 'present')

ggplot(df.viral.count,aes(x=reorder(diagnosis,-n),y=log10(n)))+
  geom_boxplot(aes(fill=diagnosis),outlier.shape = NA)+
  geom_jitter(width = 0.2)+
  labs(x=" ", y="log10(CPM)", title = 'Picobirnavirus abundance in IBD group', fill = 'IBD group')+
  scale_fill_brewer(palette="Dark2")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
  theme(axis.title.y = element_text(face="bold",size=12),
        axis.text.x = element_text(face="bold",size=12))

wilcox.test(df.viral.count%>%
              filter(diagnosis == 'CD')%>%
              pull(n),
            df.viral.count%>%
              filter(diagnosis == 'UC')%>%
              pull(n),
            alternative = 't',
            paired = F)

library(rstatix)
dunn_test<-df.viral.count%>%
  dunn_test(n~diagnosis,p.adjust.method='holm')%>%
  add_significance()
```

```{r}
#abundance association with flare

flare.plot<-df.viral.count%>%
  filter(diagnosis =='UC'|diagnosis=='CD',flare == 'yes'| flare =='no')%>%
  ggplot(aes(x=reorder(flare,-n),y=log10(n)))+
  geom_boxplot(aes(fill=flare),outlier.shape = NA)+
  geom_jitter(width = 0.2)+
  facet_grid(~diagnosis,scales = 'free',space = 'free_x')+
  theme(strip.background =element_rect(fill="black"))+
  theme(strip.text = element_text(colour = 'white',face = 'bold'))+
  labs(x=" ", y="log10(CPM)", title = ' ', color = 'IBD group')+
  scale_fill_brewer(palette="Dark2")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
  theme(axis.title.y = element_text(face="bold",size=12),
        axis.text.x = element_text(face="bold",size=10))

flare.plot

#test between flare type in CD
wilcox.test(df.viral.count%>%
              filter(diagnosis == 'CD',flare == 'yes')%>%
              pull(n),
            df.viral.count%>%
              filter(diagnosis == 'CD',flare == 'no')%>%
              pull(n),
            alternative = 't',
            paired = F)

#test between flare type in UC
wilcox.test(df.viral.count%>%
              filter(diagnosis == 'UC',flare == 'yes')%>%
              pull(n),
            df.viral.count%>%
              filter(diagnosis == 'UC',flare == 'no')%>%
              pull(n),
            alternative = 't',
            paired = F)


```

```{r}
#abundance association with inflammation

inflammation.plot<-df.viral.count%>%
  filter(diagnosis =='UC'|diagnosis=='CD',inflammation == 'active'| inflammation =='inactive')%>%
  ggplot(aes(x=reorder(inflammation,-n),y=log10(n)))+
  geom_boxplot(aes(fill=inflammation),outlier.shape = NA)+
  geom_jitter(width = 0.2)+
  facet_grid(~diagnosis,scales = 'free',space = 'free_x')+
  theme(strip.background =element_rect(fill="black"))+
  theme(strip.text = element_text(colour = 'white',face = 'bold'))+
  labs(x=" ", y="log10(CPM)", title = ' ', color = 'IBD group')+
  scale_fill_brewer(palette="Dark2")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
  theme(axis.title.y = element_text(face="bold",size=12),
        axis.text.x = element_text(face="bold",size=10))

inflammation.plot

#test between flare type in CD
wilcox.test(df.viral.count%>%
              filter(diagnosis == 'CD',inflammation == 'active')%>%
              pull(n),
            df.viral.count%>%
              filter(diagnosis == 'CD',inflammation == 'inactive')%>%
              pull(n),
            alternative = 't',
            paired = F)

#test between flare type in UC
wilcox.test(df.viral.count%>%
              filter(diagnosis == 'UC',inflammation == 'active')%>%
              pull(n),
            df.viral.count%>%
              filter(diagnosis == 'UC',inflammation == 'inactive')%>%
              pull(n),
            alternative = 't',
            paired = F)

```
