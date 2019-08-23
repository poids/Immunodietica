#Import Packages
library('dplyr')
library('ggplot2')

#Set Working Directory
setwd('~/Bioinformatics/Immunodietica/')

#Import Tissue-Specificity
tissue=read.csv(unzip('Data/normal_tissue.tsv.zip'), header = TRUE, sep = "\t")

system=read.csv('Data/tissuesystem.csv', header = TRUE, sep = ",")
colnames(system)=c("Tissue", "Organ", "System")

#FUNCTION: If epitope is found in all systems it changes it so it just reads as ubiquitous instead
ubiquitous = function(epitope_df) {
  for (epi in unique(epitope_df$Gene.name)) {
    if (sum(epitope_df$Gene.name==epi)==11) { #If gene is found in every tissue system,
      epitope_df=epitope_df[!(epitope_df$Gene.name==epi),] #Removes all rows where gene appears
      epitope_df=add_row(epitope_df, Gene.name = epi, System = 'ubiquitous')#And inserts row where gene is counted as ubiquitous instead
    }
  }
  return(epitope_df)
}

#Add Tissue/Organ System row to tissue dataframe
tis_sys <- merge(tissue, system, by='Tissue') %>%
  select(Gene, Gene.name, Organ, System, Tissue, Cell.type, Level, Reliability)

#Count expression level
tissue %>%
  filter(Gene.name=='SOX4') %>%
  filter(Level!='Not detected') %>%
  filter(Reliability!='Uncertain') %>%
  group_by(Level) %>%
  count(Tissue)

#NOTES:
#Tried to threshold expression-level with only 'high' level matches, but many genes did not show up


#FOR TESTING
uniprot=read.csv('Data/protname_gene.csv', header = TRUE, sep = ",")
genes=array(uniprot$Gene)

#genes=c('Aqp4', 'FTL',	'HBB2',	'KCNJ10', 'MAG',	'MBP',	'MOG', 'Mag',	'Mog',	'PLP1',	'RTN4R',	'TALDO1',	'TKT',	'TUBB1',	'Taldo1',	'Tubb1',	'gpmA1',	'gpmA2')
#genes=sample(unique(tis_sys$Gene.name), 70)


#System
test_system<-tis_sys %>%
  filter(Gene.name %in% genes) %>%
  filter(Level!='Not detected') %>%
  filter(Reliability!='Uncertain') %>%
  select(Gene.name, System) %>%
  #distinct() %>%
  count(System)

#Organ
test_organ<-tis_sys %>%
  filter(Gene.name %in% genes) %>%
  filter(Level!='Not detected') %>%
  filter(Reliability!='Uncertain') %>%
  select(Gene.name, Organ) %>%
  #distinct()%>%
  count(Organ)

#Ubiquitous function (Not Needed?)
test <- ubiquitous(test) %>%
    count(System)


  


mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")

ggplot(data=test, aes(x=Tissue, y=Level)) +
  geom_tile(aes(fill=n)) +
  theme(axis.text.x = element_text(angle = 45,
                                    hjust = 1,
                                    vjust = 1)) +
  ggtitle('SOX4 Tissue Expression') + ylab('Epitope Count') +
  scale_fill_gradientn(colours = mycol)

#System
ggplot(data=test_system, aes(x=System, y=n)) +
  geom_col(aes(fill=n)) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1)) +
  ggtitle('Epitope Tissue-Expression') + ylab('Epitope Count') +
  scale_fill_gradientn(colours = mycol)

#Organ
ggplot(data=test_organ, aes(x=Organ, y=n)) +
  geom_col(aes(fill=n)) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1)) +
  ggtitle('Epitope Tissue-Expression') + ylab('Epitope Count') +
  scale_fill_gradientn(colours = mycol)
