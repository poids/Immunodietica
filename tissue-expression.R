#Import Packages
library('dplyr')
library('ggplot2')
library('stringr')

#Set Working Directory
setwd('~/Bioinformatics/Immunodietica/')

#Import Tissue-Specificity
tissue=read.csv(unzip('Data/normal_tissue.tsv.zip'), header = TRUE, sep = "\t")

system=read.csv('Data/tissuesystem.csv', header = TRUE, sep = ",")
colnames(system)=c("Tissue", "Organ", "System")

#FUNCTION: If epitope is found in all systems it changes it so it just reads as ubiquitous instead
ubiquitous = function(epitope_df, org_no) {
  #org_no=length(unique(epitope_df$Organ)) # of Organs in dataframe
  #epitope_df=add_row(epitope_df, Gene.name = "", Organ = 'Ubiquitous')#Inserts blank placeholder row for ubiquitous so it shows up on graph even if there are no hits
  for (epi in unique(epitope_df$Gene.name)) {
    if (sum(epitope_df$Gene.name==epi)==org_no) { #If gene is found in every tissue system,
      epitope_df=epitope_df[!(epitope_df$Gene.name==epi),] #Removes all rows where gene appears
      epitope_df=add_row(epitope_df, Gene.name = epi, Organ = 'Ubiquitous')#And inserts row where gene is counted as ubiquitous instead
    }
  }
  return(epitope_df)
}

#Add Tissue/Organ System row to tissue dataframe
tis_expr <- merge(tissue, system, by='Tissue') %>%
  select(Gene, Gene.name, Organ, System, Tissue, Cell.type, Level, Reliability)

#Remove entries for pituitary gland and eye (Not enough data form Human Protein Atlas to be accurate)
tis_expr=tis_expr[!(tis_expr$Organ=='Pituitary' | tis_expr$Organ=='Eye' | tis_expr$Organ=='Foot' | tis_expr$Organ=='Thymus'),]

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
uniprot=read.csv('Data/protname_gene.csv', header = TRUE, sep = ",", stringsAsFactors = F)

#PIG_TEST
#colnames(uniprot)<-c('Epitope', 'ACC', 'Gene_Symbol')

uniprot$Gene_Symbol<-toupper(uniprot$Gene_Symbol)#Make all proteins uppercase

#Options to filter based on Organism and AutoImmune disease
# disease_epitopes <- read.csv("~/Bioinformatics/Immunodietica/Data/disease_epitopes.csv", sep="")
# uniprot=merge(uniprot, disease_epitopes, by.x = 'Epitope', by.y='description')
# uniprot=uniprot[,1:4]
# colnames(uniprot)=c('Epitope','UniProt_Name', 'Gene_Symbol', 'Disease')
# uniprot=filter(uniprot, Disease=='ulcerative colitis,')

#Iterates through entries that have multiple protein aliases and replaces it with the one usesd by Human Protein Atlas
protein_aliases<-uniprot %>% filter(str_detect(Gene_Symbol, "//")) #Finds all rows wher e
for (i in unique(protein_aliases$Gene_Symbol)) {
  for (j in strsplit(i, "//")) {
    for (k in j) {
      if (k %in% unique(tis_expr$Gene.name)) {
        uniprot$Gene_Symbol[uniprot$Gene_Symbol==i] <- k
        #print(paste0(i, " <- ", k))
        next
      }
    }
  }
}

genes=array(uniprot$Gene) #Make list of proteins

prot_freq=as.data.frame(table(genes)) #Count frequency of each protein in list
colnames(prot_freq) = c('Gene.name', 'Freq') #Format column names so it can be merged with tissue_expression data

#genes=c('Aqp4', 'FTL',	'HBB2',	'KCNJ10', 'MAG',	'MBP',	'MOG', 'Mag',	'Mog',	'PLP1',	'RTN4R',	'TALDO1',	'TKT',	'TUBB1',	'Taldo1',	'Tubb1',	'gpmA1',	'gpmA2')
#genes=sample(unique(tis_sys$Gene.name), 70)


###TESTING




#System
test_system<-tis_expr %>%
  filter(Gene.name %in% genes) %>% #SHould this be other way around??
  filter(Level!='Not detected') %>%
  filter(Reliability!='Uncertain') %>%
  select(Gene.name, System) %>%
  distinct() %>%
  merge(., prot_freq, by='Gene.name') %>%
  count(System, wt=Freq)

#Organ
test_organ<-tis_expr %>%
  filter(Gene.name %in% genes) %>%
  filter(Level!='Not detected') %>%
  filter(Reliability!='Uncertain') %>%
  select(Gene.name, Organ) %>%
  distinct() %>%
  ubiquitous(., length(unique(tis_expr$Organ))) %>% #Ubiquitous function
  merge(., prot_freq, by='Gene.name') %>%
  count(Organ, wt=Freq)
  
# test_organ<-data.frame(lapply(test_organ, rep, test_organ$Freq)) %>%
#   select(Gene.name, Organ) %>%
#   count(Organ)
    
  
#%>%
  #distinct()%>%
#  count(Organ)



mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")

#System
ggplot(data=test_system, aes(x=System, y=n)) +
  geom_col(colour="black", aes(fill=n)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1)) +
  ggtitle('Epitope Tissue-Expression (System Level)') + ylab('Epitope Count') +
  scale_fill_gradient(high = "firebrick", low = "#ffce00", name = "# of Epitopes")
  #scale_fill_gradientn(colours = mycol)

#Organ
ggplot(data=test_organ, aes(x=Organ, y=n)) +
  geom_col(colour="black",aes(fill=n)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1)) +
  ggtitle('Epitope Tissue-Expression (Organ Level)') + ylab('Epitope Count') +
  scale_fill_gradient(high = "firebrick", low = "#ffce00", name = "# of Epitopes")
  #scale_fill_gradientn(colours = mycol)
