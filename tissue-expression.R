#Import Packages
library('dplyr')
library('ggplot2')
library('stringr')
library('svglite')
library('grDevices')

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

#FUNCTION: Saves plot as PNG and adjusts plot width based on number of columns
save_plot=function(p, plot_type){
  no_col=nrow(p$data)#Number of columns in graph
  p_height=500 #height of plot
  min_p_width=450 #min width of plot
  max_p_width=800 #max width of plot
  aspect_ratio=(max_p_width-min_p_width)/34 #Aspect Ratio
  if (no_col<=12) {
    p_width=min_p_width #assign plot width to minimum value
  } else {
    p_width=450+round(no_col*aspect_ratio) #Adujusts width of plot based on number of columns
  }
  
  #Converts plot width-height to inches by pixel/inch ratio (320="retina")
  ppi=100
  p_height=p_height/ppi
  p_width=p_width/ppi
  
  #Saves plot
  ggsave(paste0(plot_type, ".svg"), plot = p, path = './Data/results',
         scale = 1, width = p_width, height = p_height, units = "in",
         dpi = "retina", limitsize = TRUE)
}

#Add Tissue/Organ System row to tissue dataframe
tis_expr <- merge(tissue, system, by='Tissue') %>%
  select(Gene, Gene.name, Organ, System, Tissue, Cell.type, Level, Reliability)

#Remove entries for pituitary gland and eye (Not enough data form Human Protein Atlas to be accurate)
tis_expr=tis_expr[!(tis_expr$Organ=='Pituitary' | tis_expr$Organ=='Eye' | tis_expr$Organ=='Foot' | tis_expr$Organ=='Thymus'),]

#Import dataframe contianing results of epitope-uniprot query
uniprot=read.csv('Data/protname_gene.csv', header = TRUE, sep = ",", stringsAsFactors = F)
uniprot$Gene_Symbol<-toupper(uniprot$Gene_Symbol)#Make all proteins uppercase

  # uniprot=read.csv('Data/acc_gene.csv', header = TRUE, sep = ",", stringsAsFactors = F)
  # uniprot$Gene_Symbol<-toupper(uniprot$Gene)#Make all proteins uppercase


#Options to filter based on Organism and AutoImmune disease
uniprot=read.csv('Data/2019_Paper_Data/protname_gene.csv', header = TRUE, sep = ",", stringsAsFactors = F)
uniprot$Gene_Symbol<-toupper(uniprot$Gene_Symbol)#Make all proteins uppercase

disease_epitopes <- read.csv("~/Bioinformatics/Immunodietica/Data/2019_Paper_Data/epitope-disease.csv", stringsAsFactors=FALSE)
disease_epitopes <- disease_epitopes[c('description', 'disease')] #remove epitopeID
# colnames(disease_epitopes)<-c("Epitope_ID", "Epitope", "Disease")

#Disease to filter by:
# disease='Behcet syndrome'
# disease='insulin-dependent diabetes mellitus' #pancreas highest??
disease='systemic lupus erythematosus' #everywhere?
# disease='rheumatoid arthritis'
# disease='multiple sclerosis' #brain highest?
# disease="ulcerative colitis"

uniprot=merge(uniprot, disease_epitopes, by.x = 'Epitope', by.y='description')
uniprot=uniprot[,1:4]
colnames(uniprot)=c('Epitope','UniProt_Name', 'Gene_Symbol', 'Disease')
uniprot=filter(uniprot, Disease==disease)

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
epitope.disease <- read.csv("~/Bioinformatics/Immunodietica/Data/epitope-disease.csv", stringsAsFactors=FALSE)
colnames(epitope.disease)<-c("Epitope_ID", "Epitope", "Disease")

#Creates dataframe containing tissue-specificity and metadata for the epitope subset used in the uniprot query
# test<-tis_expr %>%
#   filter(Gene.name %in% genes) %>% #SHould this be other way around??
#   filter(Level!='Not detected') %>%
#   filter(Reliability!='Uncertain') %>%
#   select(-Tissue, -Cell.type) %>%
#   merge(uniprot, ., by.x='Gene_Symbol', by.y='Gene.name') %>%
#   merge(epitope.disease, ., by='Epitope') %>%
  #write.csv(., "DELETE.csv", row.names = FALSE, quote = TRUE)
  
#write.csv(test, "../Immunodietica (Personal)/full_human_autoimmune_data.csv", row.names = FALSE, col.names = TRUE, quote = TRUE, sep = ",")

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

#Capitalize letters in disease
disease<-gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", disease, perl=TRUE)

mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")

#System
p<-ggplot(data=test_system, aes(x=System, y=n)) +
  geom_col(colour="black", width=0.8, aes(fill=n)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1)) +
  ggtitle('Epitope Tissue-Expression (System Level)') + ylab('Epitope Count') +
  scale_fill_gradient(high = "firebrick", low = "dodgerblue3", name = "# of Epitopes") +
  scale_y_continuous(expand = c(0, 0))
  #scale_fill_gradientn(colours = mycol)

#Organ
p<-ggplot(data=test_organ, aes(x=Organ, y=n)) +
  geom_col(colour="black", width=0.8, aes(fill=n)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1)) +
  ggtitle(paste(disease, '-', 'Epitope Tissue-Expression (Organ Level)')) + ylab('Epitope Count') +
  # ggtitle(paste('All Autoimmune Diseases - Epitope Tissue-Expression (Organ Level)')) + ylab('Epitope Count') +
  scale_fill_gradient(high = "firebrick", low = "dodgerblue3", name = "# of Epitopes") +
  scale_y_continuous(expand = c(0, 0))
  #scale_fill_gradientn(colours = mycol)

p

save_plot(p, 'lupus_organ-expression')