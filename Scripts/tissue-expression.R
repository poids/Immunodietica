#Script imports datframe containing results of epitope-uniprot query, creates plots, and exports as png

#Written by Vasco Morais for Immunodietica
#Cell: (415)-845-2118

#Import Packages
library('dplyr')
library('ggplot2')
library('stringr')

#Sets working directory
args = commandArgs(trailingOnly=TRUE)
run_path=as.character(args[1])
setwd(run_path)

###FUNCTIONS###
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
  ggsave(paste0(plot_type, ".png"), plot = p, path = './results/',
         scale = 1, width = p_width, height = p_height, units = "in",
         dpi = "retina", limitsize = TRUE)
}


#START OF SCRIPT#
#Import Tissue-Specificity from Human Protein Atlas 
tissue=read.csv('data/normal_tissue.tsv', header = TRUE, sep = "\t")

#Import comparative labels for tissue at system and organ level
system=read.csv('data/tissuesystem.csv', header = TRUE, sep = ",")
colnames(system)=c("Tissue", "Organ", "System")

#Add Tissue/Organ System row to tissue dataframe
tis_expr <- merge(tissue, system, by='Tissue') %>%
  select(Gene, Gene.name, Organ, System, Tissue, Cell.type, Level, Reliability)

#Remove unused entries (Not enough data form Human Protein Atlas to be accurate)
tis_expr=tis_expr[!(tis_expr$Organ=='Pituitary' | tis_expr$Organ=='Eye' | tis_expr$Organ=='Foot' | tis_expr$Organ=='Thymus'),]

#Import dataframe contianing results of epitope-uniprot query
if (as.character(args[2])=='prot') {
  uniprot=read.csv('results/protname_gene.csv', header = TRUE, sep = ",", stringsAsFactors = F)
  uniprot$Gene_Symbol<-toupper(uniprot$Gene_Symbol)#Make all proteins uppercase
} else if (as.character(args[2])=='acc') {
  uniprot=read.csv('results/acc_gene.csv', header = TRUE, sep = ",", stringsAsFactors = F)
  uniprot$Gene_Symbol<-toupper(uniprot$Gene)#Make all proteins uppercase
} else {
  print("No dataframe for epitope-uniprot query results found!")
}

#Options to filter based on Organism and AutoImmune disease
# uniprot=uniprot[,1:4]
# colnames(uniprot)=c('Epitope','UniProt_Name', 'Gene_Symbol', 'Disease')
# uniprot=filter(uniprot, Disease=='ulcerative colitis,')

#Iterates through entries that have multiple protein aliases and replaces it with the one usesd by Human Protein Atlas
protein_aliases<-uniprot %>% filter(str_detect(Gene_Symbol, "//")) #Finds all rows where aliases are used
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


#Make list of proteins
genes=array(uniprot$Gene)
#Count frequency of each protein in list
prot_freq=as.data.frame(table(genes))
#Format column names so it can be merged with tissue_expression data
colnames(prot_freq) = c('Gene.name', 'Freq')

#Writes dataframe containing tissue-specificity and metadata for the epitope subset used in the uniprot query
tis_expr %>%
  filter(Gene.name %in% genes) %>% #SHould this be other way around??
  filter(Level!='Not detected') %>%
  filter(Reliability!='Uncertain') %>%
  select(-Tissue, -Cell.type) %>%
  merge(uniprot, ., by.x='Gene_Symbol', by.y='Gene.name') %>%
  merge(epitope.disease, ., by='Epitope') %>%
  write.csv(., "./results/epitope_tissue-specificity.csv", row.names = FALSE, quote = TRUE)

#Creates tibble_df for System-Level Expression
tis_system<-tis_expr %>%
  filter(Gene.name %in% genes) %>%
  filter(Level!='Not detected') %>%
  filter(Reliability!='Uncertain') %>%
  select(Gene.name, System) %>%
  distinct() %>%
  merge(., prot_freq, by='Gene.name') %>%
  count(System, wt=Freq)

#Creates tibble_df for Organ-Level Expression
tis_organ<-tis_expr %>%
  filter(Gene.name %in% genes) %>%
  filter(Level!='Not detected') %>%
  filter(Reliability!='Uncertain') %>%
  select(Gene.name, Organ) %>%
  distinct() %>%
  ubiquitous(., length(unique(tis_expr$Organ))) %>% #Ubiquitous function
  merge(., prot_freq, by='Gene.name') %>%
  count(Organ, wt=Freq)


#Colors for plot
mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")

#Plots System-Level Expression
p<-ggplot(data=tis_system, aes(x=System, y=n)) +
  geom_col(colour="black", width=0.8, aes(fill=n)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1)) +
  ggtitle('Epitope Tissue-Expression (System Level)') + ylab('Epitope Count') +
  scale_fill_gradient(high = "firebrick", low = "dodgerblue3", name = "# of Epitopes") +
  scale_y_continuous(expand = c(0, 0))

#Saves PNG of Tissue-System Plot
save_plot(p, "system-expression")

#Plots Organ-Level Expression
p<-ggplot(data=tis_organ, aes(x=Organ, y=n)) +
  geom_col(colour="black", width=0.8, aes(fill=n)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1)) +
  ggtitle('Epitope Tissue-Expression (Organ Level)') + ylab('Epitope Count') +
  scale_fill_gradient(high = "firebrick", low = "dodgerblue3", name = "# of Epitopes") +
  scale_y_continuous(expand = c(0, 0))

#Saves PNG of Organ Plot
save_plot(p, "organ-expression")