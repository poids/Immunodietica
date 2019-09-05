'''
Code written by Vasco Morais (August 2019)
(415)-845-2118

DESCRIPTION: Takes dataframe with list of epitopes
and UniProt Acession ID's (Obtained from SwissProt);
Runs query on UniProt to obtain matching gene name
for each entry in the dataframe.
'''

###IMPORT AND CLEAN-UP UNIPROT ACESSION ID'S FOR QUERY
import pandas as pd
import sys
import os

#Import pipeline directory path and set working directoy
run_path=str(sys.argv[1])
os.chdir(run_path)

#Reads in dataframe with Epitope and ACC ID
df = pd.read_csv(r'results/uniprot_ACC.tsv', header=None, delimiter="\t")

#Appends column names to dataframe
df.columns=['Epitope', 'ACC']

#Makes ACC IDs into List
ACC=list(df.ACC)

#Removes NaN values from list
cleanedList = [x for x in ACC if str(x) != 'nan']

#Reformats list into string so it can be used for UniProt Query
ACC=' '.join(cleanedList)

###UNIPROT QUERY (Following code chunk taken from UniProt)
import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uploadlists/'

##Query Parameters
params = {
'from': 'ACC+ID', # Queries using UniProt Accession ID#'s
'to': 'GENENAME', # Converts to Gene Name
'format': 'tab', # Output dataframe is tab delimited
'query': ACC # Input is variable for string created from dataframe above
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()
print(response.decode('utf-8')) #Prints query out to terminal (Comment out line if unwanted)


###POST-PROCESSING AND CLEAN-UP OF QUERY RESULTS
from pandas.compat import StringIO

#Creates tab-delimited dataframe out of query results
df1=response.decode('utf-8')

df1=pd.read_csv(StringIO(df1), delimiter='\t')

#Assigns column headers
df1.columns=['ACC', 'Gene']

#Merges query input from above with query results
output=pd.merge(df, df1, how='inner', on='ACC')

#Drops duplicate entries in dataframe
#output=output.drop_duplicates()

#Writes output to csv file
output.to_csv(r'results/acc_gene.csv', sep=',', encoding='utf-8', index=False, header=True)