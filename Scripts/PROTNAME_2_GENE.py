'''
Code written by Vasco Morais (August 2019)
(415)-845-2118

DESCRIPTION: Takes dataframe with list of epitopes
and UniProt Protein Names (Obtained from SwissProt);
Runs query on db2db to obtain matching gene name
for each entry in the dataframe.
'''

###IMPORT AND CLEAN-UP UNIPROT PROTEIN NAMES FOR QUERY
import pandas as pd
import urllib.request, json
import sys
import os

#Import pipeline directory path and set working directoy
run_path=str(sys.argv[1])
os.chdir(run_path)

#Reads in dataframe with Epitope and UniProt Protein Name
df = pd.read_csv(r'results/uniprot_names.tsv', header=None, delimiter="\t")

#Appends column names to dataframe
df.columns=['Epitope', 'UniProt_Name']

#Makes Protein Names into list of unique names (removes duplicate search terms to speed up query)
uniprot=list(df.UniProt_Name.unique())

'''
Removes all unmatched (NaN) values from list
Converts non-alphanumeric characters in UniProt Name
to their hexcode equivalent so it is URL compatible
'''
cleanedList = [urllib.parse.quote(x) for x in uniprot if str(x) != 'nan']

'''
Can only query a certain max amount of entries at a time
Divides list of protein names into sublists for iteration
'''
uniprotList=[]
while (len(cleanedList)>50):
    subList=cleanedList[0:50]
    cleanedList=cleanedList[50:]
    uniprotList.append(list(subList))
else:
    uniprotList.append(list(cleanedList))
    
#Creates empty pandas dataframe for query results to be loaded into
sub_df=pd.DataFrame()

###db2db QUERY (Iterates through each SubList)

for uniprotSubList in uniprotList:
    #Formats UniProt Names list so it is compliant with db2db API Query
    uniprotNames=','.join(uniprotSubList)

    #Parameters for db2db Query (Change if necessary)
    method='db2db'
    format_type='row'
    input_type='uniprotproteinname'
    inputValues=uniprotNames
    outputs='genesymbol'
    taxonId='9606'

    json_url = "https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method{method}&format={format_type}&input={input_type}&inputValues={inputValues}&outputs={outputs}&taxonId={taxonId}".format(method=method, format_type=format_type, input_type=input_type, inputValues=inputValues, outputs=outputs, taxonId=taxonId)


    #Results imported as JSON
    with urllib.request.urlopen(json_url) as url:
        data = json.loads(url.read().decode())
        

    ###POST-PROCESSING AND CLEAN-UP OF QUERY RESULTS

    #Converts JSON to Pandas Dataframe
    json_df = pd.io.json.json_normalize(data)

    #Sets Column Headers
    json_df.columns=['Gene_Symbol', 'UniProt_Name']

    #Merges sub-query results with previous query results
    sub_df=pd.concat([sub_df, json_df], axis=0, ignore_index=True)

    #Drops duplicate entries in dataframe
    sub_df=sub_df.drop_duplicates()

#Merges query input from above with all query results
output=pd.merge(df, sub_df, how='inner', on='UniProt_Name')

#Drops duplicate entries in dataframe
#output=output.drop_duplicates()

#Writes output to csv file
output.to_csv(r'results/protname_gene.csv', sep=',', encoding='utf-8', index=False, header=True)