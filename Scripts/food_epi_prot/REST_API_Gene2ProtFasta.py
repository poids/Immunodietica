# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 18:40:44 2019

@author: Vasco Morais for Immunodietica
"""
import requests, sys

proteins=["ACAT1", "ACO2", "ACSL6"]
taxid=9913
organism="Bos taurus"

#Write FASTA File
file = open("{organism}_proteins.fa".format(organism=organism.replace(" ", "_")), "a")

for protein in proteins:

    requestURL = "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=1&gene={gene}&taxid={taxid}".format(gene=protein, taxid=taxid)
    
    r = requests.get(requestURL, headers={ "Accept" : "text/x-fasta"})
    
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    
    responseBody = r.text
    file.write(responseBody)
    print(responseBody)

file.close()