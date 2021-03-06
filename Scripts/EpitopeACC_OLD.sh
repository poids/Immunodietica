#!/bin/bash

chmod 775 $1

EPITOPES=`awk '{print $1}' $1`

echo 'Building table with peptide sequence & UniProt Accession ID'
for epitope in ${EPITOPES};
do
        #Runs grep search for UniProt name
        PROTEIN=`grep -m1 -B3 ${epitope} swissprot | grep '>' | cut -d '>' -f 2 | cut -d '.' -f 1 | head -1`

        VALID_GREP=`grep -c ${epitope} swissprot`
		
        #Checks to see if grep function works, if not it will try and run blast
        if [ ${VALID_GREP} == 0 ]
        then
                #Creates temporary fasta file for missing epitope
               echo ">"${epitope} > temp_fasta.fa
               echo ${epitope} >> temp_fasta.fa

                #Runs blast search on unmatched epitope
                PROTEIN=`blastp -query ./temp_fasta.fa -db swissprot -outfmt 5 -max_target_seqs 50 | grep 'Hit_def>' | awk -F 'Hit_def>' '{ print $2 }' | cut -d '.' -f 1 | head -1`
        fi > temp_fasta.fa

        #Removes temporary fasta file
        rm ./temp_fasta.fa

        #If no protein is found, add epitope to error_log file
        if [ -z "$PROTEIN" ]
        then
                echo ${epitope}
        fi >> unmatched_epitopes.txt

        paste <(printf %s "${epitope}") <(printf %s "$PROTEIN")
done > uniprot_ACC.tsv

echo 'Building dataframe with corresponding Gene/Protein Names'

#Runs python script to build dataframe containing ACC-ID and corresponding Gene/Protein Name
python3 ACC_2_GENE.py

#Summary Statistics
EPITOPE_NO=`awk '{print $2}' uniprot_ACC.tsv | wc -l`
UNIPROT=`awk '{print $2}' uniprot_ACC.tsv | grep ^[A-Z] | wc -l`
UNMATCHED=`wc -l unmatched_epitopes.txt | awk {'print $1'}`
GENE_NAME=`wc -l acc_gene.csv | awk {'print $1'}`

echo "Number of epitopes in original list: ${EPITOPE_NO}"
echo "Matching UniProt Accession IDs found in query: ${UNIPROT}"
echo "Number of Epitopes Unmatched: ${UNMATCHED}"
echo "Number of matching Protein/Gene Names ${GENE_NAME}"
