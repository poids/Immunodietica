#!/bin/bash

#PATHS
SCRIPTS="$( cd "$(dirname "$0")" ; pwd -P )"
MAIN_DIR=${SCRIPTS}/..
DATA=${SCRIPTS}/../data
RESULTS=${SCRIPTS}/../results

#Changes permissions so list of epitopes can be accessed
chmod 775 $1

#Import list of epitopes
EPITOPES=`awk '{print $1}' $1`

echo 'Building table with peptide sequence & UniProt Accession ID'
for epitope in ${EPITOPES};
do
        #Runs grep search for UniProt name
        PROTEIN=`grep -m1 -B3 ${epitope} ${DATA}/swissprot | grep '>' | cut -d '>' -f 2 | cut -d '.' -f 1 | head -1`

        #Creates temporary fasta file for missing epitope
        echo ">SEQUENCE" > temp_fasta.fa
        echo ${epitope} >> temp_fasta.fa


        #Checks to see if grep function works, if not it will try and run blast
        if [ -z "$PROTEIN" ]
        then
                #Runs blast search on unmatched epitope
                PROTEIN=`blastp -query ./temp_fasta.fa -db ${DATA}/swissprot -outfmt 5 -max_target_seqs 50 | grep 'Hit_def>' | awk -F 'Hit_def>' '{ print $2 }' | cut -d '.' -f 1 | head -1`
        fi

        #Removes temporary fasta file
        rm ./temp_fasta.fa

        #If no protein is found, add epitope to error_log file
        if [ -z "$PROTEIN" ]
        then
                echo ${epitope}
        fi >> ${RESULTS}/unmatched_epitopes.txt

        paste <(printf %s "${epitope}") <(printf %s "$PROTEIN")
done > ${RESULTS}/uniprot_ACC.tsv

echo 'Building dataframe with corresponding Gene/Protein Names'
echo 'finished!'

#Runs python script to build dataframe containing ACC-ID and corresponding Gene/Protein Name
python3 ${SCRIPTS}/ACC_2_GENE.py ${MAIN_DIR}

#Summary Statistics
EPITOPE_NO=`awk '{print $2}' uniprot_ACC.tsv | wc -l`
UNIPROT=`awk '{print $2}' uniprot_ACC.tsv | grep ^[A-Z] | wc -l`
UNMATCHED=`wc -l ${RESULTS}/unmatched_epitopes.txt | awk {'print $1'}`
GENE_NAME=`wc -l ${REUSLTS}/acc_gene.csv | awk {'print $1'}`

echo "Number of epitopes in original list: ${EPITOPE_NO}" > ${RESULTS}/ACC_stats.txt
echo "Matching UniProt Accession IDs found in query: ${UNIPROT}" >> ${RESULTS}/ACC_stats.txt
echo "Number of Epitopes Unmatched: ${UNMATCHED}" >> ${RESULTS}/ACC_stats.txt
echo "Number of matching Protein/Gene Names ${GENE_NAME}" >> ${RESULTS}/ACC_stats.txt

#Print summary statistics to console
cat ${RESULTS}/ACC_stats.txt

#Run RScript to produce final dataframes and barplots (Results Written to Results Subdirectory)
Rscript ${SCRIPTS}/tissue-expression.R ${MAIN_DIR} acc