#!/bin/bash

chmod 775 $1

EPITOPES=`awk '{print $1}' $1`

echo 'Building table with peptide sequence & UniProt names'
for epitope in ${EPITOPES};
do
	#Runs grep search for UniProt name
	PROTEIN=`grep -m1 -B20 ${epitope} swissprot | grep '>' | tail -1 | cut -d '=' -f 2 | cut -d ";" -f 1`

	VALID_GREP=`grep -c ${epitope} swissprot`

	#Checks to see if grep function works, if not it will try and run blast
	if [ ${VALID_GREP} == 0 ]
	then
		#Creates temporary fasta file for missing epitope
		echo ">"${epitope}
		echo ${epitope}

		#Runs blast search on unmatched epitope
		PROTEIN=`blastp -query ./temp_fasta.fa -db swissprot -outfmt 5 -max_target_seqs 5 | grep 'Full=' | awk -F 'Full=' '{ print $2 }' | awk -F '</Hit_def>' '{ print $1 }' | cut -d ';' -f 1 | head -1`
	fi > temp_fasta.fa

	#Removes temporary fasta file
	rm ./temp_fasta.fa

	#If no protein is found, add epitope to error_log file
	if [ -z "$PROTEIN" ]
	then
		echo ${epitope}
	fi >> unmatched_epitopes.txt

	paste <(printf %s "${epitope}") <(printf %s "$PROTEIN")
done > uniprot_names.tsv

echo 'Building dataframe with corresponding Gene/Protein Names'

#Runs python script to build dataframe containing ACC-ID and corresponding Gene/Protein Name
python3 PROTNAME_2_GENE.py

#Summary Statistics
EPITOPE_NO=`awk '{print $2}' uniprot_names.tsv | wc -l`
UNIPROT=`awk '{print $2}' uniprot_names.tsv | grep ^[A-Z0-9] | wc -l`
UNMATCHED=`wc -l unmatched_epitopes.txt | awk {'print $1'}`
GENE_NAME=`wc -l acc_gene.csv | awk {'print $1'}`

echo "Number of epitopes in original list: ${EPITOPE_NO}"
echo "Matching UniProt Protein Names found in query: ${UNIPROT}"
echo "Number of Epitopes Unmatched: ${UNMATCHED}"
echo "Number of matching Protein/Gene Names ${GENE_NAME}"