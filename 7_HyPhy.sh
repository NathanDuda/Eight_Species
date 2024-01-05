

#!/bin/bash


# align proteins of orthogroups
for orthogroup_file in ./Grouped_Fastas/One_to_Two_Orthogroups/Protein_Fastas/*; do
 filename=$(basename "$orthogroup_file")
 muscle -in "$orthogroup_file" -out "./Grouped_Fastas/One_to_Two_Orthogroups/Protein_Alignments/${filename}"
done


# get codon alignment of orthogroups with pal2nal
for orthogroup_file in ./Grouped_Fastas/One_to_Two_Orthogroups/Protein_Alignments/*; do
 filename=$(basename "$orthogroup_file")
 pal2nal.pl "./Grouped_Fastas/One_to_Two_Orthogroups/Protein_Alignments/${filename}" "./Grouped_Fastas/One_to_Two_Orthogroups/Nucleotide_Fastas/${filename}" -output fasta > "./Grouped_Fastas/One_to_Two_Orthogroups/Codon_Alignments/${filename}"
done


# format tree files of orthogroups from OrthoFinder
for file in ./OrthoFinder_Output/Results_Jan01/Gene_Trees/*; do
 sed -i 's/d[^_]*_prot_//g' "$file"
done


# run hyphy busted on orthogroups 
for orthogroup_file in ./Grouped_Fastas/One_to_Two_Orthogroups/Codon_Alignments/*; do
 filename=$(basename "$orthogroup_file")
 filename="${filename%.*}"
 hyphy busted --alignment "./Grouped_Fastas/One_to_Two_Orthogroups/Codon_Alignments/${filename}.fa" --tree "./OrthoFinder_Output/Results_Jan01/Gene_Trees/${filename}_tree.txt" --output "./Evolutionary_Rate/HyPhy_bt_Orthogroup_Output/${filename}.txt" &
done
stop

