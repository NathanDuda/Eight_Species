#!/bin/bash


### HyPhy: 

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


# run hyphy busted on orthogroups (with at least 4 sequences bc then tree was made for it)

max_parallel=5
current_jobs=0

for orthogroup_file in ./Grouped_Fastas/One_to_Two_Orthogroups/Codon_Alignments/*; do
 if [[ $(grep -c '>' "$orthogroup_file") -ge 4 ]]; then
  filename=$(basename "$orthogroup_file")
  filename="${filename%.*}"
  hyphy busted --alignment "./Grouped_Fastas/One_to_Two_Orthogroups/Codon_Alignments/${filename}.fa" --tree "./OrthoFinder_Output/Results_Jan01/Gene_Trees/${filename}_tree.txt" --output "./Evolutionary_Rate/HyPhy_busted_Orthogroup_Output/${filename}.txt" &
 
  ((current_jobs++))
  if ((current_jobs >= max_parallel)); then
    wait
    current_jobs=0
  fi
 fi 
done
wait




### MEGA:

# align proteins of one to one orthogroups 
for orthogroup_file in ./Grouped_Fastas/One_to_One_Orthogroups/Protein_Fastas/*; do
 filename=$(basename "$orthogroup_file")
 muscle -in "$orthogroup_file" -out "./Grouped_Fastas/One_to_One_Orthogroups/Protein_Alignments/${filename}"
done


# get codon alignments of one to one orthogroups 
for orthogroup_file in ./Grouped_Fastas/One_to_One_Orthogroups/Protein_Alignments/*; do
 filename=$(basename "$orthogroup_file")
 pal2nal.pl "./Grouped_Fastas/One_to_One_Orthogroups/Protein_Alignments/${filename}" "./Grouped_Fastas/One_to_One_Orthogroups/Nucleotide_Fastas/${filename}" -output fasta > "./Grouped_Fastas/One_to_One_Orthogroups/Codon_Alignments/${filename}"
done


# get list of all files to run MEGA on 
rename 's/\.fa$/\.meg/' ./Grouped_Fastas/One_to_One_Orthogroups/Codon_Alignments/*.fa

for file in ./Grouped_Fastas/One_to_One_Orthogroups/Codon_Alignments/*.meg; do
  echo -e "#MEGA\nTITLE" | cat - "$file" > temp && mv temp "$file"
  sed -i 's/>/#/g' "$file"
  unix2dos "$file"
done

ls -d /mnt/c/Users/17735/Downloads/Eight_Species/Grouped_Fastas/One_to_One_Orthogroups/Codon_Alignments/*.meg > ./MEGA_Phylogeny/file_list.txt
sort ./MEGA_Phylogeny/file_list.txt -n > ./MEGA_Phylogeny/sorted_file_list.txt

sed 's|/mnt/c/Users/17735/Downloads/Eight_Species/Grouped_Fastas/One_to_One_Orthogroups/Codon_Alignments/OG0002391.meg||' ./MEGA_Phylogeny/sorted_file_list.txt > ./MEGA_Phylogeny/file_list.txt
sed 's|/mnt/c/Users/17735/Downloads/Eight_Species/Grouped_Fastas/One_to_One_Orthogroups/Codon_Alignments/OG0003946.meg||' ./MEGA_Phylogeny/file_list.txt > ./MEGA_Phylogeny/sorted_file_list.txt
sed 's|/mnt/c/Users/17735/Downloads/Eight_Species/Grouped_Fastas/One_to_One_Orthogroups/Codon_Alignments/OG0004744.meg||' ./MEGA_Phylogeny/sorted_file_list.txt > ./MEGA_Phylogeny/file_list.txt
sed 's|/mnt/c/Users/17735/Downloads/Eight_Species/Grouped_Fastas/One_to_One_Orthogroups/Codon_Alignments/OG0005168.meg||' ./MEGA_Phylogeny/file_list.txt > ./MEGA_Phylogeny/sorted_file_list.txt

tac ./MEGA_Phylogeny/sorted_file_list.txt > ./MEGA_Phylogeny/file_list.txt

rm ./MEGA_Phylogeny/sorted_file_list.txt

# run MEGA
megacc -a ./MEGA_Phylogeny/infer_NJ_amino_acid.mao -d ./MEGA_Phylogeny/file_list.txt -o ./MEGA_Phylogeny/MEGA_Output/MEGA_output


# concatenate all output newicks into one file 
cat ./MEGA_Phylogeny/MEGA_Output/*_consensus.nwk > ./MEGA_Phylogeny/MEGA_consensus_Newicks.txt
sed -E 's/(YOgn[A-Z][A-Z])[0-9]+:/\1:/g' ./MEGA_Phylogeny/MEGA_consensus_Newicks.txt > ./MEGA_Phylogeny/MEGA_consensus_Newicks_temp.txt
sed -E 's/(YOgn[A-Z][A-Z])[0-9]+,/\1:/g' ./MEGA_Phylogeny/MEGA_consensus_Newicks_temp.txt > ./MEGA_Phylogeny/MEGA_consensus_Newicks.txt
sed -E 's/(YOgn[A-Z][A-Z])[0-9]+/\1:/g' ./MEGA_Phylogeny/MEGA_consensus_Newicks.txt > ./MEGA_Phylogeny/MEGA_consensus_Newicks_temp.txt
sed 's/:/,/g' ./MEGA_Phylogeny/MEGA_consensus_Newicks_temp.txt > ./MEGA_Phylogeny/MEGA_consensus_Newicks.txt
sed 's/,)/)/g' ./MEGA_Phylogeny/MEGA_consensus_Newicks.txt > ./MEGA_Phylogeny/MEGA_consensus_Newicks_temp.txt
mv ./MEGA_Phylogeny/MEGA_consensus_Newicks_temp.txt ./MEGA_Phylogeny/MEGA_consensus_Newicks.txt


cat ./MEGA_Phylogeny/MEGA_Output/*\).nwk > ./MEGA_Phylogeny/MEGA_Newicks.txt
sed -E 's/(YOgn[A-Z][A-Z])[0-9]+:/\1:/g' ./MEGA_Phylogeny/MEGA_Newicks.txt > ./MEGA_Phylogeny/MEGA_Newicks_temp.txt
sed -E 's/(YOgn[A-Z][A-Z])[0-9]+,/\1:/g' ./MEGA_Phylogeny/MEGA_Newicks_temp.txt > ./MEGA_Phylogeny/MEGA_Newicks.txt


#rm ./MEGA_Phylogeny/MEGA_consensus_Newicks_temp.txt


