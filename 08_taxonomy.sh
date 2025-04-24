#!/bin/bash

##################################
#     Bioinformatic pipeline     #
#     16s Bacteria #
#     ----------------------     #
#     taxonomy assignment     #
#                                #
#      Zachary Noel         #
##################################

module load vsearch 
#mkdir ./taxonomy

# Assign taxonomy using SINTAX algorithm
vsearch -sintax clustered/otus.fasta -db ~/noel_research/db_bacteria/SILVA_138.2_SSURef_NR99_tax_silva_newheaders.fasta -tabbedout taxonomy/16s_taxonomy.txt -strand both -sintax_cutoff 0.8


# Print the first and fourth columns of the sintax output, replace the commas with tabs, get rid of the taxonomic prefixes, save it in a text file
awk -F'\t' '{print $1, $4}' taxonomy/16s_taxonomy.txt | tr ',' '\t' |  sed 's/d://; s/p://; s/c://; s/o://; s/f://; s/g://; s/s://' > taxonomy/16s_taxonomy2.txt

# Create column headers
echo -e "OTU\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" > taxonomy/headers.txt

# append the taxonomy to the headers.txt
cat taxonomy/16s_taxonomy2.txt >> taxonomy/headers.txt

# rename the output 
mv taxonomy/headers.txt taxonomy/16s_taxonomy_SINTAX.txt