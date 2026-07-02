#!/bin/bash

#Read a list showing taxonomic identification of MAGs from archaea from the second line
sed 1d mags_archaea_id.txt | while read -r line; do (

	#Get MAG name
	mag_name=$(echo $line | awk '{print $1}')
	#Get the taxonomic classification of each MAG
	tax_id=$(echo $line | awk -F " " '{$1=""; print}')
	#Remove the last character of taxonomic classification of each MAG if it is a white space
	tax_id=$(echo $tax_id | sed 's/ *$//')
	#Format sp taxonomic identifications
	tax_id="${tax_id// sp/ sp }"
	#Replace whitespaces by _
	tax_id_without_spaces="${tax_id// /_}"
	#Remove sp codes for searching in NCBI database
	tax_id_without_spaces="${tax_id_without_spaces//_sp_*/_sp}"
	#Remove additional one-letter codes for searching in NCBI database
	tax_id_without_spaces="${tax_id_without_spaces//_?_/}"
	#Remove additional two-letter codes for searching in NCBI database
	tax_id_without_spaces="${tax_id_without_spaces//_??_/}"
	#Remove additional one-letter codes at the end for searching in NCBI database
	tax_id_without_spaces=$(echo $tax_id_without_spaces | sed -r 's|_.?$||g')
	
	#Get taxa id codes from ncbi
	echo "Searching for " $tax_id_without_spaces " in the NCBI web portal..."
    tax_id_link=$(curl -L https://www.ncbi.nlm.nih.gov/taxonomy/?term="${tax_id_without_spaces}"&report=taxid&format=text)
    echo "NCBI code was retrieved..."
	tax_id_link=${tax_id_link//*wwwtax.cgi?id=/}
	#Remove web content before taxa id
    tax_id_link=${tax_id_link//ref=*/} 
	#Remove web content after taxa id
    tax_id_link=${tax_id_link//'"'*/}
	
	
	
	
	#If no NCBI code is found, re-format the input (conditional 1)
	if echo "$tax_id_link" | grep -q "xml"; then
	#Print an error message
	echo "No NCBI code was found, re-formatting the input..."
	#Get the name of the problematic taxon and remove the "sp"
	problematic_taxon_1="${tax_id_without_spaces//_sp/}"	
	#Add spaces to the problematic taxon to search in the complete taxonomic identification results
	problematic_taxon_1="${problematic_taxon_1//_/ }"
	#Get complete taxonomic identification of the problematic taxon
	complete_tax_id_1=$(grep "$mag_name" gtdbtk_output/gtdbtk.bac120.summary.tsv | cut -f2)
	#Get taxonomic ranks other than problematic taxon
	complete_tax_id_1="${complete_tax_id_1//"${problematic_taxon_1}"*/}"
	#Format lines where taxonomic identification at species level was not reach (-r should be indicated to use .?)
	complete_tax_id_1=$(echo "$complete_tax_id_1" | sed -r 's|;.?__$||g')
	#Show taxonomic identification of each MAG at non-problematic corresponding level
	complete_tax_id_1=$(echo "$complete_tax_id_1" | sed 's|d__.*__||g')
	#Replace whitespaces by _ (not recommended to search taxa like ER4 sp, CAG-217 sp and VUNA01 sp)
	complete_tax_id_1="${complete_tax_id_1// /_}"
	#Get taxa id codes from ncbi
	echo "Searching for " $complete_tax_id_1 " instead of " $problematic_taxon_1 " in the NCBI web portal... (first conditional)"
    tax_id_link=$(curl -L https://www.ncbi.nlm.nih.gov/taxonomy/?term="${complete_tax_id_1}"&report=taxid&format=text)
    echo "NCBI code was retrieved..."
	tax_id_link=${tax_id_link//*wwwtax.cgi?id=/}
	#Remove web content before taxa id
    tax_id_link=${tax_id_link//ref=*/} 
	#Remove web content after taxa id
    tax_id_link=${tax_id_link//'"'*/}
	
	
	
	
	#If no NCBI code is found, re-format the input (conditional 2)
	if echo "$tax_id_link" | grep -q "xml"; then
	#Print an error message
	echo "No NCBI code was found, re-formatting the input..."
	#Get the name of the problematic taxon and remove the "sp"
	problematic_taxon_2="${complete_tax_id_1//_sp/}"	
	#Add spaces to the problematic taxon to search in the complete taxonomic identification results
	problematic_taxon_2="${problematic_taxon_2//_/ }"
	#Get complete taxonomic identification of the problematic taxon
	complete_tax_id_2=$(grep "$mag_name" gtdbtk_output/gtdbtk.bac120.summary.tsv | cut -f2)
	#Get taxonomic ranks other than problematic taxon
	complete_tax_id_2="${complete_tax_id_2//"${problematic_taxon_2}"*/}"
	#Format lines where taxonomic identification at species level was not reach (-r should be indicated to use .?)
	complete_tax_id_2=$(echo "$complete_tax_id_2" | sed -r 's|;.?__$||g')
	#Show taxonomic identification of each MAG at non-problematic corresponding level
	complete_tax_id_2=$(echo "$complete_tax_id_2" | sed 's|d__.*__||g')
	#Replace whitespaces by _ (not recommended to search taxa like ER4 sp, CAG-217 sp and VUNA01 sp)
	complete_tax_id_2="${complete_tax_id_2// /_}"
	#Get taxa id codes from ncbi
	echo "Searching for " $complete_tax_id_2 " instead of " $problematic_taxon_2 " in the NCBI web portal... (second conditional)"
    tax_id_link=$(curl -L https://www.ncbi.nlm.nih.gov/taxonomy/?term="${complete_tax_id_2}"&report=taxid&format=text)
    echo "NCBI code was retrieved..."
	tax_id_link=${tax_id_link//*wwwtax.cgi?id=/}
	#Remove web content before taxa id
    tax_id_link=${tax_id_link//ref=*/} 
	#Remove web content after taxa id
    tax_id_link=${tax_id_link//'"'*/}
	
	
	
	
	#If no NCBI code is found, re-format the input (conditional 3)
	if echo "$tax_id_link" | grep -q "xml"; then
	#Print an error message
	echo "No NCBI code was found, re-formatting the input..."
	#Get the name of the problematic taxon and remove the "sp"
	problematic_taxon_3="${complete_tax_id_2//_sp/}"	
	#Add spaces to the problematic taxon to search in the complete taxonomic identification results
	problematic_taxon_3="${problematic_taxon_3//_/ }"
	#Get complete taxonomic identification of the problematic taxon
	complete_tax_id_3=$(grep "$mag_name" gtdbtk_output/gtdbtk.bac120.summary.tsv | cut -f2)
	#Get taxonomic ranks other than problematic taxon
	complete_tax_id_3="${complete_tax_id_3//"${problematic_taxon_3}"*/}"
	#Format lines where taxonomic identification at species level was not reach (-r should be indicated to use .?)
	complete_tax_id_3=$(echo "$complete_tax_id_3" | sed -r 's|;.?__$||g')
	#Show taxonomic identification of each MAG at non-problematic corresponding level
	complete_tax_id_3=$(echo "$complete_tax_id_3" | sed 's|d__.*__||g')
	#Replace whitespaces by _ (not recommended to search taxa like ER4 sp, CAG-217 sp and VUNA01 sp)
	complete_tax_id_3="${complete_tax_id_3// /_}"
	#Get taxa id codes from ncbi
	echo "Searching for " $complete_tax_id_3 " instead of " $problematic_taxon_3 " in the NCBI web portal... (third conditional)"
    tax_id_link=$(curl -L https://www.ncbi.nlm.nih.gov/taxonomy/?term="${complete_tax_id_3}"&report=taxid&format=text)
    echo "NCBI code was retrieved..."
	tax_id_link=${tax_id_link//*wwwtax.cgi?id=/}
	#Remove web content before taxa id
    tax_id_link=${tax_id_link//ref=*/} 
	#Remove web content after taxa id
    tax_id_link=${tax_id_link//'"'*/}
	
	
	
	
	#If no NCBI code is found, re-format the input (conditional 4)
	if echo "$tax_id_link" | grep -q "xml"; then
	#Print an error message
	echo "No NCBI code was found, re-formatting the input..."
	#Get the name of the problematic taxon and remove the "sp"
	problematic_taxon_4="${complete_tax_id_3//_sp/}"	
	#Add spaces to the problematic taxon to search in the complete taxonomic identification results
	problematic_taxon_4="${problematic_taxon_4//_/ }"
	#Get complete taxonomic identification of the problematic taxon
	complete_tax_id_4=$(grep "$mag_name" gtdbtk_output/gtdbtk.bac120.summary.tsv | cut -f2)
	#Get taxonomic ranks other than problematic taxon
	complete_tax_id_4="${complete_tax_id_4//"${problematic_taxon_4}"*/}"
	#Format lines where taxonomic identification at species level was not reach (-r should be indicated to use .?)
	complete_tax_id_4=$(echo "$complete_tax_id_4" | sed -r 's|;.?__$||g')
	#Show taxonomic identification of each MAG at non-problematic corresponding level
	complete_tax_id_4=$(echo "$complete_tax_id_4" | sed 's|d__.*__||g')
	#Replace whitespaces by _ (not recommended to search taxa like ER4 sp, CAG-217 sp and VUNA01 sp)
	complete_tax_id_4="${complete_tax_id_4// /_}"
	#Get taxa id codes from ncbi
	echo "Searching for " $complete_tax_id_4 " instead of " $problematic_taxon_4 " in the NCBI web portal... (fourth conditional)"
    tax_id_link=$(curl -L https://www.ncbi.nlm.nih.gov/taxonomy/?term="${complete_tax_id_4}"&report=taxid&format=text)
    echo "NCBI code was retrieved..."
	tax_id_link=${tax_id_link//*wwwtax.cgi?id=/}
	#Remove web content before taxa id
    tax_id_link=${tax_id_link//ref=*/} 
	#Remove web content after taxa id
    tax_id_link=${tax_id_link//'"'*/}
	
	
	
	
	#End of fourth conditional
	fi
	#End of third conditional
	fi
	#End of second conditional
	fi
	#End of first conditional
	fi
	
	#Print taxa name and NCBI code per iteration
	echo "MAG:" $mag_name "Taxa:" $tax_id_without_spaces "NCBI:" $tax_id_link > temp_"${mag_name}"_archaea_ncbi_codes

); done