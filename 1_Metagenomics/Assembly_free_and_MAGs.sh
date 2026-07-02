
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

#TORMES pipeline

#Paired-end metagenomic reads are in the "raw_reads" folder
#Metadata file "metadata.txt" is needed

#Define variables indicating different database folders
$METAPHLAN_DBDIR
$KRAKEN_DB1DIR
$KRAKEN_DB2DIR
$HMMER_DBDIR
$NCBIFAM_DBDIR
$CAZY_DBDIR
$TORMES_ENV

#Run TORMES
nohup tormes --metadata metadata.txt --output tormes_results --assembler megahit --threads 48 &> nohup_TORMES.out &
#Check assembly quality
mkdir fastqc_shotgun
fastqc --noextract --nogroup -o fastqc_shotgun tormes_results/megahit_*/genomes/*/*.fastq.gz 

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

#Assembly-free pipeline

#Concatenate R1 and R2
mkdir cat_reads
nohup bash -c 'for clean_read_r1 in tormes_results/megahit_*/cleaned_*/*ok_1*; do (clean_read_r2="${clean_read_r1//ok_1/ok_2}" && sample_name="${clean_read_r1//"/"*/}" && sample_name="${sample_name//megahit_/}" && zcat $clean_read_r1 $clean_read_r2 > cat_reads/"${sample_name}".fastq.gz); done' &> nohup_concatenate.out &
#Clean file extensions
rename fastq.gz fastq cat_reads/*fastq*
#Functional and taxonomic analysis using Humann 3.0
mkdir humann3_out
nohup parallel -j 9 'humann --threads 4 --input {} --output humann3_out/{/.} --metaphlan-options "-x mpa_v30_CHOCOPhlAn_201901 --bowtie2db $METAPHLAN_DBDIR/metaphlan_db"' ::: cat_reads/*fastq &> nohup_humann.out &
#Merge results from individual samples
conda activate biobakery3
mkdir humann3_final_out
humann_join_tables -s --input humann3_out/ --file_name pathabundance --output humann3_final_out/humann3_pathabundance.tsv
humann_join_tables -s --input humann3_out/ --file_name pathcoverage --output humann3_final_out/humann3_pathcoverage.tsv
humann_join_tables -s --input humann3_out/ --file_name genefamilies --output humann3_final_out/humann3_genefamilies.tsv
#Re-normalize gene families and metabolic pathways abundances
humann_renorm_table --input humann3_final_out/humann3_pathabundance.tsv --units cpm --output humann3_final_out/humann3_pathabundance_cpm.tsv
humann_renorm_table --input humann3_final_out/humann3_genefamilies.tsv --units cpm --output humann3_final_out/humann3_genefamilies_cpm.tsv
#Split taxonomic contributions
humann_split_stratified_table --input humann3_final_out/humann3_pathabundance_cpm.tsv --output humann3_final_out
humann_split_stratified_table --input humann3_final_out/humann3_genefamilies_cpm.tsv --output humann3_final_out
humann_split_stratified_table --input humann3_final_out/humann3_pathcoverage.tsv --output humann3_final_out
#Extract taxonomic composition data
mkdir metaphlan3_out
#Copy metaphlan outputs for individual samples
cp humann3_out/*/*/*metaphlan_bugs_list.tsv metaphlan3_out/
#Merge metaphlan tables
merge_metaphlan_tables.py metaphlan3_out/*metaphlan_bugs_list.tsv > metaphlan3_merged.txt
#Clean merged metaphlan tables
sed -i 's/_metaphlan_bugs_list//g' metaphlan3_merged.txt

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

#Kraken pipeline

#Copy TORMES outputs
mkdir cleaned_reads
cp megahit_*/cleaned_*/*ok_* cleaned_reads/

#Taxonomic identification with kraken
#parallel -j 9 'sample={} && sample="${sample//_1.*/}" && sample="${sample//*"/"/}" && echo "Sample: " "${sample}" " First iteration: R1 " cleaned_reads/"${sample}"_1.fastq.gz " R2: " cleaned_reads/"${sample}"_2.fastq.gz' ::: cleaned_reads/*_1*.fastq.gz
#Database 1: /opt/home/carlos/bases_datos_jose
mkdir kraken_out
nohup parallel -j 9 'sample={} && sample="${sample//.ok_1.*/}" && sample="${sample//*"/"/}" && kraken2 --db $KRAKEN_DB1DIR \
  cleaned_reads/"${sample}".ok_1.fastq.gz cleaned_reads/"${sample}".ok_2.fastq.gz \
  --paired --use-names --memory-mapping \
  --confidence 0.5 \
  --minimum-hit-groups 5 \
  --minimum-base-quality 20 \
  --threads 4 \
  --report kraken_out/"${sample}"_kraken_report \
  --classified-out kraken_out/"${sample}"_kraken_class# \
  --unclassified-out kraken_out/"${sample}"_kraken_unclass# \
  --output kraken_out/"${sample}"_kraken_output' ::: cleaned_reads/*.ok_1.fastq.gz &> nohup_kraken.out &

#Database 2 (maxikraken)
mkdir maxikraken_out
nohup parallel -j 9 'sample={} && sample="${sample//.ok_1.*/}" && sample="${sample//*"/"/}" && kraken2 --db $KRAKEN_DB2DIR \
  cleaned_reads/"${sample}".ok_1.fastq.gz cleaned_reads/"${sample}".ok_2.fastq.gz \
  --paired --use-names --memory-mapping \
  --confidence 0.5 \
  --minimum-hit-groups 5 \
  --minimum-base-quality 20 \
  --threads 4 \
  --report maxikraken_out/"${sample}"_kraken_report \
  --classified-out maxikraken_out/"${sample}"_kraken_class# \
  --unclassified-out maxikraken_out/"${sample}"_kraken_unclass# \
  --output maxikraken_out/"${sample}"_kraken_output' ::: cleaned_reads/*.ok_1.fastq.gz &> nohup_maxikraken.out &

#Bracken estimations

#Get results - strain level
mkdir kraken_bracken_out_strain
mkdir maxikraken_bracken_out_strain
#Jose database
parallel -j 6 'bracken -d $TORMES_ENV/db/minikraken-DB/minikraken_8GB_20200312 -i {} -o kraken_bracken_out_strain/{/.}.bracken -l S1' ::: kraken_out/*_kraken_report 
#Maxikraken database
parallel -j 6 'bracken -d $TORMES_ENV/db/minikraken-DB/minikraken_8GB_20200312 -i {} -o maxikraken_bracken_out_strain/{/.}.bracken -l S1' ::: maxikraken_out/*_kraken_report 
#Generate abundance output - collate the Bracken2 results - species level 
mkdir kraken_bracken_out_merged_strain
mkdir maxikraken_bracken_out_merged_strain
combine_bracken_outputs.py --files kraken_bracken_out_strain/*.bracken -o kraken_bracken_out_merged_strain/kraken_strain.tsv
combine_bracken_outputs.py --files maxikraken_bracken_out_strain/*.bracken -o maxikraken_bracken_out_merged_strain/maxikraken_strain.tsv

#Get results - species level
mkdir kraken_bracken_out_species
mkdir maxikraken_bracken_out_species
#Jose database
parallel -j 6 'bracken -d $TORMES_ENV/db/minikraken-DB/minikraken_8GB_20200312 -i {} -o kraken_bracken_out_species/{/.}.bracken -l S' ::: kraken_out/*_kraken_report 
#Maxikraken database
parallel -j 6 'bracken -d $TORMES_ENV/db/minikraken-DB/minikraken_8GB_20200312 -i {} -o maxikraken_bracken_out_species/{/.}.bracken -l S' ::: maxikraken_out/*_kraken_report 
#Generate abundance output - collate the Bracken2 results - species level 
mkdir kraken_bracken_out_merged_species
mkdir maxikraken_bracken_out_merged_species
combine_bracken_outputs.py --files kraken_bracken_out_species/*.bracken -o kraken_bracken_out_merged_species/kraken_species.tsv
combine_bracken_outputs.py --files maxikraken_bracken_out_species/*.bracken -o maxikraken_bracken_out_merged_species/maxikraken_species.tsv

#Get results - genus level
mkdir kraken_bracken_out_genus
mkdir maxikraken_bracken_out_genus
#Jose database
parallel -j 6 'bracken -d $TORMES_ENV/db/minikraken-DB/minikraken_8GB_20200312 -i {} -o kraken_bracken_out_genus/{/.}.bracken -l G' ::: kraken_out/*_kraken_report 
#Maxikraken database
parallel -j 6 'bracken -d $TORMES_ENV/db/minikraken-DB/minikraken_8GB_20200312 -i {} -o maxikraken_bracken_out_genus/{/.}.bracken -l G' ::: maxikraken_out/*_kraken_report 
#Generate abundance output - collate the Bracken2 results - species level 
mkdir kraken_bracken_out_merged_genus
mkdir maxikraken_bracken_out_merged_genus
combine_bracken_outputs.py --files kraken_bracken_out_genus/*.bracken -o kraken_bracken_out_merged_genus/kraken_genus.tsv
combine_bracken_outputs.py --files maxikraken_bracken_out_genus/*.bracken -o maxikraken_bracken_out_merged_genus/maxikraken_genus.tsv

#Get results - family level
mkdir kraken_bracken_out_family
mkdir maxikraken_bracken_out_family
#Jose database
parallel -j 6 'bracken -d $TORMES_ENV/db/minikraken-DB/minikraken_8GB_20200312 -i {} -o kraken_bracken_out_family/{/.}.bracken -l F' ::: kraken_out/*_kraken_report 
#Maxikraken database
parallel -j 6 'bracken -d $TORMES_ENV/db/minikraken-DB/minikraken_8GB_20200312 -i {} -o maxikraken_bracken_out_family/{/.}.bracken -l F' ::: maxikraken_out/*_kraken_report 
#Generate abundance output - collate the Bracken2 results - family level 
mkdir kraken_bracken_out_merged_family
mkdir maxikraken_bracken_out_merged_family
combine_bracken_outputs.py --files kraken_bracken_out_family/*.bracken -o kraken_bracken_out_merged_family/kraken_family.tsv
combine_bracken_outputs.py --files maxikraken_bracken_out_family/*.bracken -o maxikraken_bracken_out_merged_family/maxikraken_family.tsv

#Move kraken outputs to a single folder
mkdir results_kraken
mv *kraken_* results_kraken/

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

#MAGs pipeline

#Copy TORMES outputs to obtain MAGs
mkdir mags
cp megahit_*/genomes/*fasta mags/

#Run bowtie2 on megahit results
mkdir mags_binning
parallel -j 36 'input_name={} && output_name="${input_name//.fa*/}" && output_name="${output_name//*"/"/}" && bowtie2-build "${input_name}" mags_binning/"${output_name}".final.contigs' ::: mags/*.fa*
#Copy genome assemblies to the same folder containing bowtie2 indexes
cp mags/*.fa* mags_binning/
#Rename assemblies to match bowtie2 indexes
for mags in mags_binning/*.fa* ; do (rename_sample="${mags//.fa/.final.contigs.fa}" && mv $mags $rename_sample); done
#Rename assemblies to match .fa format
for mags in mags_binning/*.fa* ; do (rename_sample="${mags//.fasta/.fa}" && mv $mags $rename_sample); done
#Map the reads against the assembly to get coverage information 
nohup parallel -j 36 'input_name={} && input_name="${input_name//.fa/}" && output_name="${input_name//*"/"/}" && bowtie2 -x "${input_name}" -1 cleaned_reads/"${output_name//.final.contigs/}".ok_1.fastq.gz -2 cleaned_reads/"${output_name//.final.contigs/}".ok_2.fastq.gz | samtools view -bS -o mags_binning/"${output_name}"_to_sort.bam' ::: mags_binning/*.fa* &> nohup_bt2_to_sort.out &
#Sort bam file
nohup parallel -j 36 'input_name={} && output_name="${input_name//_to_sort.bam/}" && output_name="${output_name##*/}" && samtools sort "${input_name}" -o mags_binning/"${output_name}"_sorted.bam' ::: mags_binning/*_to_sort.bam &> nohup_bt2_sorted.out &
#Generate index for bam file
nohup parallel -j 36 'input_name={} && output_name="${input_name//_sorted.bam/}" && output_name="${output_name##*/}" && samtools index "${input_name}"' ::: mags_binning/*_sorted.bam &> nohup_bt2_index.out &
#Run metabat for binning (no option for cores)
nohup parallel -j 36 'input_name={} && output_name="${input_name//.final.contigs.fa/}" && output_name="${output_name##*/}" && metabat -m 1500 -i "${input_name}" -o mags_binning/"${output_name}"_sorted.bam' ::: mags_binning/*.final.contigs.fa &> nohup_metabat.out &

#Copy raw MAGs to a new folder
mkdir raw_mags
cp mags_binning/*_sorted.bam.*.fa raw_mags/
#Rename raw MAGs
rename \sorted.bam. mag_ raw_mags/*

#Quality control of MAGs using checkm
#Run checkm standard lineage workflow
nohup checkm lineage_wf -t 36 -x fa raw_mags checkm_out/ &> nohup_checkm.out &

#Store mag quality results in a variable
#Get bin id, marker organism, completness and contamination from quality results
#cols 1,3,26,27,28,29
mag_quality=$(cut -f 1,3,26,27,28,29 -d " " checkm_out/storage/bin_stats_ext.tsv)

#Remove commas from assembly results to get numerical values
mag_quality="${mag_quality//,/}"
echo "$mag_quality" > raw_mags_stats.tsv

#Filter incomplete (< 50%) and contaminated (> 5%) mags
#Select mags showing completness higher than 50% (col 4)
completeness_contamination=$(echo "$mag_quality" | cut -f 1,4,6 -d " ")
completeness_contamination="${completeness_contamination//{\'marker /}"
echo "$completeness_contamination" | awk 'BEGIN {FS=" "}; $2>50' | awk 'BEGIN {FS=" "}; $3<5' > good_quality_mags_stats.tsv

#Copy good quality MAGs to a new folder
mkdir good_quality_mags

#Iterate quality MAGs list to copy these files to a new folder
while read line
do
mag_name=$(echo $line | awk '{print $1}') && cp raw_mags/"${mag_name}".fa good_quality_mags/ ;
done < good_quality_mags_stats.tsv

#Check MAGs's quality parameters
#Run quast
python quast-5.2.0/quast.py good_quality_mags/*.fa

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

#Functional and taxonomic annotation of MAGs

#Determine ORFs using prodigal to perform functional domain annotation
#Annotate sequences using prodigal
nohup parallel -j 36 'prodigal -f gff -a {/.}_aaORFs.fasta -i {} -o {/.}_ORFs.gff -d {/.}_ORFs.gff' ::: good_quality_mags/*.fa &> nohup_prodigal.out &
#Store results in different folders
mkdir prodigal_results
mv *ORFs* prodigal_results/

#Annotate general protein functional domains using HMMER and ORFs determined by prodigal
#HMMER database
#ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release
#HMMER analysis
nohup parallel -j 36 'hmmsearch --cpu 1 -E 1e-10 --tblout {/.}_pfam_annotation_tbl -o {/.}_pfam_annotation $HMMER_DBDIR/Pfam-A.hmm {}' ::: prodigal_results/*_aaORFs.fasta &> nohup_hmmsearch.out &
#Store hmmer outputs in a different folder
mkdir hmmer_results
mkdir hmmer_results/annotation
mkdir hmmer_results/tbl
mv *_pfam_annotation hmmer_results/annotation/
mv *_pfam_annotation_tbl hmmer_results/tbl/

#Additional HMMER analysis using NCBIFAM
nohup parallel -j 36 'hmmsearch --cpu 1 -E 1e-10 --tblout {/.}_ncbifam_annotation_tbl -o {/.}_ncbifam_annotation $NCBIFAM_DBDIR/ncbifam.hmm {}' ::: prodigal_results/*_aaORFs.fasta &> nohup_hmmsearch_ncbifam.out &
#Store hmmer outputs in a different folder
mkdir ncbifam_results
mkdir ncbifam_results/annotation
mkdir ncbifam_results/tbl
mv *_ncbifam_annotation ncbifam_results/annotation/
mv *_ncbifam_annotation_tbl ncbifam_results/tbl/

#Annotate CAZy functional domains present in MAGs
nohup parallel -j 20 'run_dbcan {} --db_dir $CAZY_DBDIR/cazy_db prok --out_dir {/.}' ::: good_quality_mags/*.fa &> nohup_cazy.out &
#Store cazy outputs in a different folder
mkdir cazy_results
mv *_mag_* cazy_results/
#Create folder to store functional domains results in csv
mkdir cazy_functional_domains
for folder in cazy_results/*/ ; do (cd "$folder" && code=$(basename "$PWD") && mv hmmer.out "${code}"_hmmer.csv); done
cp cazy_results/*/*_hmmer.csv cazy_functional_domains

#Taxonomic annotation of MAGs
gtdbtk classify_wf --genome_dir good_quality_mags --out_dir gtdbtk_output --cpus 36 --extension fa

#Inspect taxonomic identification of MAGs from bacteria
mag_bacteria_tax=$(cut -f1-2 gtdbtk_output/gtdbtk.bac120.summary.tsv)
#Format lines where taxonomic identification at species level was not reach (-r should be indicated to use .?)
mag_bacteria_tax=$(echo "$mag_bacteria_tax" | sed -r 's|;.?__$||g')
#Show taxonomic identification of each MAG at the corresponding level
mag_bacteria_tax=$(echo "$mag_bacteria_tax" | sed 's|d__.*__||g')
#Create a list showing taxonomic identification of MAGs from bacteria
echo "$mag_bacteria_tax" > mags_bacteria_id.txt

#Inspect taxonomic identification of MAGs from archaea
mag_archaea_tax=$(cut -f1-2 gtdbtk_output/gtdbtk.ar53.summary.tsv)
#Format lines where taxonomic identification at species level was not reach (-r should be indicated to use .?)
mag_archaea_tax=$(echo "$mag_archaea_tax" | sed -r 's|;.?__$||g')
#Show taxonomic identification of each MAG at the corresponding level
mag_archaea_tax=$(echo "$mag_archaea_tax" | sed 's|d__.*__||g')
#Create a list showing taxonomic identification of MAGs from archaea
echo "$mag_archaea_tax" > mags_archaea_id.txt

#Run Bacteria_ncbi_codes.sh to get bacteria ncbi codes
nohup bash Bacteria_ncbi_codes.sh &> nohup_ncbi_bacteria.out &
#Concatenate temporary files
cat temp_*bacteria* > mags_bacteria_ncbi_codes.txt
rm -rf temp_*bacteria*

#Get the names of MAGs from bacteria
bacteria_mags=$(cat mags_bacteria_id* | awk '{print $1}' | sed 1d)
#Add file extension
bacteria_mags=$(echo "$bacteria_mags" | sed -e 's/$/.fa/')
#Add the directory containing MAGs from bacteria
bacteria_mags=$(echo "$bacteria_mags" | sed -e 's/^/good_quality_mags\//')

#Execute prokka annotation pipeline
#Annotate MAGs from bacteria
nohup parallel -j 36 'prokka --outdir {/.}_prokka_out {} --prefix {/.} --addgenes --kingdom Bacteria' ::: echo "$bacteria_mags" &> nohup_prokka_bacteria.out &
#Store prokka outputs in a different folder
mkdir prokka_results
mv *prokka_out prokka_results/
ls prokka_results/*/*.gbk | wc -l

#Add NCBI taxa id to gbk files corresponding to MAGs from bacteria
while read -r line; do (
	#Get MAG name
	bacteria_mag_name=$(echo $line | awk '{print $2}')
	#Get NCBI taxa it corresponding to each MAG
	bacteria_mag_ncbi=$(echo $line | awk '{print $6}')
	#Add NCBI taxa id to gbk file generated by prokka
	sed -i -e 's/strain="strain"/db_xref="taxon:'${bacteria_mag_ncbi}'"/g' prokka_results/*/"${bacteria_mag_name}".gbk 
); done < mags_bacteria_ncbi_codes*

#Store megahit results in a single folder
mkdir megahit_results
mv *megahit* megahit_results/

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

