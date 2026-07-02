
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

#Genome assembly

#Paired-end genomic reads are in the "raw_reads" folder

#Define variables indicating different database folders
$HMMER_DBDIR
$NCBIFAM_DBDIR
$CAZY_DBDIR

#Download additional reference genomes
#Lactiplantibacillus pentosus DSM 20314
#https://www.dsmz.de/collection/catalogue/details/culture/DSM-20314
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001433755.1/
#Levilactobacillus brevis DSM 20054
#https://www.dsmz.de/collection/catalogue/details/culture/DSM-20054
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001433855.1/

#DNA quality filtering
mkdir clean_reads
for r1 in raw_reads/*_1.fastq*; do (r2="${r1//_1/_2}" && fastq-mcf -k 0 --qual-mean 20 -w 5 -q 20 -l 150 n/a $r1 $r2 -o filtered_"${r1//"raw_reads/"/}" -o filtered_"${r2//"raw_reads/"/}");done
mv filtered_* clean_reads/
 
#Genome assembly
nohup sh -c 'for r1 in clean_reads/*_1.fastq*; do (strain_name="${r1//*filtered_/}" && strain_name="${strain_name//_*/}" && r2="${r1//_1./_2.}" && spades.py -t 36 --disable-gzip-output -m 100 -k21,33,55,77,99,127 --isolate -1 $r1 -2 $r2 -o spades_out_"${strain_name}"); done' &> nohup_spades.out &

#Quality control using quast
#Store spades outputs in a different folder
mkdir spades_results
mv spades_out_* spades_results/
#Copy assembled contigs to a different folder
mkdir assembled_strain_genomes
for contigs in spades_results/*/contigs.fasta; do (sample="${contigs//*_/}" && sample="${sample//"/"*/}" && cp $contigs assembled_strain_genomes/"${sample}".fasta); done
#Run quast
python quast-5.2.0/quast.py assembled_strain_genomes/*.fasta

#Add DSM genomes to continue the analysis
mkdir dsm_assemblies
#Upload reference DSM assemblies 
mkdir assembled_strain_genomes_and_dsm
cp assembled_strain_genomes/*.fa* assembled_strain_genomes_and_dsm/
cp dsm_assemblies/* assembled_strain_genomes_and_dsm/

#Antibiotic resistance and mobile element analysis
#Copy assembly genomes to a new folder
mkdir tormes_input/
cp assembled_strain_genomes_and_dsm/* tormes_input/
for genome in tormes_input/*.f*; do (
    genome_path=$(realpath $genome)
    genome_name="${genome//*"/"/}" && genome_name="${genome_name//.f*/}"
    echo -e "$genome_name\tGENOME\t$genome_path" > temp_tormes_"${genome_name}"
); done
#Generate combined metadata file
echo -e "Samples\tRead1\tRead2" > temp_header.txt
cat temp_header.txt temp_tormes_* > tormes_metadata.txt 
#Remove temporary files generated to create tormes metadata file
rm -rf temp_*
#Antibiotic resistance and mobile element analysis using TORMES
nohup tormes --metadata tormes_metadata.txt --output tormes_results --threads 36 &> nohup_tormes.out &

#Taxonomic identification
nextflow run CAMI-challenge/CAMITAX -profile singularity --db /opt/home/carlos/databases/CAMITAX/ --i tormes_results/genomes/ --x fasta

#Determine ORFs using prodigal to perform functional domain annotation
for sequence in assembled_strain_genomes_and_dsm/*; do (sequence_name="${sequence//*"/"/}" && sequence_name="${sequence_name//.*/}" && prodigal -f gff -a "${sequence_name}"_aaORFs.fasta -i $sequence -o "${sequence_name}"_ORFs.gff -d "${sequence_name}"_ORFs.gff); done
#Store results in different folders
mkdir prodigal_results
mv *ORFs* prodigal_results/

#HMMER analysis
nohup parallel -j 6 'hmmsearch --cpu 1 -E 1e-10 --tblout {/.}_pfam_annotation_tbl -o {/.}_pfam_annotation $HMMER_DBDIR/Pfam-A.hmm {}' ::: prodigal_results/*_aaORFs.fasta &> nohup_hmmer.out &
#Store HMMER outputs in a different folder
mkdir hmmer_results
mkdir hmmer_results/annotation
mkdir hmmer_results/tbl
mv *_pfam_annotation hmmer_results/annotation/
mv *_pfam_annotation_tbl hmmer_results/tbl/

#Additional HMMER analysis with NCBIFAM
#Annotate general protein functional domains using HMMER and ORFs determined by prodigal
#Activate conda environment
conda activate run_dbcan
#HMMER database
#ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release
#/opt/home/carlos/databases/hmmer_db/Pfam-A.hmm
#mkdir /opt/home/carlos/databases/NCBIfam_db

#Additional HMMER analysis using NCBIFAM
nohup parallel -j 36 'hmmsearch --cpu 1 -E 1e-10 --tblout {/.}_ncbifam_annotation_tbl -o {/.}_ncbifam_annotation $NCBIFAM_DBDIR/ncbifam.hmm {}' ::: prodigal_results/*_aaORFs.fasta &> nohup_hmmsearch_ncbifam.out &
#Store hmmer outputs in a different folder
mkdir ncbifam_results
mkdir ncbifam_results/annotation
mkdir ncbifam_results/tbl
mv *_ncbifam_annotation ncbifam_results/annotation/
mv *_ncbifam_annotation_tbl ncbifam_results/tbl/

#Annotate carbohydrate active enzymes
nohup parallel -j 20 'run_dbcan {} --db_dir $CAZY_DBDIR/cazy_db prok --out_dir cazy_out_{/.}' ::: assembled_strain_genomes_and_dsm/* &> nohup_cazy.out &
#Store cazy outputs in a different folder
mkdir cazy_results
mv cazy_out_* cazy_results/
#Create folder to store functional domains results in csv
mkdir cazy_functional_domains
for folder in cazy_results/*/ ; do (cd "$folder" && code=$(basename "$PWD") && mv hmmer.out "${code}"_hmmer.csv); done
cp cazy_results/*/*_hmmer.csv cazy_functional_domains
#Check all files were correctly copied (n=5)
ls cazy_functional_domains/*_hmmer.csv | wc -l
#Clean names from beneficial bacteria functional domains - csv format incompatibilities between linux and windows
cd cazy_functional_domains
#Delete cariage return from functional domain filenames (i.e. Agathobaculum'$'\r''_GCA_000701665.csv)
for domains in ./*_hmmer.csv; do
    #https://unix.stackexchange.com/questions/60779/how-do-i-remove-carriage-returns-from-directory-names
    rename 's/\r//' $domains
	echo $domains
done
#Return to the main folder
cd ../

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

#Finish
