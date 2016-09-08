#pipeline for mapping EMS mutants
# require for mac (10.11.6) and Linux (centOS 6.7) Java 1.7 (7u79; http://www.oracle.com/technetwork/java/javase/downloads/jdk7-downloads-1880260.html)


# 1. Download and unpack the mapper package from where-do-we-host-it?
# 2. Place the mapper folder in your home directory.
# 3. Rename your fastq files as follow. For the mutant and WT bulk, the names should start with mut. and wt. respectively (note the dot); for single or paired-end you should then have R1 or R1 and R2, respectively and end with .fastq. For example, if your mutant bulk was sequence in a paired-end format and the WT as single-end you should rename the three files as follow: mut.R1.fastq, mut.R2.fastq and wt.R1.fastq
# 4. Place the renamed fastq files in the fastq folder within the mapper folder.
# 5. Decompress the relevant (your species) reference file in the ref folder. That will create a folder with genome.fa and knownsmps.vcf files. Transfer these two files to the ref folder (one directory up). If the species you are working with has no files in this folder you should be able to download the FASTA and knownsnsps (a VCF file) from Ensembl. For example, for corn, go to the following link: http://plants.ensembl.org/Zea_mays/Info/Index. The FASTA file is under the genome assembly section; download the toplevel DNA file. The VCF file is under the variation section at the same webpage (zea_mays.vcf.gz). Unpack the files in the ref folder and name them genome.fa for the FASTA file and knownsnsps.vcf for the VCF file. If your species does not have a VCF file, unpack the empty.vcf.tar.gz file that is already located in the ref folder and place it in the ref folder.
# 6. Open the folder programs and then the folder snpEFF. Open the snpEFF.database.xlxs file and find your species in the second column (e.g., Rice); most species will have more than one entry. Chose the latest annotation as shown in the first column. For example, for rice the latest annotation would be rice7. Copy this latest annotation name.
# 7. Open the file mapping_pipeline_SPE_paper located in the scripts folder within the mapper folder and paste the genome annotation name you just copied to replace the text “paste_the_snpEff_genome_annotation_here. Save the file.
# 8. Open the Terminal application
# 9. Type: cd ~/mapper. Press return.
# 10. Type: chmod +x ./scripts/mapping_pipeline_SPE_paper. Press return.
# 11. Type: ./scripts/mapping_pipeline_SPE_paper. Press return.
# 12. The last command will execute the program.

# 13. The script will run for a few hours to a couple of days, Depending on the size of your fastq files and the size of the genome you are working with. You will know it finished once the prompt is back (the $ sign) and the file Rplot.pdf is located in the output folder.