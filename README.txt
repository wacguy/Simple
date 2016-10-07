#pipeline for mapping EMS mutants

# runs on mac (10.11.6) or Linux (centOS 6.7)
#requires Java 1.7 (7u79; http://www.oracle.com/technetwork/java/javase/downloads/jdk7-downloads-1880260.html)
#requires R and the following packages: ggplot2 and ggrepel


# 1. Download and unpack the Simple package from the following link: https://github.com/wacguy/EMS_Simple by pressing the green link: "Clone or download"
# 2. Place the Simple folder in your home directory.
# 3. Download the GATK executable from the following link: https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.5-0-g36282e4, unpack it and copy the GenomeAnalysisTK.jar file to the programs directory.
# 4. Rename your fastq files as follow. For the mutant and WT bulk, the names should start with mut. and wt. respectively (note the dot); for single or paired-end you should then have R1 or R1 and R2, respectively and end with .fastq. For example, if your mutant bulk was sequence in a paired-end format and the WT as single-end you should rename the three files as follow: mut.R1.fastq, mut.R2.fastq and wt.R1.fastq.
# 5. Place the renamed fastq files in the fastq folder located in the Simple folder.
# 6. Open the folder scripts inside Simple; open the data_base.txt file. Locate your species in the first column and copy it. 
# 7. Open the file simple.sh inside the folder scripts with a text editor and paste the species name you've just copied to replace Arabidopsis_thaliana as the species name (e.g., this line should look like: my_species=Arabidopsis_thaliana or my_species=Oryza_sativa_Japonica)
# 8. Open the Terminal application.
# 9. Type: cd ~/Simple. Press return.
# 10. Type: chmod +x ./scripts/simple.sh. Press return.
# 11. Type: ./scripts/simple.sh. Press return.
# 12. The last command will execute the program.
# 13. The script will run for a few hours up to a couple of days, depending on the size of your fastq files and the size of the genome you are working with. You will know it finished once the prompt says “Simple is done”.

