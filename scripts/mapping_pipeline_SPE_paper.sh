#!/bin/bash
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
# 10. Type: chmod +x ./scripts/mapping_pipelinre_SPE_paper. Press return.
# 11. Type: ./scripts//mapping_pipelinre_SPE_paper. Press return.
# 12. The last command will execute the program.

# 13. The script will run for a few hours to a couple of days, Depending on the size of your fastq files and the size of the genome you are working with. You will know it finished once the prompt is back (the $ sign) and the file Rplot.pdf is located in the output folder.

#exec 3>&1 4>&2
#trap 'exec 2>&4 1>&3' 0 1 2 3
#exec 1>log.out 2>&1

#input files
mut_files=fastq/mut*
wt_files=fastq/wt*

#names
line=EMS
mut=EMS_mut
wt=EMS_wt

#install programs bwa and samtools
cd programs/bwa-0.7.12
make clean 
make

cd ../samtools-0.1.19/
make clean
make

cd ../../


#reference input files that are necessary to run the prograns
fa=refs/genome.fa
knownsnps=refs/knownsnps.vcf #ftp://ftp.ensemblgenomes.org/pub/plants/release-31/vcf/arabidopsis_thaliana/arabidopsis_thaliana.vcf.gz
snpEffDB=TAIR10.29 #paste the snpEff annotated genome name

#reference files that need to generated
#creating .fai file
programs/samtools-0.1.19/samtools faidx $fa
#creating bwa index files
programs/bwa-0.7.12/bwa index -p genome.fa -a is $fa
mv genome* refs/

#generating dict file for GATK
java -Xmx2g -jar programs/picard-tools-1.119/CreateSequenceDictionary.jar R=$fa O=refs/genome.dict
#mapping w/ BWA
programs/bwa-0.7.12/bwa mem -t 2 -M $fa ${mut_files[*]} > output/$mut.sam &
programs/bwa-0.7.12/bwa mem -t 2 -M $fa ${wt_files[*]} > output/$wt.sam
wait


#due to old samtools version this step is probably necessary
programs/samtools-0.1.19/samtools view -bSh output/$mut.sam > output/$mut.bam &
programs/samtools-0.1.19/samtools view -bSh output/$wt.sam > output/$wt.bam
wait

rm -r output/*.sam

#this step is probably needed only when you have paired-end; in any case it should come before coordinate sorting (next step) on name-sorted files
programs/samtools-0.1.19/samtools fixmate output/$mut.bam output/$mut.fix.bam &
programs/samtools-0.1.19/samtools fixmate output/$wt.bam output/$wt.fix.bam
wait

#sort by coordinates
programs/samtools-0.1.19/samtools sort output/$mut.fix.bam output/$mut.sort &
programs/samtools-0.1.19/samtools sort output/$wt.fix.bam output/$wt.sort
wait

java -Xmx2g -jar programs/picard-tools-1.119/MarkDuplicates.jar INPUT=output/$mut.sort.bam OUTPUT=output/$mut.sort.md.bam METRICS_FILE=output/$mut.matrics.txt ASSUME_SORTED=true &
java -Xmx2g -jar programs/picard-tools-1.119/MarkDuplicates.jar INPUT=output/$wt.sort.bam OUTPUT=output/$wt.sort.md.bam METRICS_FILE=output/$wt.matrics.txt ASSUME_SORTED=true
wait

#this part is just to add header for further gatk tools
java -Xmx2g -jar programs/picard-tools-1.119/AddOrReplaceReadGroups.jar INPUT=output/$mut.sort.md.bam OUTPUT=output/$mut.sort.md.rg.bam RGLB=$mut RGPL=illumina RGSM=$mut RGPU=run1 SORT_ORDER=coordinate &
java -Xmx2g -jar programs/picard-tools-1.119/AddOrReplaceReadGroups.jar INPUT=output/$wt.sort.md.bam OUTPUT=output/$wt.sort.md.rg.bam RGLB=$wt RGPL=illumina RGSM=$wt RGPU=run1 SORT_ORDER=coordinate
wait

java -Xmx2g -jar programs/picard-tools-1.119/BuildBamIndex.jar INPUT=output/$mut.sort.md.rg.bam &
java -Xmx2g -jar programs/picard-tools-1.119/BuildBamIndex.jar INPUT=output/$wt.sort.md.rg.bam
wait


#Variant calling using GATK HC extra parameters
java -Xmx2g -jar programs/GenomeAnalysisTK.jar -T HaplotypeCaller -R $fa -I output/$mut.sort.md.rg.bam -I output/$wt.sort.md.rg.bam -o output/$line.hc.vcf -minReadsPerAlignStart 7 -gt_mode DISCOVERY -out_mode EMIT_ALL_SITES -writeFullFormat -stand_emit_conf 10 -stand_call_conf 10 -nct 2 -variant_index_type LINEAR -variant_index_parameter 128000

############prepering for R#########################
#Exclude indels from a VCF
java -Xmx2g -jar programs/GenomeAnalysisTK.jar -R $fa -T SelectVariants --variant output/$line.hc.vcf -o output/$line.selvars.vcf --selectTypeToInclude SNP

#now make it into a table
java -jar programs/GenomeAnalysisTK.jar -R $fa -T VariantsToTable -V output/$line.selvars.vcf -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF GQ -o output/$line.table


####################################################################################################################################################
########################################now let's find the best candidates##########################################################################
java -jar programs/snpEff/snpEff.jar -c programs/snpEff/snpEff.config $snpEffDB -s output/snpEff_summary.html output/$line.selvars.vcf > output/$line.se.vcf

#and finally, get only the ones that are wt or hets in the wt bulk
grep -v '^##' output/$line.se.vcf | awk 'BEGIN{FS=" "; OFS=" "} $1~/#CHROM/ || $10~/^1\/1/ && ($11~/^1\/0/ || $11~/^0\/0/ || $11~/^0\/1/) && length($4==1) && length($5)==1 && $1~/^[0-9]*$/ && /splice_acceptor_variant|splice_donor_variant|start_lost|stop_gained|missense_variant/ {$3=$7=""; print $0}' | sed 's/  */ /g' | awk '{split($9,a,":"); split(a[2],b,","); if (b[1]>b[2] || $1~/#CHROM/) print $0}' > output/$line.cands2.txt
awk 'FNR==NR{a[$1$2];next};!($1$2 in a) || $1~/#CHROM/' $knownsnps output/$line.cands2.txt > output/$line.cands3.txt


#getting things a bit more organized and only the relevant data from cands3
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "chr" "pos" "ref" "alt" "mutation_effect" "gene" "At_num" "CDS_change" "protein_change" "$mut.ref" "$mut.alt" "$wt.ref" "$wt.alt" > output/cands4.txt
awk 'BEGIN{OFS="\t"} NR>1 {split($6,a,"|");split($8,b,":"); split(b[2],c,","); split($9,d,":"); split(d[2],e,","); gsub("c.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, a[2], a[4], a[5], a[10], a[11], c[1], c[2], e[1], e[2]}' output/$line.cands3.txt | awk '$0!~/\./ && (($10+$11)>4) && (($12+$13)>4)' >> output/cands4.txt 
####################################################################################################################################################
####################################################################################################################################################

####################################################################################################################################################
####################################################################################################################################################
#this command will make it ready to run w/ R to produce Manhatten plot
printf "%s\t" "CHR" "POS" "REF" "ALT" "mut_GT" "mut.ref" "mut.alt" "mut.DP" "mut.GQ" "wt.GT" "wt.ref" "wt.alt" "wt.DP" "wt.GQ" > output/$line.plot.txt; printf "\n" >> output/$line.plot.txt
awk 'length($3)==1 && length($4)==1 && $1~/^[0-9]*$/ && $5~/^[AGCT]/ && $9~/^[AGCT]/ && $0 !~ /NA/ {gsub(/\,/, "\t"); print}' output/$line.table | awk '$6+$11>0 && $8>3 && $13>3' >> output/$line.plot.txt

#and finally, just get rid of known snps
awk 'FNR==NR{a[$1$2];next};!($1$2 in a)' $knownsnps output/$line.plot.txt > output/$line.plot.no_known_snps.txt

#get the all snps in SnpEff format
awk 'FNR==NR{a[$1$2];next};($1$2 in a)' output/$line.plot.no_known_snps.txt output/$line.se.vcf > output/$line.plot2.txt
awk '{$3=$7=""; print $0}' output/$line.plot2.txt | sed 's/  */ /g' > output/$line.plot3.txt
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "chr" "pos" "ref" "alt" "mutation_effect" "gene" "At_num" "CDS_change" "protein_change" "mut.ref" "mut.alt" "wt.ref" "wt.alt" > output/plot4.txt
awk 'BEGIN{OFS="\t"} {split($6,a,"|");split($8,b,":"); split(b[2],c,","); split($9,d,":"); split(d[2],e,","); gsub("c.", "", a[10]); gsub("p\\.", "", a[11]); print $1, $2, $3, $4, a[2], a[4], a[5], a[10], a[11], c[1], c[2], e[1], e[2]}' output/$line.plot3.txt >> output/plot4.txt

####################################################################################################################################################
####################################################################################################################################################
Rscript ./scripts/analysis3.R






