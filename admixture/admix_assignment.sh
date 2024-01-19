#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10gb
#SBATCH --time=1:00:00
#SBATCH --account=UniKoeln
#SBATCH --output=/home/group.kurse/qbcbp009/Group2/admixture.out
#SBATCH --error=/home/group.kurse/qbcbp009/Group2/admixture.err

###
# Admixture Analysis
###

# Generate the input file in plink format
FILE=a.thaliana_admix
plink --double-id --vcf /home/group.kurse/qbcbp009/Group2/group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz --make-bed --out /home/group.kurse/qbcbp009/Group2/$FILE --allow-extra-chr

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' /home/group.kurse/qbcbp009/Group2/$FILE.bim > /home/group.kurse/qbcbp009/Group2/$FILE.bim.tmp
mv /home/group.kurse/qbcbp009/Group2/$FILE.bim.tmp /home/group.kurse/qbcbp009/Group2/$FILE.bim

#Now, we are ready to run ADMIXTURE. We will run it with cross-validation (the default is 5-fold CV, for higher, choose e.g. cv=10) and K=2.

#Let’s now run it in a for loop with K=2 to K=10 and direct the output into log files

for i in {1..10}
do
 admixture --cv /home/group.kurse/qbcbp009/Group2/$FILE.bed $i > /home/group.kurse/qbcbp009/Group2/log${i}.out
done

###
#To identify the best value of k clusters which is the value with lowest cross-validation error, we need to collect the cv errors. 
#Below are three different ways to extract the number of K and the CV error for each corresponding K. 
#Like we said at the start of the course, there are many ways to achieve the same thing in bioinformatics!
###

FILE=a.thaliana_admix
awk '/CV/ {print $3,$4}' /home/group.kurse/qbcbp009/Group2/*out | cut -c 4,7-20 > /home/group.kurse/qbcbp009/Group2/$FILE.cv.error
cat /home/group.kurse/qbcbp009/Group2/a.thaliana_admix.cv.error

grep -h CV /home/group.kurse/qbcbp009/Group2/log*.out>/home/group.kurse/qbcbp009/Group2/cross_validation.txt
