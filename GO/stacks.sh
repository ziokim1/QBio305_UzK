#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --partition=devel-rh7
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40gb
#SBATCH --time=1:00:00
#SBATCH --account=UniKoeln
#SBATCH --error=/home/group.kurse/qbcbp009/Group2/stacks-%x-%j.err
#SBATCH --output=/home/group.kurse/qbcbp009/Group2/stacks-%x-%j.out

module load stacks/2.65

# Input directory # don'T forget to change this directory where you have saved 
# "input.vcf" and "pop.txt" files

input_dir="/home/group.kurse/qbcbp009/Group2/"

# Output directory, don't forget to change this to your stacks output directory
output_dir="/home/group.kurse/qbcbp009/Group2/stacks"

populations -V "${input_dir}group_2_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz" --popmap "${input_dir}pop.txt" --fstats -t 10 --structure -O "${output_dir}"

