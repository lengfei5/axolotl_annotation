#!/bin/bash
 
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --job-name="Blast"
#SBATCH --output="slurm-%x_%j.out"
#SBATCH --error="slurm-%x_%j.err"
#SBATCH --qos=short
#SBATCH --partition=c
 
module load blast+/2.8.1-foss-2018b
 
blastn -db /groups/tanaka/Databases/Blast/Am_genome/AmexG_v6.0.DD/AmexG_v6.DD -query “filename”.fa -out “filename”.outfmt6 -outfmt 6 -max_target_seqs 5 -num_threads 8


