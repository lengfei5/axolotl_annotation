#!/bin/bash
 
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --job-name="Blast"
#SBATCH --output="slurm-%x_%j.out"
#SBATCH --error="slurm-%x_%j.err"
#SBATCH --qos=short
#SBATCH --partition=m
 
module load blast+/2.8.1-foss-2018b
filename='results/AmHoxCs_Predicted_AK'

blastn -db /groups/tanaka/Databases/Blast/Am_genome/AmexG_v6.0.DD/AmexG_v6.DD -query ${filename}.fa -out ${filename}.outfmt6 -outfmt 6 -max_target_seqs 5 -num_threads 16


