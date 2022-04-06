## run liftOver python script
python ../scripts/liftOver_axolotl.py --chr /groups/tanaka/People/current/jiwang/scripts/schloissnig_axolotl/data/ambMex60DD.scaffolds_structure.txt.gz --sizes /groups/tanaka/People/current/jiwang/scripts/schloissnig_axolotl/data/ambMex60DD.contigs.chrom.sizes --bed Clusters_HoxAs_contigs_coordinate.bed > HoxA_cluster_chr.bed

python ../scripts/liftOver_axolotl.py --chr /groups/tanaka/People/current/jiwang/scripts/schloissnig_axolotl/data/ambMex60DD.scaffolds_structure.txt.gz --sizes /groups/tanaka/People/current/jiwang/scripts/schloissnig_axolotl/data/ambMex60DD.contigs.chrom.sizes --bed Clusters_HoxBs_contigs_coordinate.bed > HoxB_cluster_chr.bed

python ../scripts/liftOver_axolotl.py --chr /groups/tanaka/People/current/jiwang/scripts/schloissnig_axolotl/data/ambMex60DD.scaffolds_structure.txt.gz --sizes /groups/tanaka/People/current/jiwang/scripts/schloissnig_axolotl/data/ambMex60DD.contigs.chrom.sizes --bed Clusters_HoxCs_contigs_coordinate.bed > HoxC_cluster_chr.bed

