#########################################
# GSEA analysis
#########################################

# Generate ranking files

cell_line=KBM7

i=DESeq2_RNAseq_CooksT_FiltT.W7_vs_DMSO.meta_score.txt

# order genes based on their Log2FC
awk -F "\t" '{if(NR>1 && $3!="NA"){print $1 "\t" $3}}' ${i} |cut -f2 -d "_" | awk -F "\t" '{if($2!=""){print $1 "\t" $2}}'  |sort -k2 -gr > ${cell_line}_${i%.txt}.Log2FC.rnk

# order genes based on their sign Log10(adjP)
awk -F "\t" '{if(NR>1 && $8!="NA"){print $1 "\t" $8}}' ${i} |cut -f2 -d "_" | awk -F "\t" '{if($2!=""){print $1 "\t" $2}}'  |sort -k2 -gr > ${cell_line}_${i%.txt}.signLog10adjP.rnk

# order genes based on their meta-score
awk -F "\t" '{if(NR>1 && $10!="NA"){print $1 "\t" $10}}' ${i} |cut -f2 -d "_" | awk -F "\t" '{if($2!=""){print $1 "\t" $2}}'  |sort -k2 -gr > ${cell_line}_${i%.txt}.meta_score.rnk



# Run GSEA

#!/bin/bash
#SBATCH --job-name='GSEA_RNA_msigdb'
#SBATCH --output='GSEA_RNA_msigdb.log'
#SBATCH --mem=50G
#SBATCH --cpus-per-task=2
#SBATCH --time=3-12:00:00
#SBATCH --partition=longq
#SBATCH --qos=longq
#SBATCH -m block
/bin/hostname
date

for i in KBM7_DESeq2_RNAseq_*W7*rnk
do
echo $i

java -cp /research/lab_winter/users/himrichova/resources/software/gsea/gsea-3.0.jar -Xmx102400m -Xms51200m xtools.gsea.GseaPreranked -rpt_label ${i%.rnk}_Msigdb_v7.5.1_nperm1k -gmx /research/lab_winter/users/himrichova/resources/software/gsea/msigdb.v7.5.1.symbols.gmt -rnk ${i} -set_min 15 -set_max 5500 -plot_top_x 20 -nperm 1000

date

done
