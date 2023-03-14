
path_to_STARS_code=/research/lab_winter/users/himrichova/resources/Brunello_CRISPR_screen/STARS_v1.3
path_to_STARS_chip=/research/lab_winter/users/himrichova/resources/Brunello_CRISPR_screen/STARS_chip.txt

sample_annotation=/research/lab_winter/users/himrichova/1_Projects/15_Molecular_Glues/7_CRISPR_Brunello_MC21_3_BRUNELLO_L4851/CRISPR_Brunello_MC21_3_sample_annotation.txt

processing_results_dir=/research/lab_winter/users/himrichova/1_Projects/15_Molecular_Glues/7_CRISPR_Brunello_MC21_3_BRUNELLO_L4851/processing_results

mkdir -p ${processing_results_dir}/STARS_analysis_rank01_1000i

cd ${processing_results_dir}/STARS_analysis_rank01_1000i



# Generate the null distribution using stars_null_v1.3.fixed.py

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name='STARS_Null_rank01'
#SBATCH --output='STARS_Null_rank01.log'
#SBATCH --mem-per-cpu=80G
#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH -m block
/bin/hostname

date
module load Python/3.9.5-GCCcore-10.3.0-bare

path_to_STARS_code=/research/lab_winter/users/himrichova/resources/Brunello_CRISPR_screen/STARS_v1.3
path_to_STARS_chip=/research/lab_winter/users/himrichova/resources/Brunello_CRISPR_screen/STARS_chip.txt

sample_annotation=/research/lab_winter/users/himrichova/1_Projects/15_Molecular_Glues/7_CRISPR_Brunello_MC21_3_BRUNELLO_L4851/CRISPR_Brunello_MC21_3_sample_annotation.txt

processing_results_dir=/research/lab_winter/users/himrichova/1_Projects/15_Molecular_Glues/7_CRISPR_Brunello_MC21_3_BRUNELLO_L4851/processing_results


python ${path_to_STARS_code}/stars_null_v1.3.fixed.py --input-file ${processing_results_dir}/counts/STARS_input_Log2FC_DMSO_rank01.txt --chip-file ${path_to_STARS_chip} --thr 10 --num-ite 1000 --use-first-pert N
date









# Run STARS analysis using stars_v1.3.fixed.py



# results for the negative directionality
python  ${path_to_STARS_code}/stars_v1.3.fixed.py --input-file ${processing_results_dir}/counts/STARS_input_Log2FC_DMSO_rank01.txt --chip-file ${path_to_STARS_chip} --thr 10 --dir N --null Null_STARSOutput8_10.txt --use-first-pert N

# results for the positive directionality
python  ${path_to_STARS_code}/stars_v1.3.fixed.py --input-file ${processing_results_dir}/counts/STARS_input_Log2FC_DMSO_rank01.txt --chip-file ${path_to_STARS_chip} --thr 10 --dir P --null Null_STARSOutput8_10.txt --use-first-pert N


