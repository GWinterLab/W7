##################################################################
# RNA_seq data analysis
##################################################################


genome=hg38
sample_annotation=RNAseq_sample_annotation.txt 
processing_results_dir=mapping_results
bsf_data_location=/scratch/lab_bsf/samples

path_to_processing_data_code=/data/groups/lab_winter/himrichova/resources/scripts/15_CRISPR_screens/RNAseq/1_RNA_seq_sample_name_XX.sh




######################################################################
# Generate directories for your RNA-seq data analysis results
######################################################################


mkdir -p ${processing_results_dir}
mkdir -p ${processing_results_dir}/fastq
mkdir -p ${processing_results_dir}/fastqc
mkdir -p ${processing_results_dir}/star_mapping
mkdir -p ${processing_results_dir}/hisat2




cd ${processing_results_dir}



# copy the original "processing data" script to the working directory
cp ${path_to_processing_data_code} .


################################################################################################
# Adjust the processing data script per each sample to run it in parallel
################################################################################################

for sample_name in $( awk -F "\t" '{if(NR>1){print $1}}' $sample_annotation)
do

echo $sample_name

flowcell=$(cat ${sample_annotation} | multigrep.sh -f 1 -w -p ${sample_name}|cut -f2)
lane=$(cat ${sample_annotation} | multigrep.sh -f 1 -w -p ${sample_name}|cut -f3)

raw_data_location=${bsf_data_location}/${flowcell}/${flowcell}_${lane}_samples

sed 's/sample_name_XX/'${sample_name}'/g' 0_RNA_seq_sample_name_XX.sh | sed 's+SPECIFY_raw_data_location+'${raw_data_location}'+g' | sed 's+SPECIFY_sample_annotation+'${sample_annotation}'+g' | sed 's+SPECIFY_processing_results_dir+'${processing_results_dir}'+g' > ${sample_name}_RNA_seq_processing_data.sh
done



################################
# Submit job per each sample
################################
for sample_name in $(awk -F "\t" '{if(NR>1){print $1}}' ${sample_annotation})
do
echo $sample_name
sbatch ${sample_name}_RNA_seq_processing_data.sh
done

