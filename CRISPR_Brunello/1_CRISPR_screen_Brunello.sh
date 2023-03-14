#!/bin/bash

##########################################################################################################
# Process the CRIPSR screen Brunello data
# - from raw bam files to count tables per each sample
##########################################################################################################


if [ $# -ne 5 ] ; then
   printf '\nUsage: %s path_to_reference_genome sample_annotation bsf_data_location processing_results_dir path_to_processing_data_code \n' "${0}";
   printf '	  path_to_reference_genome:	 Specify genome - path_to_reference_Brunello_genome.\n';
   printf '	  sample_annotation:         Specify location of the sample annotation file.\n';
   printf '	  bsf_data_location:         Specify location of the BSF data.\n';
   printf '	  processing_results_dir:    Specify location where a new directory for generating results should be created.\n';
   printf '	  path_to_processing_data_code:   Specify location of the processing data code (sample_name_XX_CRISPR_screen_processing_data.sh).\n';
   exit 1;
fi



module load Bowtie2/2.4.4-GCC-10.3.0
module load Sambamba/0.7.1
module load HTSlib/1.12-GCC-10.3.0
module load SAMtools/1.12-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0


path_to_reference_genome="${1}"
sample_annotation="${2}"
bsf_data_location="${3}"
processing_results_dir="${4}"
path_to_processing_data_code="${5}"



mkdir ${processing_results_dir}

mkdir -p ${processing_results_dir}/fastq
mkdir -p ${processing_results_dir}/FASTQC
mkdir -p ${processing_results_dir}/aligned
mkdir -p ${processing_results_dir}/counts


cd ${processing_results_dir}


export LC_ALL="en_US.UTF-8"
export LC_CTYPE="en_US.UTF-8"


# copy the original "processing data" script to the working directory
cp ${path_to_processing_data_code} .

# adjust the processing data script per each sample to run it in parallel
for sample_name in $( awk -F "\t" '{if(NR>1){print $1}}' $sample_annotation)
do

echo $sample_name
flowcell=$(cat ${sample_annotation} | grep -w  ${sample_name}|cut -f2)
lane=$(cat ${sample_annotation} | grep -w  ${sample_name}|cut -f3)

raw_data_location=${bsf_data_location}/${flowcell}/${flowcell}_${lane}_samples

lentiGuide_adaptor=$(cat ${sample_annotation} | grep -w  ${sample_name}|cut -f4)

sed 's/sample_name_XX/'${sample_name}'/g' sample_name_XX_CRISPR_screen_processing_data.sh | sed 's+SPECIFY_raw_data_location+'${raw_data_location}'+g' | sed 's+SPECIFY_sample_annotation+'${sample_annotation}'+g' | sed 's+SPECIFY_processing_results_dir+'${processing_results_dir}'+g' | sed 's+SPECIFY_lentiGuide_adaptor+'${lentiGuide_adaptor}'+g' | sed 's+SPECIFY_path_to_reference_genome+'${path_to_reference_genome}'+g'> ${sample_name}_CRISPR_screen_processing_data.sh
done


# submit job per each sample
for sample_name in $(awk -F "\t" '{if(NR>1){print $1}}' ${sample_annotation})
do
echo $sample_name
sbatch ${sample_name}_CRISPR_screen_processing_data.sh
done
