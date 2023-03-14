#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --job-name='sample_name_XX_CRISPR_screen'
#SBATCH --output='sample_name_XX_CRISPR_screen.log'
#SBATCH --mem-per-cpu=80G
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH -m block
/bin/hostname

module load Python/3.9.5-GCCcore-10.3.0-bare
module load  HTSlib/1.12-GCC-10.3.0
module load FastQC
module load BEDTools/2.30.0-GCC-10.3.0
module load bio/BamTools/2.5.2-GCC-10.3.0
module load bio/cutadapt/3.4-GCCcore-10.2.0



# path to the directory with raw (unmapped) bam files
raw_data_location=SPECIFY_raw_data_location

# path to the annotation file with sample names
sample_annotation=SPECIFY_sample_annotation

# path to the mapping results directory
processing_results_dir=SPECIFY_processing_results_dir

#sample_name=S_1_CM20_37_DMSO_S64262
sample_name=sample_name_XX


lentiGuide_adaptor=SPECIFY_lentiGuide_adaptor
path_to_reference_genome=SPECIFY_path_to_reference_genome

flowcell=$(cat ${sample_annotation} | grep -w ${sample_name}|cut -f2)
lane=$(cat ${sample_annotation} | grep -w ${sample_name}|cut -f3)


# *** run the job ***
date


#Target to produce: fq file
bamtools convert -in ${raw_data_location}/${flowcell}_${lane}#${sample_name}.bam -format fastq > ${processing_results_dir}/fastq/${sample_name}.fq


#Target to produce: FASTQC file
fastqc -o ${processing_results_dir}/FASTQC -t 4 --nogroup ${processing_results_dir}/fastq/${sample_name}.fq


### Removing adapter sequences before sgRNA insert:

#Target to produce: `fastq/5trim_DMSO.fq`

cutadapt -g ${lentiGuide_adaptor} -o ${processing_results_dir}/fastq/5trim_${sample_name}.fq --minimum-length=10 ${processing_results_dir}/fastq/${sample_name}.fq


#Target to produce: `fastq/3trim_DMSO.fq`

module load FASTX-Toolkit/0.0.14-GCC-10.3.0
module load Bowtie2/2.4.4-GCC-10.3.0

fastx_trimmer -l 20 -i ${processing_results_dir}/fastq/5trim_${sample_name}.fq -o ${processing_results_dir}/fastq/3trim_${sample_name}.fq


#Target to produce: `aligned/bowtie2_completed.flag`
bowtie2 -x ${path_to_reference_genome} -U ${processing_results_dir}/fastq/3trim_${sample_name}.fq -N1 -S ${processing_results_dir}/aligned/${sample_name}.sam --no-hd

touch ${processing_results_dir}/aligned/${sample_name}_bowtie2_completed.flag


#Target to produce: `counts/counts_completed.flag`
cut -f 3 ${processing_results_dir}/aligned/${sample_name}.sam | sort | uniq -c > ${processing_results_dir}/counts/${sample_name}.txt

touch ${processing_results_dir}/counts/${sample_name}.txt



date


