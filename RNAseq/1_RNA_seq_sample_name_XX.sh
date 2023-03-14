#!/bin/bash
#SBATCH --job-name='RNA_seq_sample_name_XX'
#SBATCH --output='RNA_seq_sample_name_XX.log'
#SBATCH --time='12:00:00'
#SBATCH --mem-per-cpu=50G
#SBATCH --partition='shortq'
#SBATCH -m block
/bin/hostname

module load star/2.5.2b
module load FastQC
module load samtools/1.8



# path to the directory with raw (unmapped) bam files
raw_data_location=SPECIFY_raw_data_location

# path to the annotation file with sample names
sample_annotation=SPECIFY_sample_annotation

# path to the mapping results directory
processing_results_dir=SPECIFY_processing_results_dir

#specify sample_name
sample_name=sample_name_XX


flowcell=$(cat ${sample_annotation} | multigrep.sh -f 1 -w -p ${sample_name}|cut -f2)
lane=$(cat ${sample_annotation} | multigrep.sh -f 1 -w -p ${sample_name}|cut -f3)

genome=hg38


echo $sample_name


date 

############################################################################################
# Generate fastq file
############################################################################################

# Target to produce: ${processing_results_dir}/fastq/${sample_name}_R1.fastq

fastq_out=$(echo ${processing_results_dir}/fastq/${sample_name}_R1.fastq)

`samtools view ${raw_data_location}/${flowcell}_${lane}#${sample_name}_*.bam | awk -v fastq_out=${fastq_out} '{ print "@"$1"\n"$10"\n+\n"$11 > fastq_out; }'`


############################################################################################
# Generate fastQC report
############################################################################################


# Target to produce fastqc report

`fastqc --noextract --outdir ${processing_results_dir}/fastqc/ ${processing_results_dir}/fastq/${sample_name}_R1.fastq`


############################################################################################
# Trimming
############################################################################################

# Target to produce: `{out_dir}/fastq/${sample_name}_R1_trimmed.fastq`

`/cm/shared/apps/java/jdk/1.7.0_80/bin/java -Xmx60000m -jar /cm/shared/apps/trimmomatic/0.32/trimmomatic-0.32-epignome.jar SE -phred33 -threads 2 ${processing_results_dir}/fastq/${sample_name}_R1.fastq ${processing_results_dir}/fastq/${sample_name}_R1_trimmed.fastq HEADCROP:13 ILLUMINACLIP:/data/groups/lab_bock/shared/resources/adapters/epignome_adapters_2_add.fa:2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:30`

`fastqc --noextract --outdir ${processing_results_dir}/fastqc/ ${processing_results_dir}/fastq/${sample_name}_R1_trimmed.fastq`


############################################################################################
# STAR mapping
############################################################################################

# Target to produce: `${processing_results_dir}/star_mapping/${sample_name}_Aligned.sortedByCoord.out.bam`

`STAR --runThreadN 4 --genomeDir /data/groups/lab_winter/reference_files/indices/STAR/ --readFilesIn ${processing_results_dir}/fastq/${sample_name}_R1_trimmed.fastq  --outFileNamePrefix ${processing_results_dir}/star_mapping/${sample_name}_ --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outStd Log`

/cm/shared/apps/samtools/1.8/bin/samtools index ${processing_results_dir}/star_mapping/${sample_name}_Aligned.sortedByCoord.out.bam


############################################################################################
# Generate a table with read counts using htseq-count
############################################################################################

`htseq-count -m intersection-nonempty -s yes -f bam -r pos --additional-attr gene_name ${processing_results_dir}/star_mapping/${sample_name}_Aligned.sortedByCoord.out.bam /data/groups/lab_winter/reference_files/genomes/hg38.gtf > ${processing_results_dir}/star_mapping/${sample_name}_read_counts.txt`

`htseq-count -m intersection-nonempty -s yes -f bam -r pos --additional-attr gene_name ${processing_results_dir}/hisat2/${sample_name}.aln_sorted.bam /data/groups/lab_winter/reference_files/genomes/hg38.gtf > ${processing_results_dir}/hisat2/${sample_name}_read_counts.txt`



############################
# Generate bigwig files:
############################

`bamCoverage -b ${processing_results_dir}/star_mapping/${sample_name}_Aligned.sortedByCoord.out.bam -o ${processing_results_dir}/star_mapping/${sample_name}_Aligned.sortedByCoord.out.RPKM.bw --normalizeUsing RPKM`



date
