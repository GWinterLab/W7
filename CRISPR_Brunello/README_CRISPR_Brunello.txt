################################
# CRISRP Brunello - README
################################


##############################
# Prepare the following files:
##############################

1) Generate a tab separated annotation file where the following information will be included. 
Keep the same order of the columns. 
###

sample_name: names of the raw *.bam samples
flowcell: flowcell name provided by the sequencing facility
lane: lane number provided by the sequencing facility (flowcell and lane are then included in the directory where the raw sequenced samples are stored by the sequencing facility)
lentiGuide_adaptor: sequence of the lentiGuide adaptor that is supposed to be cut off
condition: condition/treatment
corresponding_control_sample: corresponding control sample to the treatment 


Example of the annotation table:

sample_name	flowcell	lane	lentiGuide_adaptor	condition	corresponding_control_sample	
S_1_CM20_37_DMSO_S64262	BSF_0767_HH3NGBBXY	8	CGAAACACCG	DMSO	S_1_CM20_37_DMSO_S64262	
S_3_CM20_37_T243_S64261	BSF_0767_HH3NGBBXY	8	CGAAACACCG	T243	S_1_CM20_37_DMSO_S64262	




# Please note that CeMM raw samples are located in /nobackup/lab_bsf/samples
# Specific sequenced samples are located in e.g.: /nobackup/lab_bsf/samples/BSF_0955_HL5LNBBXY/BSF_0955_HL5LNBBXY_3_samples/
# therefore:
bsf_data_location=/nobackup/lab_bsf/samples
flowcell= BSF_0955_HL5LNBBXY
lane=3



2) copy the "per sample processing data code" 
sample_name_XX_CRISPR_screen_processing_data.sh


3) copy the brunello reference genome, STARS_chip - file with all perturbations, control plasmid library
- directory called "brunello_index"
- STARS_chip.txt
- MJ_17_61_p216library_N705.txt

4) copy 1_CRISPR_screen_Brunello.sh script and make sure it is executable, if not make it executable by running the following:
chmod a+x 1_CRISPR_screen_Brunello.sh


##############################
# 1st part -  bash script
##############################

5) Run the 1_CRISPR_screen_Brunello.sh to process raw data

Usage: 1_CRISPR_screen_Brunello.sh path_to_reference_genome sample_annotation bsf_data_location processing_results_dir path_to_processing_data_code

Specify the following 5 paths:
	  sample_annotation bsf_data_location processing_results_dir path_to_processing_data_code 
	  path_to_reference_genome:	 Specify genome - path_to_reference_Brunello_genome.
	  sample_annotation:         Specify location of the sample annotation file.
	  bsf_data_location:         Specify location of the BSF data.
	  processing_results_dir:    Specify location where a new directory for generating results should be created.
	  path_to_processing_data_code:   Specify location of the processing data code (sample_name_XX_CRISPR_screen_processing_data.sh).

Examples:
sample_annotation=/research/lab_winter/himrichova/1_Projects/15_CRISPR_screens/7_CRISPR_Brunello_MC21_3_BRUNELLO_L4851/CRISPR_Brunello_MC21_3_sample_annotation.txt

processing_results_dir=/research/lab_winter/himrichova/1_Projects/15_CRISPR_screens/7_CRISPR_Brunello_MC21_3_BRUNELLO_L4851/processing_results

path_to_reference_genome=/research/lab_winter/users/himrichova/resources/Brunello_CRISPR_screen/brunello_index/brunello

bsf_data_location=/nobackup/lab_bsf/samples

path_to_processing_data_code=/research/lab_winter/users/himrichova/resources/Brunello_CRISPR_screen/sample_name_XX_CRISPR_screen_processing_data.sh


Run the code with the following command:
./1_CRISPR_screen_Brunello.sh ${path_to_reference_genome} ${sample_annotation} ${bsf_data_location} ${processing_results_dir} ${path_to_processing_data_code}


Required software/tools:

module load BEDTools/2.27.1-foss-2018b
module load bio/HTSlib/1.9-foss-2018b 
module load bio/Bowtie2/2.3.4.2-foss-2018b 
module load bio/Sambamba/0.7.1
module load Python/3.6.6-foss-2018b
module load FastQC
module load bio/BamTools/2.5.1-foss-2018b 
module load bio/cutadapt/2.1-foss-2018b-Python-3.6.6
module load SAMtools/1.9-foss-2018b
module load FASTX-Toolkit/0.0.14-foss-2016b


##############################
# 2nd part - R script
##############################

6) Use the R script "2_CRISPR_generate_tables.R" to generate tables of:
- raw reads
- normalized reads
- log2FC values (over DMSO)

Needed R packages: 
writexl

############################################################
# 3rd part - 3_CRISPR_screen_Brunello_STARS.sh
############################################################

7) STARS data analysis

- copy the STARS scripts from the provided "resources" - the directory called STARS_v1.3 or download from here 
https://portals.broadinstitute.org/gpp/public/software/stars

-> Use the commands in the 3_CRISPR_screen_Brunello_STARS.sh script

The following variables have to be adapted:

processing_results_dir=/research/lab_winter/himrichova/1_Projects/15_CRISPR_screens/7_CRISPR_Brunello_MC21_3_BRUNELLO_L4851/processing_results
path_to_STARS_code=/research/lab_winter/users/himrichova/resources/Brunello_CRISPR_screen/STARS_v1.3
path_to_STARS_chip=/research/lab_winter/users/himrichova/resources/Brunello_CRISPR_screen/STARS_chip.txt


7.1. Generate the null distribution using stars_null_v1.3.fixed.py
(can take around 6h for 1000 iterations)

7.2. Run STARS analysis using stars_v1.3.fixed.py
- generate results for the negative and/or positive directionality


Needed modules: 
module load Python/3.6.6-foss-2018b



Please note that the original stars code has to be fixed to be compatible with the python version Python/3.6.6-foss-2018b
=> Parentheses has to be added
- the following parentheses has to be added in stars_null_v1.3.py

print('Please enter a relevant direction; P for positive and N for negative')
print(rand)
print("Analyzing ...")


- the following parentheses has to be added in stars_v1.3.py

print('Please enter a relevant direction; P for positive and N for negative')
print('Analyzing ...')
print("No hit genes found for ")+c



