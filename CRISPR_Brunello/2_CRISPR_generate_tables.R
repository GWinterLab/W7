rm(list=ls()) 

require(writexl)

# Set working directory - the directory including counts per each sample:
setwd("processing_results/counts/")
list.files()

# specify sample names
sample_names<-c("S_1_MC21_3_DMSO_S83453","S_6_MC21_3_W7_S83456")



# Specify sample name to be used as control (dmso)
dmso_sample_name="S_1_MC21_3_DMSO_S83453"

# load the STARS_chip file (list of all perturbations and gene names)

STARS_chip<-read.table("/Volumes/research/lab_winter/users/himrichova/resources/Brunello_CRISPR_screen/STARS_chip.txt", header=F, dec=".", fill = TRUE)
head(STARS_chip)
#> head(STARS_chip)
#V1   V2
#1 A1BG_CATCTTCTTTCACCTGAACG A1BG
#2 A1BG_CTCCGGGGAGAACTCCGGCG A1BG
#3 A1BG_TCTCCATGGTGCATCAGCAC A1BG

dim(STARS_chip)
#[1] 76441     2


# load the reference plasmid library file (from MJ_17_61 experiment)

MJ_library<-read.table("/Volumes/research/lab_winter/users/himrichova/resources/Brunello_CRISPR_screen/MJ_17_61_p216library_N705.txt", header=F, dec=".", fill = TRUE)
head(MJ_library)
#> head(MJ_library)
#V1                        V2
#1 2215 A1BG_CATCTTCTTTCACCTGAACG
#2 1678 A1BG_CTCCGGGGAGAACTCCGGCG
#3    4 A1BG_TCTCCATGGTGCATCAGCAC

dim(MJ_library)
#[1] 76340     2





# Merge STARS_chip with the library

STARS_chip_and_MJ_library<-merge(STARS_chip,MJ_library,by.x = "V1", by.y= "V2", all.x=TRUE)
head(STARS_chip_and_MJ_library)
dim(STARS_chip_and_MJ_library)
colnames(STARS_chip_and_MJ_library)<-c("STARS_chip","Gene","Library")




####################################################################################################################
# Merge the STARS_chip_Library file with each count table (per each sample), perform normalization
####################################################################################################################


all_samples_raw_reads_table<-STARS_chip_and_MJ_library
all_samples_norm_reads_table<-STARS_chip_and_MJ_library
all_samples_raw_reads_plus1_table<-STARS_chip_and_MJ_library

for(sample_id in sample_names)
{
  print(sample_id)
  A <- read.table(paste(sample_id,".txt",sep=""), header=F, dec=".", fill = TRUE) # read the table with counts
  head(A)
  A_assigned <- A[-c(which(A[,2]=="*")),] # remove the row with * instead of the perturbation name (i.e. unassigned reads)

  total_assigned_reads<-sum(A_assigned[,1]) # calculate total number of assigned reads in the sample
  colnames(A_assigned)<-c("raw_reads","perturbations")
  
  # merge the sample count table with the STARS_chip_Library
  Merge_STARS_lib_A <- merge(STARS_chip_and_MJ_library,A_assigned,by.x = "STARS_chip", by.y= "perturbations", all.x=TRUE)

  A_raw_reads<-Merge_STARS_lib_A[,4]
  A_raw_reads[is.na(A_raw_reads)]<-0 # change all NAs to 0
  A_raw_reads_plus1<-A_raw_reads+1 # add 1 to each number to avoid 0s and related dividing issues
  
  A_normalized <- A_raw_reads_plus1/(total_assigned_reads/1000000) # normalization: counts per million
  
  
  Merge_STARS_lib_A_plus1_normalized<-cbind(Merge_STARS_lib_A,A_raw_reads_plus1,A_normalized)

colnames(Merge_STARS_lib_A_plus1_normalized)<-c("STARS_chip","Gene","Library",paste(sample_id,"_raw_reads",sep=""),paste(sample_id,"_norm_reads_plus1",sep=""),paste(sample_id,"_norm_reads",sep=""))

  # Write a table per sample includin: STARS library, plasmid library, raw reads, raw plus1 reads, normalized reads
  write.table(Merge_STARS_lib_A_plus1_normalized,paste(sample_id,".normalized.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)
  write_xlsx(Merge_STARS_lib_A_plus1_normalized,paste(sample_id,".normalized.xlsx",sep=""))
  
  #extract only the normalized reads:
  A_norm_reads<-Merge_STARS_lib_A_plus1_normalized[,c(6)]
  assign(paste("norm_reads_",sample_id,sep=""),  A_norm_reads)

  A_raw_reads<-Merge_STARS_lib_A_plus1_normalized[,c(4)]
  A_raw_readsplus1<-Merge_STARS_lib_A_plus1_normalized[,c(5)]
  
  all_samples_raw_reads_table<-cbind(all_samples_raw_reads_table,A_raw_reads)
  all_samples_norm_reads_table<-cbind(all_samples_norm_reads_table,A_norm_reads)
  all_samples_raw_reads_plus1_table<-cbind(all_samples_raw_reads_plus1_table,A_raw_readsplus1)
  
  }


colnames(all_samples_raw_reads_table)<-c("STARS_chip","Gene","Library",sample_names)
colnames(all_samples_norm_reads_table)<-c("STARS_chip","Gene","Library",sample_names)
colnames(all_samples_raw_reads_plus1_table)<-c("STARS_chip","Gene","Library",sample_names)

write_xlsx(all_samples_raw_reads_table,"all_samples_raw_reads_table.xlsx")
write.table(all_samples_raw_reads_table,"all_samples_raw_reads_table.txt",quote=FALSE,sep="\t",row.names=F,col.names=T)


write_xlsx(all_samples_norm_reads_table,"all_samples_norm_reads_table.xlsx")
write.table(all_samples_norm_reads_table,"all_samples_norm_reads_table.txt",quote=FALSE,sep="\t",row.names=F,col.names=T)



####################################################################################################################
# Calculate files with log2FC (treatment/dmso) 
# and input file for STARS that includes values between 0-1 based on log2FC ranking
####################################################################################################################

number_of_rows<-nrow(all_samples_norm_reads_table) #76441

STARS_perturbations<-all_samples_norm_reads_table[,1]
Gene<-all_samples_norm_reads_table[,2]

Log2FC_DMSO<-STARS_perturbations
Log2FC_DMSO<-as.data.frame(Log2FC_DMSO)

STARS_input_Log2FC_DMSO_rank01<-STARS_perturbations
STARS_input_Log2FC_DMSO_rank01<-as.data.frame(STARS_input_Log2FC_DMSO_rank01)

for(sample_id in sample_names)
{
  print(sample_id)
  sample_log2FC_dmso<-log2(all_samples_norm_reads_table[,sample_id]/all_samples_norm_reads_table[,dmso_sample_name])
  assign(paste(sample_id,"_log2fc",sep=""), sample_log2FC_dmso)
  Log2FC_DMSO<-cbind(Log2FC_DMSO,sample_log2FC_dmso)
  
  
  sample_log2FC_rank<-rank(sample_log2FC_dmso)/number_of_rows
  assign(paste(sample_id,"_log2fc_rank",sep=""), sample_log2FC_rank)
  STARS_input_Log2FC_DMSO_rank01<-cbind(STARS_input_Log2FC_DMSO_rank01,sample_log2FC_rank)
  
  }

Log2FC_DMSO<-cbind(Log2FC_DMSO,Gene)
head(Log2FC_DMSO)
colnames(Log2FC_DMSO)<-c("STARS_chip",sample_names,"Gene")

head(Log2FC_DMSO)

write.table(Log2FC_DMSO,"all_samples_Log2FC_DMSO.txt",quote=FALSE,sep="\t",row.names=F,col.names=T)
write_xlsx(Log2FC_DMSO,"all_samples_Log2FC_DMSO.xlsx")



STARS_input_Log2FC_DMSO_rank01<-cbind(STARS_input_Log2FC_DMSO_rank01,Gene)
head(STARS_input_Log2FC_DMSO_rank01)
colnames(STARS_input_Log2FC_DMSO_rank01)<-c("STARS_chip",sample_names,"Gene")


write.table(STARS_input_Log2FC_DMSO_rank01,"STARS_input_Log2FC_DMSO_rank01.txt",quote=FALSE,sep="\t",row.names=F,col.names=T)

write_xlsx(STARS_input_Log2FC_DMSO_rank01,"STARS_input_Log2FC_DMSO_rank01.xlsx")



####################################################################################################################
# Calculate MEDIAN log2FC value per each gene
####################################################################################################################

MEDIAN_log2FC_all_samples<-aggregate(get(sample_names[1]) ~ Gene, data = Log2FC_DMSO, median)[,1]
MEDIAN_log2FC_all_samples<-as.data.frame(MEDIAN_log2FC_all_samples)

for(sample_id in sample_names)
{
  print(sample_id)
  sample_MEDIAN_log2FC <- aggregate(get(sample_id) ~ Gene, data = Log2FC_DMSO, median)
head(sample_MEDIAN_log2FC)
assign(paste(sample_id,"_MEDIAN_log2FC",sep=""), sample_MEDIAN_log2FC) 

MEDIAN_log2FC_all_samples<-data.frame(MEDIAN_log2FC_all_samples,sample_MEDIAN_log2FC)
}



head(MEDIAN_log2FC_all_samples)
dim(MEDIAN_log2FC_all_samples)
MEDIAN_log2FC_per_Gene<-MEDIAN_log2FC_all_samples[,c(1,3,5)] # specify number of columns according to the number of samples analyzed
head(MEDIAN_log2FC_per_Gene)
colnames(MEDIAN_log2FC_per_Gene)<-c("Gene",sample_names)

write.table(MEDIAN_log2FC_per_Gene,"MEDIAN_log2FC_all_samples.txt",quote=FALSE,sep="\t",row.names=F,col.names=T)
write_xlsx(MEDIAN_log2FC_per_Gene,"MEDIAN_log2FC_all_samples.xlsx")


