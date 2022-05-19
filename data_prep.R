library(openxlsx)
library(tidyverse)
library(HGNChelper)
library(ggplot2)
library(dplyr)


##reading/writing data cBioportal
mutations_df <- read.csv("https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/msk_met_2021/data_mutations.txt", sep = "\t", check.names = FALSE)
colnames(mutations_df)[17] <- "patient_id"
mutations_df <- mutations_df %>% dplyr::rename_all(make.names)

cnv_df <- read.csv("https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/msk_met_2021/data_cna.txt", sep = "\t", check.names = FALSE)
library(reshape2)
a <- melt(cnv_df)
colnames(a)[2] <- "patient_id"
colnames(a)[3] <- "CNA"
cnv_df <- a
cnv_df %>% dplyr::rename_all(make.names)

# clinical data - samples
clin_df <- read.csv("https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/msk_met_2021/data_clinical_sample.txt", sep = "\t", check.names = FALSE)
colnames(clin_df) <-  clin_df[4,]
clin_df <- clin_df[-(1:4), , drop = FALSE]


# clinical data - patients
clinp_df <- read.csv("https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/msk_met_2021/data_clinical_patient.txt", sep = "\t", check.names = FALSE)
colnames(clinp_df) <-  clinp_df[4,]
clinp_df <- clinp_df[-(1:4), , drop = FALSE]

##filter to only colon cancer samples
coad_read <- clin_df[clin_df$ONCOTREE_CODE %in% c("COAD", "READ"),]
##dist of primary and metastatic cancers
ggplot(coad_read) + geom_bar(aes(SAMPLE_TYPE), stat="count", fill = "#2986e2") + labs(title = "Distribution of patients by sample type", x = "Sample Type",  y = "Patient Count") +  
  theme(plot.title = element_text(size = 25, face = "italic"), axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), panel.background = element_blank())
                                                                                                                                                


## select only primary
coad_read_p <- coad_read[coad_read$SAMPLE_TYPE == "Primary",]


mutations_df <- mutations_df[mutations_df$patient_id %in% coad_read_p$SAMPLE_ID,]
cnv_df <- cnv_df[cnv_df$patient_id %in% coad_read_p$SAMPLE_ID,]

clinp_df <- clinp_df[clinp_df$PATIENT_ID %in% coad_read_p$PATIENT_ID, ]

# Load the cancer related list and limit to cancer-related genes 
cancer_list <-read.csv("C:/Users/susia/AUA/CMHealthcare/final_project/MolecularSubtypes/Data/version_2_mpt_cancer_gene.csv")
cancer_genes <- cancer_list %>% pull(gene)


write_csv(coad_read, file = "./data/clinical_raw.csv")
write_csv(coad_read_p, file = "./data/clinical_primary_raw.csv")
write_csv(mutations_df, file = "./data/mut_raw.csv")
write_csv(cnv_df, file = "./data/cnv_raw.csv")
write_csv(clinp_df, file = "./data/clinical_patients.csv")

# mutations_df <- read.csv("./data/mut_raw.csv")
# cnv_df <- read.csv("./data/cnv_raw.csv")
# coad_read <- read.csv("./data/clinical_raw.csv")/\
# coad_read_p <- read.csv("./data/clinical_primary_raw.csv")


##preprocessing

## remove non needed types of mutations
mut_type2remove = c("RNA", "Silent","Splice_Region","Splice_Site")
mutations_df <- mutations_df %>% 
    filter(!(Variant_Classification %in% mut_type2remove))


# Remove CNVs with copy number equal to 0 (LOSS OF HETEROZIGOUS)
cnv_df <- cnv_df  %>% filter(!(CNA == 0))
# keep only cnvs with high copy numbers (2 or -2 according cBioPortal definition)
cnv_df <- cnv_df  %>% filter(CNA == 2 | CNA == -2)


# summary statistics/plots for cnvs and mutations after filtering

# Map incorrectly labeled gene symbols to correct HUGO symbols using HGNChelper
mut_genes_orig <- unique(mutations_df$Hugo_Symbol)
mut_genes_map <- checkGeneSymbols(mut_genes_orig) %>%
  pull(Suggested.Symbol)

cnv_genes_orig <- unique(cnv_df$Hugo_Symbol)
cnv_genes_map <- checkGeneSymbols(cnv_genes_orig) %>%
  pull(Suggested.Symbol)

mut_map_df <- data.frame(old=mut_genes_orig, new=mut_genes_map)
mutations_df$Corrected_Hugo <- mut_map_df$new[match(mutations_df$Hugo_Symbol,mut_map_df$old)]

cnv_map_df <- data.frame(old=cnv_genes_orig, new=cnv_genes_map)
cnv_df$Corrected_Hugo <- cnv_map_df$new[match(cnv_df$Hugo_Symbol,cnv_map_df$old)]

mutations_df <- mutations_df %>%
  transmute(patient_id,
            Hugo_Symbol = Corrected_Hugo) %>%
  drop_na()

cnv_df <- cnv_df %>%
  transmute(patient_id,
            Hugo_Symbol = Corrected_Hugo,
            CNA) %>%
  drop_na()

# Remove samples which have mutations in DNA repair genes
dna_repair_genes <- c("PMS2", "MLH1", "MSH2", "MSH3", "MSH6")
to_remove <- mutations_df %>%
    filter(Hugo_Symbol %in% dna_repair_genes) %>%
    distinct(patient_id) %>%
    pull(patient_id)
  
mutations_df <- mutations_df %>%
    filter(!(patient_id %in% to_remove))
  
cnv_df <- cnv_df %>%
    filter(!(patient_id %in% to_remove))


# Switch from recording each mutation to recording presence of mutation in a gene
# i.e. where samples have multiple mutations in the same gene, we just keep the
# information that there is at least one
mutations_df <- mutations_df %>%
  distinct()

#Shouldn't have any effect for CNVs but just in case!
cnv_df <- cnv_df %>%
  distinct()



# filter out all the genes not in cancer related gene list
mutations_df <- mutations_df %>%
  filter(Hugo_Symbol %in% cancer_genes)
cnv_df <- cnv_df %>%
  filter(Hugo_Symbol %in% cancer_genes)



# # And extra bit of processing for CNVs to label and AMP or HOMDEL
cnv_mapper <- function(input_df){
  output_df <- input_df
  output_df$CNA <- lapply(input_df$CNA, function(input_string){
    if(as.numeric(input_string) < 0){
      output_string <- "HOMDEL"
    }
    else{
      output_string <- "AMP"
    }
    output_string
  })
  
  output_df %>%
    mutate(CNA = map_chr(CNA, toString))
}

cnv_df_ <- cnv_df %>%
  cnv_mapper() %>%
  unite('Hugo_Symbol_CNA', c(Hugo_Symbol, CNA), sep="_")


mut_df_onehot <- mutations_df %>%
  distinct() %>%
  add_column(present = 1) %>%
  pivot_wider(
    names_from = Hugo_Symbol,
    values_from = present,
    values_fill = 0
  ) %>%
  arrange(desc(patient_id))


# Convert the CNV dataset to onehot format
cnv_df_onehot <- cnv_df_ %>%
  distinct() %>%
  add_column(present = 1) %>%
  pivot_wider(
    names_from = Hugo_Symbol_CNA,
    values_from = present,
    values_fill = 0
  ) %>%
  arrange(desc(patient_id))


write_csv(mut_df_onehot, file = "./data/mut_onehot.csv")
write_csv(cnv_df_onehot, file = "./data/cnv_onehot.csv")


##merging cnv to mut

gene_list <- cancer_list
genes <- c()
for (i in colnames(cnv_df_onehot_)[-1]) {
  gene_full <- i
  gene <- strsplit(gene_full, "_")
  gene_type <- gene_list[gene_list$gene == gene[[1]][1],]$type
  if((gene_type == "oncogene") && (gene[[1]][2] == "AMP")) {
    genes <- c(genes, gene_full)
  } else if ((gene_type == "tsg") && (gene[[1]][2] == "HOMDEL")) {
    genes <- c(genes, gene_full)
  } 
}

cnv_filtered <- cnv_df_onehot[, c(colnames(cnv_df_onehot[1]),genes)]

combined <- merge(x = mut_df_onehot, y = cnv_filtered, by = "patient_id", all.x = TRUE) 
combined[is.na(combined)] <- 0   

colnames(combined) <- gsub("_AMP", "", colnames(combined))
colnames(combined) <- gsub("_HOMDEL", "", colnames(combined))

combined_summed <- t(apply((combined)[-1],1, function(x) tapply(x,colnames(combined)[-1],sum)))

df_combined <- cbind(combined[1], as.data.frame(combined_summed))

write.csv(df_combined, "./data/mut_cnv_onehot.csv", row.names = FALSE)


##EDA on clinical data

##age
age <-  read.csv("./data/Age_at_Sequencing.txt", sep = '\t')
unique(roundage$Age.at.Sequencing)
ggplot(age, aes(x=Age.at.Sequencing)) +
  geom_histogram(bins = 20, fill = "#2986e2") +  scale_x_continuous(breaks = seq(0, 100, 5)) + labs(title = "Distribution of patients by age", x = "Age Group",  y = "Patient Count") +  
  theme(plot.title = element_text(size = 25, face = "italic", hjust = 0.5), axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), panel.background = element_blank())


##mutation count
mut_count <-  read.csv("./data/Mutation_Count.txt", sep = '\t')
mut_count_ <- mut_count %>%  mutate(Mutation.Count = case_when(Mutation.Count > 15 ~ ">15",
                                                               Mutation.Count <=2 ~ "<=2",
                                                               TRUE ~ as.character(Mutation.Count))) %>%
                            mutate(Mutation.Count = factor(Mutation.Count, 
                                                       levels = c("<=2", "3", "4", "5", "6","7", "8", 
                                                                  "9","10", "11", "12", "13", "14", "15",
                                                                  ">15")))
ggplot(mut_count_, aes(x=Mutation.Count)) +
  geom_bar(fill = "#2986e2") + labs(title = "Distribution of patients by mutation counts", x = "Mutation Count",  y = "Patient Count") +  
  theme(plot.title = element_text(size = 25, face = "italic"), axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"), panel.background = element_blank())


##sex
sex <- read.csv("./data/Sex.txt", sep = '\t')
sex_ <- group_by(sex, Sex) %>% count() %>% mutate(prop = n/length(sex$Study.ID)*100)
ggplot(sex_, aes(x="", y=n, fill = Sex)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual("Sex", values = c("Female"= "#452075", "Male"= "#2986e2")) +
  geom_text(aes(label = paste0(round(sex_$prop,2), "%")),
            position = position_stack(vjust = 0.5), size = 7, colour = "white") +
  theme(legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) 

##race
race <- read.csv("./data/Race_Category.txt", sep = '\t')
race_ <- as.data.frame(group_by(race, Race.Category) %>% count() %>% mutate(prop = n/length(race$Study.ID)*100))
group.colors <-  c("Asian-far east/indian subcont"='#0911ac',  "White" = '#8C0EAB',   "Pt refused to answer" = '#34a2f8', "Black or african american"  = '#fbaa30', 
                   "Other" ='#ffef69', "Native hawaiian or pacific isl"= '#911eb4', "Native american-am ind/alaska" ='#46f0f0', "No value entered"  = '#f032e6', "Unknown"  ='#bcf60c')
ggplot(race_, aes(x=reorder(Race.Category, -n) , y = n, fill = Race.Category)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = group.colors) +
  labs( x = "Race",  y = "Patient Count") +
  theme(legend.title=element_text(size=15), 
         legend.text=element_text(size=15),
         axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_blank()) 


  





