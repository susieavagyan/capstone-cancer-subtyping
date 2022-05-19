###Spectral Clustering
library(openxlsx)
library(readr)
library(hash)
library(dplyr)
library(reshape)
library(tidyr)
library(ggplot2)

mut_cnv_df <- read_csv("./data/mut_cnv_onehot.csv", col_types=c("patient_id"="c")) %>%
  mutate(patient_id = as.character(patient_id))

mut_cnv_df <- mut_cnv_df[rowSums(mut_cnv_df[, -1])>2,] ##remove patients with 0 or 1 mutations



##weighting by the frequency of mutation
mut_weights <- read.csv("./data/Mutated_Genes.txt", sep = "\t") %>% filter(Freq != "<0.1%")
colnames(mut_weights)[1] <- "gene"
mut_weights$weighted_score <- parse_number(mut_weights$Freq)/100


mutations_colnames <- mut_cnv_df$patient_id

mutations_mat <- mut_cnv_df %>%
  select(!c(patient_id)) %>%
  as.matrix() %>%
  t()

colnames(mutations_mat) <- mutations_colnames

##similarity matrix
similarities <- Spectrum::CNN_kernel(mutations_mat, NN = 4, NN2 = 6)
drop <- rowSums(similarities)
drop <- names(drop[is.na(drop)])
similarities <- similarities[!rownames(similarities) %in% drop, ]
similarities <- similarities[, !colnames(similarities) %in% drop ]

r <- Spectrum::cluster_similarity(similarities, k = 5, 
                                  clusteralg = "GMM")


cluster_assignments <- tibble("patient_id" = colnames(similarities), "cluster" = r)

# First add a column for cluster assignments to the base df
df_with_clusters <- mut_cnv_df %>%
    inner_join(cluster_assignments, by = "patient_id")
  

num_members <- length(unique(df_with_clusters %>% pull(!!as.name("patient_id"))))
  
# Calculate cluster size, record members, calculate frequency of variants within each cluster and in the 
# global population, set up "stacks" for plotting, do Fisher's exact test to compare counts of variants in
# a cluster with outside the cluster

cluster_counts <- df_with_clusters %>%
  pivot_longer(!c(cluster, !!as.name("patient_id")), names_to = "gene", values_to = "count") %>%
  ungroup() %>%
  filter(count != 0)


cluster_stats <- cluster_counts %>%
  select(-c(count)) %>%
  unique() %>%
  group_by(cluster, gene) %>%
  mutate(cluster_freq = n()) %>%
  ungroup() %>%
  group_by(cluster) %>%
  mutate(clust_size = n_distinct(!!as.name("patient_id")), cluster_members = list(unique(!!as.name("patient_id")))) %>%
  select(-c(!!as.name("patient_id"))) %>%
  unique() %>%
  pivot_wider(names_from = "gene", values_from = "cluster_freq", values_fill = 0) %>%
  pivot_longer(!c(cluster, clust_size, cluster_members), names_to = "gene", values_to = "cluster_freq") %>%
  ungroup() %>%
  group_by(gene) %>%
  mutate(total_freq = sum(cluster_freq)) %>%
  ungroup() %>%
  mutate(ooc_freq = total_freq - cluster_freq) %>%
  group_by(gene, cluster) %>%
  mutate(cluster_prop = cluster_freq / clust_size, total_prop = total_freq / num_members, ooc_prop = ooc_freq / (num_members - clust_size),
         clust_tag = paste(cluster, ": (n = ", clust_size, ")", sep = "")) %>%
  unique() %>%
  ungroup() %>%
  rowwise() %>%
  mutate(fisher_p_value = stats::fisher.test(x = rbind(c(cluster_freq, ooc_freq), c(clust_size-cluster_freq, num_members - clust_size - ooc_freq)))$p.value) %>%
  mutate(stack_3 = ifelse(cluster_prop < ooc_prop, cluster_prop, ooc_prop),
         stack_2 = ifelse(cluster_prop < ooc_prop, ooc_prop - cluster_prop, 0),
         stack_1 = ifelse(cluster_prop < ooc_prop, 0, cluster_prop - ooc_prop)) %>%
  arrange(desc(total_freq))

# create the list of gene mutations and cnvs frequencies
most_frequent_genes <- cluster_stats %>% ungroup() %>%
  select(gene, total_freq) %>% unique()  %>% arrange(desc(total_freq))
top_most_frequent_genes <- most_frequent_genes %>% slice_head(n = 10)
selected_genes <- top_most_frequent_genes$gene


# Save inputs and outputs of this run of the clustering algorithm.
# create output file name
output_name <- "Spectrum_clust.RData"


save_dir <- file.path("./results/spectrum")

print(save_dir)

if (!dir.exists(save_dir)){
  dir.create(save_dir, recursive = TRUE)
} 

save(similarities,
     cluster_assignments,
     cluster_stats,
     mut_weights,
     file = file.path(save_dir, output_name))

write.csv(cluster_assignments, file = "./results/spectrum/spec_cluster_assignments_5.csv", row.names = F)

plot_cluster_stats <- function(input_omics = cluster_stats, num_genes = 50, gene_weights = mut_weights, clusters = NULL){
  

##Summary statistics  
  # Get the list of most frequently mutated genes in the dataset
  top_genes_list <- input_omics %>%
                               select(gene, total_freq) %>%
                               unique() %>%
                               arrange(desc(total_freq)) %>%
                               slice_head(n=num_genes) %>%
                               pull(gene)
                           
  
  # If gene weights are provided, pull genes with top `num_genes` weights
  if(!(length(gene_weights) == 0)){
    top_genes_list_weight <- gene_weights %>%
      arrange(desc(weighted_score)) %>%
      slice_head(n=num_genes) %>%
      pull(gene)
    
    # Replace the mutations to plot with this list, or join with the most frequent
    # Currently replacing
    #top_genes_list[[1]] <- union(top_genes_list_weight, top_genes_list[[1]])
    top_genes_list <- top_genes_list_weight
    
  }
  
  # Setting up "stacks" for plotting bar charts
  omics_plot <-   input_omics %>%
      pivot_longer(c(stack_1, stack_2, stack_3), names_to = "stack", values_to = "stack_height") %>%
      filter(gene %in% top_genes_list) %>%
      arrange(desc(cluster_prop))
  
  
  
  if(!is.null(clusters)){
    omics_plot <- omics_plot %>%
      filter(cluster %in% clusters)
  }
  
  omics_plot_summary <- omics_plot %>%
    select(cluster, clust_size, gene, cluster_freq, total_freq, ooc_freq, cluster_prop, total_prop, ooc_prop, fisher_p_value) %>%
    unique() %>%
    mutate(sig_1 = fisher_p_value < 0.01, sig_2 = fisher_p_value < 0.05) %>%
    mutate(significance = sig_1 + sig_2, pos_neg = ifelse(cluster_prop>ooc_prop, "+", "-")) %>%
    select(!c(sig_1, sig_2)) %>%
    filter(significance > 0) %>%
    group_by(cluster, significance) %>%
    #filter(significance == 2) %>%
    # ungroup() %>%
    arrange(cluster) %>%
    print() %>%
    select(cluster, gene, pos_neg, clust_size, cluster_prop, ooc_prop, fisher_p_value) %>%
    group_by(cluster) %>%
    unite("gene", c(gene, pos_neg), sep = "") %>%
    summarise(cluster_summary = paste(gene, collapse="/"), cluster_size = clust_size) %>%
    unique()
  
  # Plotting "stacks" to show gain or loss of mutation frequency compared with samples outside the cluster
  omics_plot <- omics_plot %>%
    ggplot(aes(stack_height, tidytext::reorder_within(gene, total_prop, cluster), fill = stack)) + #, fill = gene)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_crossbar(stat = "identity", aes(x=ooc_prop, 
                                         y=tidytext::reorder_within(gene, total_prop, cluster),
                                         ymin=tidytext::reorder_within(gene, total_prop, cluster), 
                                         ymax=tidytext::reorder_within(gene, total_prop, cluster)),
                  width=0.002, 
                  colour="black",
                  size=1.45,
                  orientation="x",
                  show.legend=F) +
    facet_wrap(~clust_tag, scales = "free_y", ncol = 2) +
    ggtitle("Cancer genes' mutation frequency in different clusters") +
    xlab("Proportion of samples in cluster carrying mutation") +
    ylab("Gene") +
    theme_bw() +
    tidytext::scale_y_reordered() + 
    theme(legend.position="none") +
    geom_text(aes(label=ifelse(fisher_p_value<0.01 & stack_height > 0 & cluster_prop > total_prop & stack != "stack_3", "**", "")), position=position_stack(vjust=0.5), color = "black") +
    geom_text(aes(label=ifelse(fisher_p_value<0.05 & fisher_p_value >= 0.01 & cluster_prop > total_prop & stack_height > 0 & stack != "stack_3", "*", "")), position=position_stack(vjust=0), color = "black") +
    scale_fill_manual(values = c("green", rgb(0.6,0.2,0.2, alpha=0.1), rgb(0,0,0, alpha=0.3))) +
    xlim(-0.01,1.01)
  
  # Return both the summary and the plot itself 
  list(omics_plot_summary, omics_plot)
  
}

cluster_stats_plot_list <- plot_cluster_stats(cluster_stats, num_genes = 50, gene_weights = mut_weights)

pdf(file = "./results/spectrum/cluster_stats_plot.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 20)
cluster_stats_plot_list[[2]]
dev.off()

write.csv(cluster_stats_plot_list[[1]], file = "./results/spectrum/cluster_stats.csv", row.names = F)


##calculating silhouette score

library(cluster)
cluster_assignments <- read.csv("./results/spectrum/spec_cluster_assignments_5.csv")
assignment_vector <- cluster_assignments$cluster
names(assignment_vector) <- cluster_assignments$patient_id
si <- silhouette(assignment_vector, similarities)
summary(si)

