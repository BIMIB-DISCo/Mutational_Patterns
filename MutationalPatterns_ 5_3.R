### MutationalPatterns_.R

library(MutationalPatterns)
library(BSgenome)
library(NMF)
library(Biostrings)
library(stringr)
library(stringdist)
library(Matching)
library(gplots)


ref_genome <- available.genomes()[59]
ref_genome <- available.genomes()[52] #64
print(ref_genome)

if(ref_genome %in% rownames(installed.packages()) == FALSE) 
{
	BiocManager::install(ref_genome, update = FALSE)
}
library(ref_genome, character.only = TRUE)




###
compress_context <- function(signature, flanking = 1) {
  #flanking <- 2
  stringlen <- 2*flanking+5
  rnames <- rownames(signature)
  stringlen_orig <- str_length(rnames[1])
  stringlen_rem <- (stringlen_orig-stringlen)/2
  if (stringlen_rem<0)
    stringlen_rem <- 0
  triseq <- unique(str_sub(rnames,stringlen_rem+1,stringlen_orig-stringlen_rem))
  df <- data.frame(matrix(ncol=length(colnames(signature)),nrow=0))
  colnames(df) <- colnames(signature)
  for (s in triseq)
  {
    str_pat <- paste0(".*",s,".*") 
    pat <- str_pat
    pat <- str_replace_all(pat,"\\[","\\\\[")
    pat <- str_replace_all(pat,"\\]","\\\\]")
    df[s,] <- colSums(signature[rnames[str_detect(rnames,pat)],])
  }
  return(df)
}


### data loading

fig_folder <- "./fig_3_vs_5/"
dir.create(fig_folder, showWarnings = FALSE)
vcf_folder <- "./vcf_organoids/"


vcf_files_20 <- list.files(path = vcf_folder, pattern = "test[^_][^B].*[^mM].vcf$", full.names = FALSE)
vcf_files_20

vcf_files_treated <- list.files(path = vcf_folder, pattern = "Ttest_full.[^B].*.vcf$", full.names = FALSE)
vcf_files_treated

vcf_files_par <- list.files(path = vcf_folder, pattern = "Ptest_full.[^B].*.vcf$", full.names = FALSE)
vcf_files_par

vcf_files <- c(vcf_files_20, vcf_files_par, vcf_files_treated)
vcf_files <- file.path(vcf_folder, vcf_files)
vcf_files

sample_names <- c(vcf_files_20, vcf_files_par, vcf_files_treated)
sample_names <- str_replace_all(str_replace_all(sample_names, ".*test.",""), "\\.+full|^full.|.vcf$", "")
sample_names

tissue <- c(rep("colon", length(c(vcf_files_20,vcf_files_par))), 
            rep("colon_treated", length(c(vcf_files_treated))))



grl_20 <- read_vcfs_as_granges(vcf_files[1:20], sample_names[1:20], ref_genome, type = "all", predefined_dbs_mbs = TRUE)

grl_4 <- read_vcfs_as_granges(vcf_files[21:24], sample_names[21:24], ref_genome, type = "all", predefined_dbs_mbs = TRUE)

#risoluzione e merging
grl_4b <- grl_4
seqinfo(grl_4b) <- Seqinfo(genome="GRCh38")
grl_4c <- trim(grl_4b)
grl <- c(grl_20,grl_4c)

### end of data loading



snv_grl <- get_mut_type(grl, type = "snv", predefined_dbs_mbs=TRUE)
#indel_grl <- get_mut_type(grl, type = "indel", predefined_dbs_mbs=TRUE)
#dbs_grl <- get_mut_type(grl, type = "dbs", predefined_dbs_mbs=TRUE)
#mbs_grl <- get_mut_type(grl, type = "mbs", predefined_dbs_mbs=TRUE)


type_context5 <- type_context(snv_grl[[1]], ref_genome, extension=2)
lapply(type_context5, head, 12)


type_context3 <- type_context(snv_grl[[1]], ref_genome, extension=1)
lapply(type_context3, head, 12)


type_occurrences <- mut_type_occurrences(snv_grl, ref_genome)
type_occurrences

i=0

p=plot_spectrum(type_occurrences)
pdf(paste0(fig_folder, i, ".pdf"), height=15, width=17)
p
dev.off() 
i=i+1


p=plot_spectrum(type_occurrences, CT = FALSE, 
              indv_points = TRUE, legend = TRUE)
pdf(paste0(fig_folder, i, ".pdf"), height=15, width=17)
p
dev.off() 
i=i+1


mut_mat_3 <- mut_matrix(snv_grl, ref_genome, extension = 1)


plot_profile_heatmap(mut_mat_3, by = sample_names)
p=plot_profile_heatmap(mut_mat_3, by = sample_names)

pdf(paste0(fig_folder, i, ".pdf"), height=65, width=17)
p
dev.off() 
i=i+1

p=plot_96_profile(mut_mat_3)
pdf(paste0(fig_folder, i, ".pdf"), height=55, width=17)
p
dev.off() 
i=i+1


mut_mat_5 <- mut_matrix(snv_grl, ref_genome, extension = 2)
head(mut_mat_5)

plot_profile_heatmap(mut_mat_5, by = sample_names)
p=plot_profile_heatmap(mut_mat_5, by = sample_names)
pdf(paste0(fig_folder, i, ".pdf"), height=65, width=17)
p
dev.off() 
i=i+1


p=plot_profile_heatmap(mut_mat_5, by = tissue)
pdf(paste0(fig_folder, i, ".pdf"), height=6, width=17)
p
dev.off() 
i=i+1


mut_mat_5_to_3 <- compress_context(mut_mat_5, flanking = 1)
print(sum(mut_mat_5_to_3-mut_mat_3))

p=plot_96_profile(mut_mat_5_to_3)
pdf(paste0(fig_folder, i, ".pdf"), height=55, width=17)
p
dev.off() 
i=i+1



mut_mat2 <- mut_mat_5 #+ 0.0001 #operation already in the code


estimate_5 <- nmf(mut_mat_5, rank = 2:4, method = "brunet", nrun = 200, seed = 123456, .opt = "v-p")

plot(estimate_5)
p=plot(estimate_5)
pdf(paste0(fig_folder, i, ".pdf"), height=5, width=6)
p
dev.off() 
i=i+1


signatures = get_known_signatures()


# 3 denovo signatures using 5 bases
nmf_res_5_3S <- extract_signatures(mut_mat_5, rank = 3 , nrun = 200, single_core = TRUE)
colnames(nmf_res_5_3S$signatures) <- paste("Signature",seq(1:3))
rownames(nmf_res_5_3S$contribution) <- colnames(nmf_res_5_3S$signatures)



p=plot_profile_heatmap(nmf_res_5_3S$signatures)
pdf(paste0(fig_folder, i, ".pdf"), height=6, width=17)
p
dev.off() 
i=i+1


nmf_res_5_3S_to_3 <- compress_context(nmf_res_5_3S$signatures, flanking = 1)
nmf_res_5_3S_to_3_rec <- compress_context(nmf_res_5_3S$reconstructed, flanking = 1)

p=plot_96_profile(nmf_res_5_3S_to_3)
pdf(paste0(fig_folder, i, ".pdf"), height=55, width=17)
p
dev.off() 
i=i+1


p=plot_original_vs_reconstructed(mut_mat_5, nmf_res_5_3S$reconstructed, 
                                 y_intercept = 0.95)
pdf(paste0(fig_folder, i, ".pdf"), height=7, width=6)
p
dev.off() 
i=i+1


p=plot_original_vs_reconstructed(nmf_res_5_3S_to_3, nmf_res_5_3S_to_3_rec, 
                                 y_intercept = 0.95)
pdf(paste0(fig_folder, i, ".pdf"), height=7, width=6)
p
dev.off() 
i=i+1

#should this be done for each sample?
plot_compare_profiles(as.matrix(mut_mat_5_to_3)[, 1],
                      nmf_res_5_3S_to_3_rec[, 1],
                      profile_names = c("Original", "Reconstructed"),
                      #condensed = TRUE
)



rename_nmf_signatures(nmf_res_5_3S_to_3, signatures, cutoff = 0.85)
colnames(nmf_res$signatures)

# 2 denovo signatures using 5 bases
nmf_res_5_2S <- extract_signatures(mut_mat_5, rank = 2 , nrun = 200, single_core = TRUE)
colnames(nmf_res_5_2S$signatures) <- paste("Signature",seq(1:2))
rownames(nmf_res_5_2S$contribution) <- colnames(nmf_res_5_2S$signatures)



p=plot_profile_heatmap(nmf_res_5_2S$signatures)
pdf(paste0(fig_folder, i, ".pdf"), height=6, width=17)
p
dev.off() 
i=i+1


nmf_res_5_2S_to_3 <- compress_context(nmf_res_5_2S$signatures, flanking = 1)

p=plot_96_profile(nmf_res_5_2S_to_3)
pdf(paste0(fig_folder, i, ".pdf"), height=55, width=17)
p
dev.off() 
i=i+1

####




p=plot_contribution_heatmap(nmf_res$contribution, 
                          cluster_samples = TRUE, 
                          cluster_sigs = TRUE)
pdf(paste0("./fig/", i, ".pdf"), height=6.6, width=7)
p
dev.off() 
i=i+1



hclust_signatures <- cluster_signatures(nmf_res$signatures, method = "complete")
signatures_order <- colnames(nmf_res$signatures)[hclust_signatures$order]
signatures_order

hclust_samples <- cluster_signatures(mut_mat2, method = "complete")
samples_order <- colnames(mut_mat2)[hclust_samples$order]
samples_order

p=plot_contribution_heatmap(nmf_res$contribution,
                          sig_order = signatures_order, sample_order = samples_order,
                          cluster_sigs = FALSE, cluster_samples = FALSE
)
pdf(paste0("./fig/", i, ".pdf"), height=10, width=12)
p
dev.off() 
i=i+1


p=plot_original_vs_reconstructed(mut_mat2, nmf_res$reconstructed, 
                               y_intercept = 0.95)
pdf(paste0("./fig/", i, ".pdf"), height=7, width=6)
p
dev.off() 
i=i+1


cos_sim_signatures <- cos_sim_matrix(nmf_res$signatures, nmf_res$signatures)
p=plot_cosine_heatmap(cos_sim_signatures, 
                    cluster_rows = TRUE, cluster_cols = TRUE)
pdf(paste0("./fig/", i, ".pdf"), height=5, width=6)
p
dev.off() 
i=i+1



p=plot_cosine_heatmap(cos_sim_matrix(mut_mat2, mut_mat2), 
                    cluster_rows = TRUE, cluster_cols = TRUE)
pdf(paste0("./fig/", i, ".pdf"), height=10, width=12)
p
dev.off() 
i=i+1


signatures = get_known_signatures()




cutree(hclust_signatures, k = 4)

merged_signatures <- merge_signatures(nmf_res$signatures, cos_sim_cutoff = 0.82)
#merged_signature <- round(merged_signatures)
colnames(merged_signatures)

p=plot_cosine_heatmap(cos_sim_matrix(merged_signatures, merged_signatures), 
                    cluster_rows = TRUE, cluster_cols = TRUE)
pdf(paste0("./fig/", i, ".pdf"), height=10, width=12)
p
dev.off() 
i=i+1


fit_res_merged <- fit_to_signatures(mut_mat2, merged_signatures)


p=plot_contribution(fit_res_merged$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
)
pdf(paste0("./fig/", i, ".pdf"), height=7, width=8.8)
p
dev.off() 
i=i+1

p=plot_original_vs_reconstructed(mut_mat2, fit_res_merged$reconstructed, 
                               y_intercept = 0.95)
pdf(paste0("./fig/", i, ".pdf"), height=7, width=5)
p
dev.off() 
i=i+1

merged_sig_name <- paste("signature_m",seq(1:length(colnames(merged_signatures))))
merged_signatures_rnamed <- merged_signatures
colnames(merged_signatures_rnamed) <- merged_sig_name


fit_res_merged_std <- fit_res_merged
fit_res_merged_std$signatures <- merged_signatures_rnamed
fit_res_merged_std$signatures <- as.data.frame(fit_res_merged_std$signatures)
fit_res_merged_std$reconstructed <- as.data.frame(fit_res_merged_std$reconstructed)
rownames(fit_res_merged_std$contribution) <- colnames(fit_res_merged_std$signatures)

fit_res_merged_reduced <- fit_res_merged_std
fit_res_merged_reduced$signatures <- compress_context(fit_res_merged_reduced$signatures)
fit_res_merged_reduced$reconstructed <- compress_context(fit_res_merged_reduced$reconstructed)
names(rename_nmf_signatures(fit_res_merged_reduced, signatures, cutoff = 0.87)$signature)



p=plot_contribution_heatmap(fit_res_merged_std$contribution, 
                            cluster_samples = TRUE, 
                            cluster_sigs = TRUE)
pdf(paste0("./fig/", i, ".pdf"), height=6.6, width=5)
p
dev.off() 
i=i+1


p=plot_profile_heatmap(as.matrix.data.frame(fit_res_merged_std$signatures), by=merged_sig_name)
pdf(paste0("./fig/", i, ".pdf"), height=10, width=12)
p
dev.off() 
i=i+1

contri_boots <- fit_to_signatures_bootstrapped(mut_mat2, 
                                               as.matrix.data.frame(fit_res_merged_std$signature),
                                               n_boots = 150,
                                               method = "strict"
)

p=plot_bootstrapped_contribution(contri_boots, plot_type = "barplot")
pdf(paste0("./fig/", i, ".pdf"), height=35, width=3)
p
dev.off() 
i=i+1


p=plot_bootstrapped_contribution(contri_boots, 
                               mode = "relative", 
                               plot_type = "dotplot")
pdf(paste0("./fig/", i, ".pdf"), height=12, width=6)
p
dev.off() 
i=i+1



### end of file -- MutationalPatterns_.R