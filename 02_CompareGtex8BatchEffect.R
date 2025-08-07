# 02: Plot out GTEx v8 batch effects for supp.
# Plot out the batch effects of GTEx v8 before and after ComBAT using fast rpca.

library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(DESeq2)
library(Rtsne)
library(rsvd)
library(ggpubr)

gtex_dir <- "../data/gtex"
in_dir <- "01_out_cached"
out_dir <- "02_out"
if (!dir.exists(out_dur)) {
  dir.create(out_dur)
}

# Read the read count files and make the uncorrected GTEx matrix
h <- readr::read_tsv(file.path(gtex_dir, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"), skip=2)
gtex_matrix <- h %>% dplyr::select(-Name) |>
  dplyr::distinct(Description, .keep_all=TRUE) |>
  tibble::column_to_rownames(var="Description") |>
  as.matrix()

# Read the annotation file
annot <- readr::read_tsv(file.path(gtex_dir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"), guess_max=15000)
annot$SUBJID <- gsub("(?<!GTEX)-.*$", "", annot$SAMPID, perl=T)

# Read the subject demographics file
subdemo <- readr::read_tsv(file.path(gtex_dir, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"))

# Get the metadata file for all RNA-seq data
i. <- annot %>% dplyr::filter(SMAFRZE == "RNASEQ")
i. <- i. %>% dplyr::left_join(subdemo)
i.$rowname <- i.$SAMPID
i. <- i. %>% column_to_rownames("rowname")
i. <- i.[colnames(gtex_matrix),]

gtex <- DESeqDataSetFromMatrix(countData = gtex_matrix,
                              colData = i.,
                              design= ~ SMTSD)

# Read the saved VST and ComBAT objects
vsd <- readRDS(file.path(in_dir, "gtexv8_vst_dds.Rds"))
vsd_tec <- readRDS(file.path(in_dir, "gtexv8_vst_combat1_dds.Rds"))
vsd_tec_tis <- readRDS(file.path(in_dir, "gtexv8_vst_combat2_dds.Rds"))
vsd_tec_bio <- readRDS(file.path(in_dir, "gtexv8_vst_combat3_dds.Rds"))
vsd_tec_bio_age <- readRDS(file.path(in_dir, "gtexv8_vst_combat4_dds.Rds"))

# We should do a PCA plot before and after each step of normalization/ batch correction
top_genes <- order(rowVars(assay(gtex)), decreasing=TRUE)[1:2000]
p0 <- rsvd::rpca(assay(gtex)[top_genes,], k=2)
p0_df <- p0$rotation %>% as.data.frame() %>% tibble::rownames_to_column() %>% tibble::as_tibble() %>% dplyr::select(SAMPID = rowname, PC1:PC2) %>% dplyr::left_join(i.)
g0 <- ggplot(data=p0_df, aes(x=PC1, y=PC2)) 
g0 <- g0 + geom_point(aes(col=SMTS), size=0.7, stroke=0) + theme_classic() + theme(aspect.ratio=1, legend.position="bottom")
g0 <- g0 + ggtitle("Raw Counts")
ggsave(plot=g0, file= file.path(out_dir, "GTExv8_batchCorrection_0.pdf"), width=5, height=5, useDingbats=F)

# Use t-SNE instead for another look
p0_tsne <- Rtsne::Rtsne(t(assay(gtex)[top_genes,]), perplexity=30, verbose=TRUE, check_duplicates=FALSE)
p0_tsne_df <- data.frame(SAMPID = i.$SAMPID, TSNE1 = p0_tsne$Y[,1], TSNE2 = p0_tsne$Y[,2]) %>% dplyr::left_join(i.)
g0_tsne <- ggplot(data=p0_tsne_df, aes(x=TSNE1, y=TSNE2))
g0_tsne <- g0_tsne + geom_point(aes(col=SMTS), size=0.7, stroke=0) + theme_classic() + theme(aspect.ratio=1, legend.position="bottom")
g0_tsne <- g0_tsne + ggtitle("Raw Counts t-SNE")
ggsave(plot=g0_tsne, file= file.path(out_dir, "GTExv8_batchCorrection_0_tsne.pdf"), width=5, height=5, useDingbats=F)

# After VST
top_genes_vsd <- order(rowVars(assay(vsd)), decreasing=TRUE)[1:2000]
p1 <- rsvd::rpca(assay(vsd), k=2)
p1_df <- p1$rotation %>% as.data.frame() %>% tibble::rownames_to_column() %>% tibble::as_tibble() %>% dplyr::select(SAMPID = rowname, PC1:PC2) %>% dplyr::left_join(i.)
g1 <- ggplot(data=p1_df, aes(x=PC1, y=PC2))
g1 <- g1 + geom_point(aes(col=SMTS), size=0.7, stroke=0) + theme_classic() + theme(aspect.ratio=1, legend.position="bottom")
g1 <- g1 + ggtitle("After Variance Stabilization")
ggsave(plot=g1, file= file.path(out_dir,"GTExv8_batchCorrection_1.pdf"), width=5, height=5, useDingbats=F)

# Use t-SNE instead
p1_tsne <- Rtsne::Rtsne(t(assay(vsd)[top_genes_vsd,]), perplexity=30, verbose=TRUE, check_duplicates=FALSE)
p1_tsne_df <- data.frame(SAMPID = i.$SAMPID, TSNE1 = p1_tsne$Y[,1], TSNE2 = p1_tsne$Y[,2]) %>% dplyr::left_join(i.)
g1_tsne <- ggplot(data=p1_tsne_df, aes(x=TSNE1, y=TSNE2))
g1_tsne <- g1_tsne + geom_point(aes(col=SMTS), size=0.7, stroke=0) + theme_classic() + theme(aspect.ratio=1, legend.position="bottom")
g1_tsne <- g1_tsne + ggtitle("After Variance Stabilization t-SNE")
ggsave(plot=g1_tsne, file= file.path(out_dir, "GTExv8_batchCorrection_1_tsne.pdf"), width=5, height=5, useDingbats=F)

# After first batch correction
top_genes_tec <- order(rowVars(assay(vsd_tec)), decreasing=TRUE)[1:2000]
p2 <- rsvd::rpca(assay(vsd_tec)[top_genes_tec, ], k=2)
p2_df <- p2$rotation %>% as.data.frame() %>% tibble::rownames_to_column()  %>% tibble::as_tibble() %>% dplyr::select(SAMPID = rowname, PC1:PC2) %>% dplyr::left_join(i.)
g2 <- ggplot(data=p2_df, aes(x=PC1, y=PC2))
g2 <- g2 + geom_point(aes(col=SMTS), size=0.7, stroke=0) + theme_classic() + theme(aspect.ratio=1, legend.position="bottom")
g2 <- g2 + ggtitle("After Stepwise Technical Batch Correction")
ggsave(plot=g2, file= file.path(out_dir,"GTExv8_batchCorrection_2.pdf"), width=5, height=5, useDingbats=F)

# Use t-SNE instead
p2_tsne <- Rtsne::Rtsne(t(assay(vsd_tec)[top_genes_tec,]), perplexity=30, verbose=TRUE, check_duplicates=FALSE)
p2_tsne_df <- data.frame(SAMPID = i.$SAMPID, TSNE1 = p2_tsne$Y[,1], TSNE2 = p2_tsne$Y[,2]) %>% dplyr::left_join(i.)
g2_tsne <- ggplot(data=p2_tsne_df, aes(x=TSNE1, y=TSNE2))
g2_tsne <- g2_tsne + geom_point(aes(col=SMTS), size=0.7, stroke=0) + theme_classic() + theme(aspect.ratio=1, legend.position="bottom")
g2_tsne <- g2_tsne + ggtitle("After Stepwise Technical Batch Correction t-SNE")
ggsave(plot=g2_tsne, file= file.path(out_dir, "GTExv8_batchCorrection_2_tsne.pdf"), width=5, height=5, useDingbats=F)


# After second batch correction
top_genes_tec_tis <- order(rowVars(assay(vsd_tec_tis)), decreasing=TRUE)[1:2000]
p3 <- rsvd::rpca(assay(vsd_tec_tis)[top_genes_tec_tis, ], k=2)
p3_df <- p3$rotation %>% as.data.frame() %>% tibble::rownames_to_column()  %>% tibble::as_tibble() %>% dplyr::select(SAMPID = rowname, PC1:PC2) %>% dplyr::left_join(i.)
g3 <- ggplot(data=p3_df, aes(x=PC1, y=PC2))
g3 <- g3 + geom_point(aes(col=SMTS), size=0.7, stroke=0) + theme_classic() + theme(aspect.ratio=1, legend.position="bottom")
g3 <- g3 + ggtitle("After Stepwise Tissue Charateristics Correction")
ggsave(plot=g3, file= file.path(out_dir,"GTExv8_batchCorrection_3.pdf"), width=5, height=5, useDingbats=F)

# Use t-SNE instead
p3_tsne <- Rtsne::Rtsne(t(assay(vsd_tec_tis)[top_genes_tec_tis,]), perplexity=30, verbose=TRUE, check_duplicates=FALSE)
p3_tsne_df <- data.frame(SAMPID = i.$SAMPID, TSNE1 = p3_tsne$Y[,1], TSNE2 = p3_tsne$Y[,2]) %>% dplyr::left_join(i.)
g3_tsne <- ggplot(data=p3_tsne_df, aes(x=TSNE1, y=TSNE2))
g3_tsne <- g3_tsne + geom_point(aes(col=SMTS), size=0.7, stroke=0) + theme_classic() + theme(aspect.ratio=1, legend.position="bottom")
g3_tsne <- g3_tsne + ggtitle("After Stepwise Tissue Charateristics Correction t-SNE")
ggsave(plot=g3_tsne, file= file.path(out_dir, "GTExv8_batchCorrection_3_tsne.pdf"), width=5, height=5, useDingbats=F)

# After third batch correction
top_genes_tec_bio <- order(rowVars(assay(vsd_tec_bio)), decreasing=TRUE)[1:2000]
p4 <- rsvd::rpca(assay(vsd_tec_bio)[top_genes_tec_bio, ], k=2)
p4_df <- p4$rotation %>% as.data.frame() %>% tibble::rownames_to_column()  %>% tibble::as_tibble() %>% dplyr::select(SAMPID = rowname, PC1:PC2) %>% dplyr::left_join(i.)
g4 <- ggplot(data=p4_df, aes(x=PC1, y=PC2))
g4 <- g4 + geom_point(aes(col=SMTS), size=0.7, stroke=0) + theme_classic() + theme(aspect.ratio=1, legend.position="bottom")
g4 <- g4 + ggtitle("After Stepwise Donor Death Charateristics Correction")
ggsave(plot=g4, file= file.path(out_dir,"GTExv8_batchCorrection_4.pdf"), width=5, height=5, useDingbats=F)

# Use t-SNE instead
p4_tsne <- Rtsne::Rtsne(t(assay(vsd_tec_bio)[top_genes_tec_bio,]), perplexity=30, verbose=TRUE, check_duplicates=FALSE)
p4_tsne_df <- data.frame(SAMPID = i.$SAMPID, TSNE1 = p4_tsne$Y[,1], TSNE2 = p4_tsne$Y[,2]) %>% dplyr::left_join(i.)
g4_tsne <- ggplot(data=p4_tsne_df, aes(x=TSNE1, y=TSNE2))
g4_tsne <- g4_tsne + geom_point(aes(col=SMTS), size=0.7, stroke=0) + theme_classic() + theme(aspect.ratio=1, legend.position="bottom")
g4_tsne <- g4_tsne + ggtitle("After Stepwise Donor Death Charateristics Correction t-SNE")
ggsave(plot=g4_tsne, file= file.path(out_dir, "GTExv8_batchCorrection_4_tsne.pdf"), width=5, height=5, useDingbats=F)

# After fourth batch correction
top_genes_tec_bio_age <- order(rowVars(assay(vsd_tec_bio_age)), decreasing=TRUE)[1:2000]
p5 <- rsvd::rpca(assay(vsd_tec_bio_age)[top_genes_tec_bio_age, ], k=2)
p5_df <- p5$rotation %>% as.data.frame() %>% tibble::rownames_to_column()  %>% tibble::as_tibble() %>% dplyr::select(SAMPID = rowname, PC1:PC2) %>% dplyr::left_join(i.)
g5 <- ggplot(data=p5_df, aes(x=PC1, y=PC2))
g5 <- g5 + geom_point(aes(col=SMTS), size=0.7, stroke=0) + theme_classic() + theme(aspect.ratio=1, legend.position="bottom")
g5 <- g5 + ggtitle("After Stepwise Age Correction")
ggsave(plot=g5, file= file.path(out_dir,"GTExv8_batchCorrection_5.pdf"), width=5, height=5, useDingbats=F)

# Use t-SNE instead
p5_tsne <- Rtsne::Rtsne(t(assay(vsd_tec_bio_age)[top_genes_tec_bio_age,]), perplexity=30, verbose=TRUE, check_duplicates=FALSE)
p5_tsne_df <- data.frame(SAMPID = i.$SAMPID, TSNE1 = p5_tsne$Y[,1], TSNE2 = p5_tsne$Y[,2]) %>% dplyr::left_join(i.)
g5_tsne <- ggplot(data=p5_tsne_df, aes(x=TSNE1, y=TSNE2))
g5_tsne <- g5_tsne + geom_point(aes(col=SMTS), size=0.7, stroke=0) + theme_classic() + theme(aspect.ratio=1, legend.position="bottom")
g5_tsne <- g5_tsne + ggtitle("After Stepwise Age Correction t-SNE")
ggsave(plot=g5_tsne, file= file.path(out_dir, "GTExv8_batchCorrection_5_tsne.pdf"), width=5, height=5, useDingbats=F)


