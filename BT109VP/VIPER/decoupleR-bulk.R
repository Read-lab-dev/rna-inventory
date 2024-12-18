## We load the required packages
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
load("~/rna/Bulk_Analysis/DIPG/rawcounts_coding.Rdata")
sample <- clipr::read_clip()
metadata <- metadata[sample,]
exprSet <- exprSet[,sample]
# Remove NAs and set row names
counts <- exprSet %>% 
  dplyr::mutate_if(~ any(is.na(.x)), ~ if_else(is.na(.x),0,.x)) %>% 
  as.matrix()

design <- data.frame(sample=rownames(metadata),condition=metadata$Treat)
# Extract t-values per gene
load("/home/hzg/rna/Bulk_Analysis/DIPG/BT109_VP_Long_GSEA.Rdata")
deg <- resdata %>%
  dplyr::select(gene, log2FoldChange, stat, pvalue) %>%
  filter(!is.na(t)) %>%
  column_to_rownames(var = "gene") %>%
  as.matrix()
head(deg)
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))
#Run wmean
sample_acts <- run_wmean(mat=counts, net=regulon, .source='tf', .target='target',
                         .mor='mor', times = 100, minsize = 5)
sample_acts

###Plot####
n_tfs <- 25

# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# Get top tfs with more variable means across clusters
tfs <- sample_acts %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)
sample_acts_mat <- sample_acts_mat[,tfs]

# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)
rownames(sample_acts_mat) <- design$condition

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap(sample_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 

# Run wmean
contrast_acts <- run_wmean(mat=deg[, 'stat', drop=FALSE], net=regulon, .source='tf', .target='target',
                           .mor='mor', times = 100, minsize = 5)
contrast_acts
# Filter norm_wmean
f_contrast_acts <- contrast_acts %>%
  filter(statistic == 'norm_wmean') %>%
  mutate(rnk = NA)

# Filter top TFs in both signs
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
tfs <- f_contrast_acts %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% c(tfs,"TEAD1",'TEAD4'))

# Plot
ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_viridis(option = "D")+ 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")


#########SPECIFIC########
tf_select <- 'MYCN'

df <- regulon %>%
  filter(tf == tf_select) %>%
  arrange(target) %>%
  mutate(ID = target, color = "3") %>%
  column_to_rownames('target')

inter <- sort(intersect(rownames(deg),rownames(df)))
df <- df[inter, ]
df[,c('logfc', 't_value', 'p_value')] <- deg[inter, ]
df <- df %>%
  mutate(color = if_else(mor > 0 & t_value > 0, '1', color)) %>%
  mutate(color = if_else(mor > 0 & t_value < 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & t_value > 0, '2', color)) %>%
  mutate(color = if_else(mor < 0 & t_value < 0, '1', color))

ggplot(df, aes(x = logfc, y = -log10(p_value), color = color, size=abs(mor))) +
  geom_point() +
  scale_colour_manual(values = c("red","royalblue3","grey")) +
  geom_label_repel(aes(label = ID, size=1)) + 
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  ggtitle(tf_select)

