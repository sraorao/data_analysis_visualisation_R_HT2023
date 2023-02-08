###R for biologists
##Irina & Rao, 22/02/2023
library(tidyverse)
library(RColorBrewer)

# LOAD and prepare data from files ####
# we first make a list of filenames of interest using the list.files() function
filenames = list.files(path = "Session3/data", pattern = "counts.txt", full.names = TRUE)

# we then use lapply(), which is similar to a for() loop
# to read the tab separated table of count data
# so that we get a list of data.frames in count_data_list
count_data_list = lapply(filenames, function(x) {
                    sample_name = substr(x, 15, 16)
                    each_sample_count = read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                                                   col.names = c("mir_name", sample_name))
                    return(each_sample_count)
                  })
head(count_data_list[[1]])
str(count_data_list)

# Optional exercise
# Use the sub() function and regular expressions to extract the sample names (e1, m1, etc.) from the filenames


# we then merge all the dataframes recursively using Reduce
norm_counts = Reduce(function(x, y) {merge(x, y, by = "mir_name")}, count_data_list)

View(norm_counts)

conditions = as.factor(substr(colnames(norm_counts)[2:7], 1, 1))

# load results table
results = read.csv("Session3/data/differentially_expressed_mirs_significant.csv", stringsAsFactors = FALSE)
# COLOURS ####
# Print available palettes
display.brewer.all()

# set colours
condition_colours = brewer.pal(2, "Accent")[conditions]
bluegreen_colours = colorRampPalette(brewer.pal(9, "GnBu"))(100)

# Draw HEATMAP ####
# install.packages("gplots)
library(gplots) # do not confuse with ggplot2

# sort results by adjusted p value and pick the top 50 microRNAs
results = arrange(results, padj)
selected_mirs = results$mir[1:50]

# filter norm_counts for only the selected 50 microRNAs
norm_counts %>%
  filter(mir_name %in% selected_mirs) -> norm_counts_selected

# shorten the miR name and assign the result to rowname
# then exclude the mir_name column (because it's in the row names now)
rownames(norm_counts_selected) <- sub("hsa-", "", norm_counts_selected$mir_name)
norm_counts_selected = norm_counts_selected[, 2:7]

# convert to matrix
norm_counts_selected = as.matrix(norm_counts_selected)

# plot heatmap
pdf("heatmap_top50.pdf", width = 8)
heatmap.2(norm_counts_selected, 
          scale = "row", 
          trace = "none", 
          col = bluegreen_colours, 
          Rowv = FALSE, 
          dendrogram = "column", 
          ColSideColors = condition_colours)
dev.off()

# why scaling is important
heatmap.2(norm_counts_selected, 
          trace = "none", 
          col = bluegreen_colours, 
          Rowv = FALSE, 
          dendrogram = "column", 
          ColSideColors = condition_colours)

# let's try plotting the counts for just one miR
norm_counts %>%
  filter(mir_name == "hsa-mir-567") %>%                       # filter for miR of interest
  gather(key = "sample", value = "count", -mir_name) %>%      # convert to long format
  mutate(condition = as.factor(substr(sample, 1, 1))) %>%     # add a column that defines the condition
  ggplot(aes(x = condition, y = count)) +                     # parent ggplot function defines x and y axes
    geom_boxplot(outlier.shape = NA) +                        # boxplot
    geom_dotplot(binaxis = "y", stackdir = "center") +        # dots overlaid on top of box plot
    theme_bw() -> p                                           # set theme to B & W
p

# Principal Component Analysis ####
# PCA plot with ggfortify
library(ggfortify)
norm_counts_transposed = data.frame(t(norm_counts[, 2:7]))
names(norm_counts_transposed) = norm_counts$mir_name
norm_counts_transposed$cell_type = sub("[123]", "", row.names(norm_counts_transposed))
pca_residues = prcomp(norm_counts_transposed[, -ncol(norm_counts_transposed)], scale. = TRUE)
plot(pca_residues$sdev^2/sum(pca_residues$sdev^2), type = "b") # scree plot
autoplot(pca_residues, data = norm_counts_transposed, colour = "cell_type") 

# Exercise: Plot a barplot for hsa-miR-372-3p (one bar for each sample; sample name on the x axis and counts on the y axis)
# hint: scale_fill_brewer(palette = "Set3") # set colours using RColorBrewer palettes
# hint: scale_fill_manual(values = unique(condition_colours)) # set colours manually


# interactive graphs with the plotly package ####
# install.packages("plotly")
library(plotly)
# bubble plot from owid_covid_newyear
# Make the plot with ggplot2 first...
owid_covid = read.csv("Session1/data/owid-covid-data.csv", header = TRUE,
                      stringsAsFactors = FALSE)
owid_covid$date = as.Date(owid_covid$date)
new_year = "2021-01-01"
owid_covid %>%
  filter(date == new_year) -> owid_covid_newyear
library(ggrepel)
(owid_covid_newyear %>%
  ggplot(aes(x = total_cases, y = total_deaths, label = location)) +
    geom_point(aes(fill = continent, size = population), shape = 21) +
    scale_size(range = c(2, 40)) +
    scale_x_log10() +
    scale_y_log10() -> p)
p
# ...and then use the ggplot2 plot as input for ggplotly()
ggplotly(p)

# Volcano plot ----
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

EnhancedVolcano(results, lab = results$mir, 
                x = 'log2FoldChange', y = 'padj')

# phylogenetic trees ----
library(ggtree)
library(treeio)

mytree = "(B, (C, D))A;"
mytree = read.tree(file = textConnection(mytree))

mytree 

ggtree(mytree) +
  geom_tippoint() +
  geom_nodepoint()


