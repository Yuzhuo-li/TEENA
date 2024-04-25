### This R script can be used to take multiple TEENA result table to plot heatmaps
### Multiple parameters need to be manually adjusted for optimized visualization

### load the required R packages: make sure they are already install
library(pheatmap)
library(xlsx)


################################################################################
##### parameters
## below are demo with five files
## Note: these parameters need to be adjusted manually for optimized visualization

## TEENA result file names
file_names <- c(
  "GATA2_TEENA_result.xlsx",
  "GATA3_TEENA_result.xlsx",
  "MSX2_TEENA_result.xlsx",
  "TFAP2A_TEENA_result.xlsx",
  "TFAP2C_TEENA_result.xlsx"
  )

## sample names: should match the file names
sample_names <- c("GATA2", "GATA3", "MSX2", "TFAP2A", "TFAP2C")

## cut-off for significantly enriched TEs
pvalue_thres <- 0.01
enrich_thres <- 2

## output figure settings
# figure name
figure_name   <- "TE_heatmap.svg"
# width
figure_width  <- 15
# height
figure_height <- 18


################################################################################
##### read and process data

## get TE information
te_data <- read.xlsx2(file = file_names[1], sheetIndex = 1)
te_data <- te_data[order(te_data$TE_name),]
te_inf <- data.frame(
  TE_name = te_data$TE_name,
  TE_class = te_data$TE_class
)
te_inf <- te_inf[te_inf$TE_class %in% c("DNA", "LINE", "SINE", "LTR"),]
te_inf <- te_inf[!duplicated(te_inf$TE_name),]

## make TE matrix
pvalue_matrix <- matrix(0, nrow = length(te_inf$TE_name), ncol = length(sample_names))
row.names(pvalue_matrix) <- te_inf$TE_name
colnames(pvalue_matrix) <- sample_names
enrich_matrix <- pvalue_matrix

for(i in 1:length(file_names)){
  file_name <- file_names[i]
  d <- read.xlsx2(file = file_name, sheetIndex = 1)
  d <- d[order(d$TE_name),]
  d <- d[d$TE_name %in% te_inf$TE_name,]
  d <- d[!duplicated(d$TE_name),]
  pvalue_matrix[,i] <- as.numeric(d$p_adj)
  enrich_matrix[,i] <- as.numeric(d$fold_enrich)
}

## only keep significantly enriched TEs
pvalue_matrix_sig <- pvalue_matrix[apply(pvalue_matrix, 1, min) < pvalue_thres & apply(enrich_matrix, 1, max) > enrich_thres, ]
enrich_matrix_sig <- pvalue_matrix[apply(pvalue_matrix, 1, min) < pvalue_thres & apply(enrich_matrix, 1, max) > enrich_thres, ]


################################################################################
##### generate heatmap: the color gradients indicate the -log10(P)

## make TE class data frame
annotation_row <- data.frame(
  row.names = te_inf$TE_name,
  TE_class = te_inf$TE_class
)

## plot the heatmap of svg format
svg(filename = figure_name, width = figure_width, height = figure_height)
pheatmap(
  -log10(pvalue_matrix_sig),
  annotation_row = annotation_row,
  color = rev(heat.colors(50))
  )
dev.off()

