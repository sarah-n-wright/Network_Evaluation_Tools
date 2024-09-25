# script.R

# get inputs
args <- commandArgs(trailingOnly=TRUE)
gofile <- args[1]
edgefile <- args[2]

min_genes <- as.numeric(args[3])
max_genes <- as.numeric(args[4])

stat <- args[5]
nfolds <- as.numeric(args[6])
outpath <- args[7]

library(EGAD)
library(data.table)

# read the go data
go_data <- fread(gofile)
go_data <- go_data[, .(GeneID, GO)]

# Get network data
net_df <- read.table(edgefile, sep='\t', header=TRUE)[, 1:2]
genelist <- make_genelist(net_df)
gene_network <- make_gene_network(net_df, genelist)

# subset to go terms with at least min-genes and at most max_genes in genelist
go_data_netset = go_data[GeneID %in% genelist]
go_counts <- go_data_netset[, .(num_genes = uniqueN(GeneID)), by = GO]
selected_go <- go_counts[num_genes >= min_genes & num_genes <= max_genes, GO]
go_data_final <- go_data_netset[GO %in% selected_go]

# create input data
goterms <- unique(go_data_final[,.(GO)][[1]])
annotations <- make_annotations(go_data_final[, .(GeneID, GO)], genelist, goterms)
# Set up inputs to neighbor voting (taken from run_GBA.R)
m <- match(rownames(gene_network), rownames(annotations))
f <- !is.na(m)
g <- m[f]

network.sub <- gene_network[f,f]
genes.labels <- filter_network_cols(annotations[g, ], min_genes, max_genes)

# Perform neighbor voting and save results
roc_sub <- neighbor_voting(genes.labels, network.sub, nFold=nfolds, output=stat)
roc_df <- as.data.frame(roc_sub)

# mutifunc
optimal <- calculate_multifunc(genes.labels) 
optimal_list <- as.numeric(optimal[,4]) 
roc_df$auc_multifunc <- auc_multifunc(genes.labels, optimal_list)

write.csv(roc_df, file=paste0(outpath,'.roc.txt'))

genes <- predictions(genes.labels, network.sub)
genes_df <- as.data.frame(genes)
write.csv(genes_df, file=paste0(outpath, '.pred.txt'))

