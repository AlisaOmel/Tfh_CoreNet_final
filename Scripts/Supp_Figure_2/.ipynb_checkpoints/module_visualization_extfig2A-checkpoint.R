library(qgraph)

#################################################################
##               Making only sig modules network               ##
#################################################################
rm(list = ls())
library(viridis)
library(readr)


#PPISig <- read.csv("/ix/djishnu/Alisa/Tfh/ForPaper/first_degree_int/pps_protein_pairs_sig_genes_df.csv")
#ref <- read.table("/ix/djishnu/Javad/RA/Data/HomoSapiens_binary_co_complex_union.txt")

PPISig <- read.csv("../../Data/CoreNet_PPI_protein_pairs_sig_genes_df.csv")
ref <- read.table("../../Data/HomoSapiens_binary_co_complex_union.txt")

ref <- ref[!((ref$V1%in%PPISig$Prot_1)&(ref$V1%in%PPISig$Prot_2)),]
set.seed(42)  # For reproducibility

#Subset reference since full network is too big to visualize
refS <- ref[sample(size=2000,dim(ref)[1]),]
colnames(refS) <- colnames(PPISig)
allPPI <- rbind(PPISig, refS)
colnames(allPPI) <- colnames(PPISig)

####Merged output oject was created by combining the hotnet2 modules that passed the threshold across the multiple deltas that were tested. 
####Duplicates were removed and it was saved into an RDS file.

#merged_data <- strsplit(read_file("../../Data/merged_output_PPS_1_edited.txt"), "\n")[[1]]
#merged_list <- lapply(merged_data, function(x) strsplit(x, "\t")[[1]])
#saveRDS(merged_list, "../../Data/merged_output.RDS")

#Modules <- readRDS("/ix/djishnu/Alisa/Tfh/ForPaper/process_network_prop_output/merged_output.RDS")
Modules <- readRDS("../../Data/merged_output.rds")
#Make a module for reference subset
Modules[[28]] <- unique(c(refS$Prot_1,refS$Prot_2))

#################################################################
##                     Decode it to numbers                     ##
#################################################################
# The qgplot doesn't perfrom well with groups when the node are in character
unique_genes <- unique(unlist(Modules))
gene_mapping <- data.frame(gene = unique_genes, number = seq_along(unique_genes))

## Making groups
grps <- lapply(Modules, function(module) {
  # Mapping genes to numbers using gene_mapping
  numbers <- gene_mapping$number[match(module, gene_mapping$gene)]
  return(numbers)
})

allPPI <- cbind(gene_mapping$number[match(allPPI$Prot_1,gene_mapping$gene)],
                gene_mapping$number[match(allPPI$Prot_2,gene_mapping$gene)])


#################################################################
##                     Plot it                                 ##
#################################################################

report <- data.frame(k=numeric(),number_nodes = numeric(),
                     number_edges = numeric())

for(i in c(seq_along(Modules))){
  
  number_nodes  <- length(Modules[[i]])
  number_edges <- nrow(allPPI[(allPPI[,1] %in% grps[[i]]) & (allPPI[,2] %in% grps[[i]]) ,])
  k <- i
  report <- rbind(report,data.frame(k=i,number_nodes=number_nodes,number_edges=number_edges))
}



# Set layout: main plot on top, legend below
layout(matrix(c(1, 2), ncol = 1), heights = c(4, 1.5))

my_palette <- c(hcl(seq(15, 375, length.out = 27), c = 100, l = 65), "#FFFFFF")

## --- Plot 1: qgraph ---
par(mar = c(1, 1, 1, 1))  # minimal margins
qgraph(allPPI,
       directed = FALSE,
       layout = "spring",
       groups = grps,
       vsize = 1,
       labels = FALSE,
       repulsion = 1,
       color = my_palette
)

## --- Plot 2: Legend ---
par(mar = c(0, 0, 0, 0))
plot.new()

num_modules <- length(my_palette)
cols <- ceiling(num_modules / 3)

legend("center",
       legend = 1:num_modules,    # or use names(grps) if available
       fill = my_palette,
       title = "Modules",
       ncol = cols,
       cex = 1.5,                  # label font size (20pt equivalent)
       title.cex = 1.5,            # title font size
       bty = "n")                # no box

png("Figure2_ppi_network_with_legend.png",
     width = 2400, height = 3000, res = 200)

# Set layout: main plot on top, legend below
layout(matrix(c(1, 2), ncol = 1), heights = c(4, 1.5))

my_palette <- c(hcl(seq(15, 375, length.out = 27), c = 100, l = 65), "#FFFFFF")

## --- Plot 1: qgraph ---
par(mar = c(1, 1, 1, 1))
qgraph(allPPI,
       directed = FALSE,
       layout = "spring",
       groups = grps,
       vsize = 1,
       labels = FALSE,
       repulsion = 1,
       color = my_palette)

## --- Plot 2: Legend ---
par(mar = c(0, 0, 0, 0))
plot.new()
num_modules <- length(my_palette)
cols <- ceiling(num_modules / 3)
legend("center",
       legend = 1:num_modules,
       fill = my_palette,
       title = "Modules",
       ncol = cols,
       cex = 1.5,
       title.cex = 1.5,
       bty = "n")

dev.off()