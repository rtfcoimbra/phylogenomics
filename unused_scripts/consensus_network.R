# sources: https://cran.r-project.org/web/packages/phangorn/vignettes/IntertwiningTreesAndNetworks.html
#          https://cran.r-project.org/web/packages/phangorn/vignettes/Networx.html

# load libraries
library(ggtree)
library(phangorn)

# set working directory
setwd("~/Downloads")
# import gene trees
gene.trees <- read.tree("./estimated_gene_trees.tree")
# import species tree
species.tree <- read.tree("./estimated_species_tree.tree")

# create consensus network
cnet <- consensusNet(gene.trees, 0.1)
# show frequencies of bipartitions found in the gene trees mapped on the corresponding edge bundles of the network
freq.bipart <- addConfidences(cnet, as.splits(gene.trees))
# find edges that are in the network but not in the tree
edge.col <- createLabel(freq.bipart, species.tree, label = "black", "edge", nomatch = "#F8766D")
# set plot parameters
par(mar = c(1, 1, 1, 1))
# plot consensus network with edge labels
cnet.freq <- plot(freq.bipart, "2D", show.edge.label = TRUE, edge.color = edge.col, col.edge.label = "#00BFC4", cex = 0.7)
# write to SplitsTree for viewing
write.nexus.networx(cnet.freq,"consensus_network.freq.nex")

# identify edge bundles (in black) in the consensus network that correspond to branches in the species tree
# create a vector of labels for the network corresponding to edges in the tree
#edge.lab <- createLabel(cnet, species.tree, species.tree$edge[,2], "edge")
# set plot parameters
#par(mfrow = c(1, 2), mar = c(1, 1, 1, 1))
# plot species tree with edge labels
#plot(species.tree, "unrooted", rotate.tree = 180, cex = 0.7)
#edgelabels(species.tree$edge[,2], col = "blue", frame = "none", cex = 0.7)
# find edges that are in the network but not in the tree
#edge.col <- createLabel(cnet, species.tree, "black", nomatch = "red")
# plot consensus network with edge labels
#plot(cnet, edge.label = edge.lab, show.edge.label = T, "2D", edge.color = edge.col, col.edge.label = "blue", cex = 0.7)

# show support values for all branches in the species tree mapped on the corresponding edge bundles of the network
#sup.values <- addConfidences(cnet, species.tree)
# find splits that are in the network but not in the tree
#split.col <- createLabel(sup.values, species.tree, label = "black", "split", nomatch = "grey")
# plot consensus network with edge labels
#cnet.sup <- plot(sup.values, "2D", show.edge.label = TRUE, split.color = split.col, col.edge.label = "blue")
# write to SplitsTree for viewing
#write.nexus.networx(cnet.sup,"consensus_network.sup.nex")
