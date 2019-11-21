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
cnet <- consensusNet(gene.trees, 0.15)
# find edges that are in the network but not in the tree
edge.col <- createLabel(cnet, species.tree, label = "black", "edge", nomatch = "grey")
# set plot margins
par(mar = c(0, 0, 0, 0))
# plot consensus network and map bipartition frequencies from a sample of trees onto corresponding network edges
cnet.freq <- plot(cnet,
                  type = "3D",
                  #use.edge.length = FALSE,
                  show.edge.label = TRUE,
                  edge.label = round(as.numeric(cnet$edge.labels) * 0.01, 2),
                  edge.color = edge.col,
                  font = 1,
                  cex = 0.7,
                  font.edge.label = 2,
                  cex.edge.label = 0.9,
                  col.edge.label = "blue")
# write nexus file for SplitsTree
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
