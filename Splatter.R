source("http://bioconductor.org/biocLite.R")
biocLite("zinbwave")
library(devtools)
install_github("Oshlack/splatter")

library(splatter)
library(monocle)
biocLite(c("DDRTree", "pheatmap"))
browseVignettes("splatter")



# Load example data
#data("sc_example_counts")
sc_example_counts <- data.matrix(read.csv("/Users/vimal/Desktop/GSE46226_SC_Expression.csv", header = T, row.names = 1))

# Estimate parameters from example data
#params <- splatEstimate(sc_example_counts)
#test_data[test_data<0] <- 0
#params1 <- splatEstimate(test_data)
params = newSplatParams()
params1 <- setParams(params, update = list( nGenes = 1000, batchCells = 31, mean.rate = 0.1, mean.shape = 0.31))
params2 <- setParams(params, update = list( nGenes = 1000, batchCells = 31, mean.rate = 0.2, mean.shape = 0.33))
params3 <- setParams(params, update = list( nGenes = 1000, batchCells = 31, mean.rate = 0.3, mean.shape = 0.37))

# Simulate data using estimated parameters
sim1 <- splatSimulate(params, batchCells = 20, dropout.present = FALSE)
sim2 <- splatSimulate(params, batchCells = 30, dropout.present = FALSE)

sim.groups <- splatSimulate(params, group.prob = c(0.3, 0.3, 0.4), method = "groups",verbose = FALSE)
plotPCA(sim.groups, colour_by = "Group", exprs_values = "counts")

sim.path1 <- splatSimulatePaths(params1, batchCells = 500, verbose = FALSE)
sim.path2 <- splatSimulatePaths(params2, batchCells = 500, verbose = FALSE)
sim.path3 <- splatSimulatePaths(params3, batchCells = 500, verbose = FALSE)

p1 <- plotPCA(sim.path1, colour_by = "Step", exprs_values = "counts")
p2 <- plotPCA(sim.path2, colour_by = "Step", exprs_values = "counts")
p3 <- plotPCA(sim.path3, colour_by = "Step", exprs_values = "counts")

grid.arrange(p1, p2, p3, ncol = 3)





comparison <- compareSCEs(list(Splat = sim1, Simple = sim2))
names(comparison)
names(comparison$Plots)