requireNamespace(plotly, quietly = TRUE)
requireNamespace(reticulate, quietly = TRUE)

setwd('/code/rgeomstats')
source('R/visualization.R')

# use_python('/usr/local/bin/python3', required = T)
gs <- import("geomstats")

set.seed(1004)

dimension = 2L
h2 <- gs$hyperbolic_space$HyperbolicSpace(dimension = dimension)
data = h2$random_uniform(n_samples=15L, bound=1.)

disk <- PoincareDisk$new()

p <- plot_ly()
p <- disk$DrawPoints(p, data)
p
