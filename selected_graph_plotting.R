# Making graphs
library(igraph)
library(scales)

# Pasting single colonies

location <- "C:\\Users\\caaguila\\Dropbox\\Neurospora_crassa_experiments_JYU\\inoculation_07112024\\processed data\\"
folder <- "results\\"
#folder <- "Hf DD2 results//"
#folder <- "Pi results//"
#folder <- "Pv results//"


colonies_ntwk <- readRDS(paste(location, folder, "colonies_ntwk.RDS", sep = ""))
spatial.data <- readRDS(paste(location, folder, "spatial.data.RDS", sep = ""))

#Extracting the names of each colony

nombres <- sapply(colonies_ntwk, function(x){unique(E(x)$name)})

# Getting those names
nombres

# Selecting which colony to plot based on a specific name
r <- which(nombres == "neurospora_1_halfpda1_fc_Merged")
rb <- colonies_ntwk[[r]]

# Adding coordinates to the nodes in the plot

V(rb)$coordinate_x <- spatial.data[[r]][,1]
V(rb)$coordinate_y <- spatial.data[[r]][,2]

# Subsetting to remove the edges connected to the root. I do this because those edges are meaningless and distort the plotting
rb_all <- rb
rb <- subgraph(rb, V(rb)[Accessibility>0])

# Plot, playing with color coding

# The most meningful color is based by width
ecol <- E(rb)$width

# Playing with color palettes:
parula_colors <- c("#352A87", "#197EC0", "#22A784", "#FDE333", "#FD9A14")

ecol2 <- cscale(ecol, palette = gradient_n_pal(parula_colors))
ecol2[E(rb)[type == "Inoculum"]] <- "gray30"

# Directly exporting the plot as png
png("Figures//tests.png",
    width = 225, height = 225, units='mm', res = 300)

par(mar=c(0,0,0,0), bg = "white")
plot(rb,
     edge.arrow.size=1,edge.curved=0,
     edge.width=2,
     #edge.width=ecol/3,
     vertex.label=NA,vertex.shape="none",
     edge.color = ecol2,
     edge.size = 150,vertex.size=0,
     layout = spatial.data[[r]]*1); 
title(unique(E(rb)$name), line = -1, col.main ="black")

dev.off()
