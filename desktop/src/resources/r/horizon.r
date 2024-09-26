abundances <- read.csv(file = Sys.getenv("HORI_CSV"))
width <- strtoi(Sys.getenv("HORI_WIDTH"))
height <- strtoi(Sys.getenv("HORI_HEIGHT"))
bandwidth <- strtoi(Sys.getenv("HORI_BANDWIDTH"))

library("latticeExtra")
png(Sys.getenv("HORI_IMAGE"), width = width, height = height)
horizonplot(ts(abundances), colorkey = TRUE, layout = c(1, ncol(abundances)),
            strip.left = FALSE, xlab = "Sample Date", ylab = list(rev(colnames(abundances)), rot = 0, cex = 1.35),
            horizonscale = bandwidth, origin = 0, col.regions = hcl.colors(14, palette = "RdYlBu"),)

dev.off()
