library(latticeExtra)

abundances <- read.csv(file = Sys.getenv("HORI_CSV"))
width <- strtoi(Sys.getenv("HORI_WIDTH"))
height <- strtoi(Sys.getenv("HORI_HEIGHT"))

png(Sys.getenv("HORI_IMAGE"), width = width, height = height)

horizonplot(ts(abundances), colorkey = TRUE, layout = c(1, ncol(abundances)),
            strip.left = FALSE, xlab = "Sample ID", ylab = list(rev(colnames(abundances)), rot = 0, cex = 1.35))

dev.off()
