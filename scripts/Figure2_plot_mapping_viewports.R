#If you get a "no package called Gvix" error message, run the following code:
# source("https://bioconductor.org/biocLite.R")
# biocLite("Gviz")

library(Gviz)

library(GenomicRanges)

options(ucscChromosomeNames=FALSE)

gtrack <- GenomeAxisTrack(chromosome="X")

dup.data <- read.csv("../data/plot_dup_data.csv", stringsAsFactors = TRUE)
dup.data$stock <- as.factor(dup.data$stock)
dup.data$fill <- ifelse(dup.data$rescues=="y", "#5555AA", "#AA5555")

dup.track <- AnnotationTrack(	chromosome=rep("X", nrow(dup.data)),
								name="Duplications",
								start=dup.data$start,
								end=dup.data$end,
								group=dup.data$stock,
								id=dup.data$stock,
								groupAnnotation="id",
								cex.group=1.0,
								fill=dup.data$fill,
								rotation=90,
								just.group="below"
								)

rescued.data <- dup.data[dup.data$rescues=="y",]
nonrescued.data <- dup.data[dup.data$rescues=="n",]


rescued.track <- AnnotationTrack(	chromosome=rep("X", nrow(rescued.data)),
								name="Rescued duplications",
								start=rescued.data$start,
								end=rescued.data$end,
								group=rescued.data$stock,
								id=rescued.data$stock,
								groupAnnotation="id",
								cex.group=1.0,
								fill=rescued.data$fill,
								rotation=90,
								just.group="below"
								)


nonrescued.track <- AnnotationTrack(	chromosome=rep("X", nrow(nonrescued.data)),
								name="Non-rescued duplications",
								start=nonrescued.data$start,
								end=nonrescued.data$end,
								group=nonrescued.data$stock,
								id=nonrescued.data$stock,
								groupAnnotation="id",
								cex.group=1.0,
								fill=nonrescued.data$fill,
								rotation=90,
								just.group="below"
								)




gene.data <- read.csv("../data/gene_data.csv", stringsAsFactors=F)
gene.track <- AnnotationTrack(chromosome=gene.data$chr, start=gene.data$start, end=gene.data$end, strand=gene.data$strand, id=gene.data$id, groupAnnotation="id", just.group="below", name="Gene models")
gene.track

pdf(file="../outputs/Figure2_Mapping.pdf", width=8, height=8)

grid.newpage()

text.pars <- gpar(cex=1.5, fontface="bold")

pushViewport(viewport(height=0.475, y=1, just="top"))
plotTracks(list(gtrack, rescued.track, nonrescued.track), from=min(dup.data$start), to=max(dup.data$end), extend.left=0.1, extend.right=0.1, add=T)
grid.text(label="A", x=0.01, y=0.98, just=c("left", "top"), gp=text.pars)
popViewport(1)

pushViewport(viewport(height=0.475, y=0, just="bottom"))
ht <- HighlightTrack(trackList=list(gtrack, rescued.track, nonrescued.track, gene.track), start=9686653, end=9762229, chromosome="X", inBackground=F, fill="#00000011")
plotTracks(list( ht), from=9677341, to=9784700, extend.left=0.1, extend.right=0.1, add=T)
grid.text(label="B", x=0.01, y=0.98, just=c("left", "top"), gp=text.pars)
popViewport(1)

dev.off()

sessionInfo()
