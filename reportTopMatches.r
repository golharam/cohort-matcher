args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cohort_matcher_results <- "cohort-matcher-results.txt"
  total_compared_file <- "total_compared.txt"
} else {
  cohort_matcher_results <- args[1]
  total_compared_file <- args[2]
}
paste("Reading", cohort_matcher_results, sep=" ")
table <- read.table(cohort_matcher_results, header=TRUE, row.names=1)
paste("Reading", total_compared_file, sep=" ")
total_compared <- read.table(total_compared_file, header=TRUE, row.names=1)

reportTopMatches <- function(x, ...) {
  f <- file("topmatches.txt", "w")
  writeLines(paste("sample", "match1", "score1", "match2", "score2", "match3", "score3", "match4", "score4", "match5", "score5", sep="\t"), f)
  for (sample in rownames(x)) {
    sample_matches <- sort(x[sample,], decreasing=TRUE)
    writeLines(paste(sample,
                     colnames(sample_matches[1]), sample_matches[1,1],
                     colnames(sample_matches[2]), sample_matches[1,2],
                     colnames(sample_matches[3]), sample_matches[1,3],
                     colnames(sample_matches[4]), sample_matches[1,4],
                     colnames(sample_matches[5]), sample_matches[1,5],
                     sep="\t"),
               f)
  }
  close(f)
}

# Plot
# http://www.phaget4.org/R/image_matrix.html
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...) {
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}

plotNumSNPsCompared <- function(x, ...) {
    heatmap.2(data.matrix(x), 
          # dendrogram control
          Rowv=NULL, Colv=NULL, 
          dendrogram='none',

          # data scaling
          scale="none",

          # mapping data to colors
          symbreaks=FALSE, 

          # colors
          col = rev(rainbow(20*10, start = 0/6, end = 4/6)), 

          # level trace
          trace='none', 

          # Row/Column Labeling
          margins=c(3,0), # ("margin.Y", "margin.X")
          labRow=NA,
          labCol=NA,
          
          # color key + density info
          keysize=1, 
          density.info='histogram', 
          denscol="black",
          symkey=FALSE, 
          key.title="# of SNPs",
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(3.5,0,3,0)),

           # plot labels
           main = NULL,
           xlab = NULL,
           ylab = NULL,
          
          # plot layout
          # lmat - visual layout: position matrix
          # lhei - column height
          # lwid - column width
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), 
          lhei=c(2.5, 5),
          lwid=c(1, 10, 1))

}

# ----- END plot function ----- #

paste("Writing top matches in topmatches.txt", sep=" ")
reportTopMatches(table)

M.table <- as.matrix(table)
pdfFile <- gsub(".txt", ".pdf", cohort_matcher_results)
paste("Writing", pdfFile, sep=" ")
pdf(pdfFile)
myImagePlot(M.table)
dev.off()

paste("Plotting total_compared.pdf")
library(gplots)
pdf("total_compared.pdf")
plotNumSNPsCompared(data.matrix(total_compared))
dev.off()
