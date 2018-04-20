library(gplots)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cohort_matcher_results <- "cohort-matcher-results.txt"
  total_compared_file <- "cohort-matcher-results.total_compared.txt"
} else {
  cohort_matcher_results <- args[1]
  total_compared_file <- args[2]
}
paste("Reading", cohort_matcher_results, sep=" ")
table <- read.table(cohort_matcher_results, header=TRUE, row.names=1)
paste("Reading", total_compared_file, sep=" ")
total_compared <- read.table(total_compared_file, header=TRUE, row.names=1)

reportTopMatches <- function(x, topMatchesFile, ...) {
  f <- file(topMatchesFile, "w")
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
plotSampleSimilarity <- function(x, ...) {
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
  #axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  #axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
  #     cex.axis=0.7)
  
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
  new <- TRUE
  if (new) {
    heatmap.2(data.matrix(x),
              # dendrogram control
              Rowv=NULL, Colv=NULL,
              dendrogram='none',

              # colors
              col = rev(rainbow(20*10, start = 0/6, end = 4/6)),

              # level trace
              trace='none',

              # Row/Column Labeling
              cexRow=0.5,
              cexCol=0.5,
              labRow=NA,
              labCol=NA,

              # color key + density info
              density.info='histogram',
              denscol="black",
              key.title="Frequency Distribution of Number of SNPs Compared between Samples",
              key.xlab = "# of SNPs",
              key.ylab = "# of Samples Compared",
              #( "bottom", "left", "top", "right" )
              key.par=list(mar=c(5.1, 4.1, 4.1, 2.1)),

              # plot labels
              main="Sample Similarity",

              # plot layout
              # lmat - visual layout: position matrix
              # lhei - column height
              # lwid - column width
              # lmat -- added 2 lattice sections (5 and 6) for padding
              # lmat is a matrix describing how the screen is to be broken up. By default,
              # heatmap.2 divides the screen into a four element grid, so lmat is a 2x2 matrix.
              # The number in each element of the matrix describes what order to plot the next
              # four plots in. Heatmap.2 plots its elements in the following order:

              # 1. Heatmap,
              # 2. Row dendrogram,
              # 3. Column dendrogram,
              # 4. Key

              # so the default lmat is:
              # > rbind(4:3,2:1)
              #      [,1] [,2]
              # [1,]    4    3
              # [2,]    2    1

              # If for example, you want to put the key underneath the heatmap you would specify:
              lmat = rbind(c(0,3),c(2,1),c(0,4)),
              # > lmat
              #      [,1] [,2]
              # [1,]    0    3
              # [2,]    2    1
              # [3,]    0    4

              # lwid and lhei are vectors that specify the height and width of each row and column.
              # The default is c(1.5,4) for both. If you change lmat you'll either have to or probably
              # want to change these as well. For the above example, if we want to keep all the other
              # elements the same size, but want a thin color key at the bottom, we might set
              lwid = c(0.1, 5.4),
              lhei = c(0.5, 3, 1.5)
    )

    # This will plot a heatmap with the column dendrogram above the heatmap, the row dendrogram
    # to the left, and the key underneath. Unfortunately the headings and the labels for the
    # key are hard coded.
  } else {
    # Original Plot
    heatmap.2(data.matrix(x),
              # dendrogram control
              Rowv=NULL, Colv=NULL,
              dendrogram='none',

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
              key.title="Frequency Distribution of Number of SNPs Compared between Samples",
              key.xlab = "# of SNPs",
              key.ylab = "# of Samples Compared",
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
              lwid=c(1, 10, 1)
    )
  }
}

# ----- END plot function ----- #

topMatchesFile <- gsub(".txt", ".topmatches.txt", cohort_matcher_results)
paste("Writing top matches in", topMatchesFile, sep=" ")
reportTopMatches(table, topMatchesFile)

pdfFile <- gsub(".txt", ".pdf", cohort_matcher_results)
paste("Writing", pdfFile, sep=" ")
pdf(pdfFile)
plotSampleSimilarity(as.matrix(table))
plotNumSNPsCompared(total_compared)
dev.off()

plotFile <- gsub(".txt", ".plot1.tiff", cohort_matcher_results)
tiff(plotFile)
plotSampleSimilarity(as.matrix(table))
dev.off()
plotFile <- gsub(".txt", ".plot2.tiff", cohort_matcher_results)
tiff(plotFile)
plotNumSNPsCompared(total_compared)
dev.off()

