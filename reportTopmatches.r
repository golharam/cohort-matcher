table <- read.table("cohort-matcher-results.txt", header=TRUE, row.names=1)
rownames(table) <- colnames(table)

f <- file("topmatches.txt", "w")
writeLines(paste("sample", "match1", "score1", "match2", "score2", "match3", "score3", "match4", "score4", "match5", "score5", sep="\t"), f)
for (sample in rownames(table)) {
	sortedTable <- table[order(-table[,sample]),]
	writeLines(paste(sample, 
		    rownames(sortedTable)[1],  sortedTable[1,sample],
		    rownames(sortedTable)[2],  sortedTable[2,sample],
		    rownames(sortedTable)[3],  sortedTable[3,sample],
		    rownames(sortedTable)[4],  sortedTable[4,sample],
		    rownames(sortedTable)[5],  sortedTable[5,sample],
		    sep="\t"),
		   f)
}
close(f)

# Plot
# http://www.phaget4.org/R/image_matrix.html
# https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
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
# ----- END plot function ----- #
M.table <- as.matrix(table)
myImagePlot(M.table)#
