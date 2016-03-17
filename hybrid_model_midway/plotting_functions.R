# plotting functions
# Function to make traceplot ############################################
make_traceplot <- function(results, param ,iter,nChains){
  test <- results[,names(results) == param]
  cut_point <- iter - 100
  df <- data.frame(
    matrix(test, ncol = nChains)
  )
  chain_names <- paste0("chain_",c(1:nChains))
  names(df) <- chain_names
  df$iter = c(1:nrow(df))
  dfm <- melt(df, id.vars = "iter")
  title <- paste("traceplot ", param, sep = "")
  traceplot <- ggplot(dfm, aes(x=iter, y = value, color = variable)) + geom_line()+ ggtitle(title) + theme_bw()
  return(traceplot)
}  

### Function to plot multiple ggplot plots on one page
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}