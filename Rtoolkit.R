### My toolkit of generic data manipulation functions

deduplicate <- function(x, sep="_"){
# Deduplicates a vector - any values X in the vector that are the nth duplicate are replaced by X{sep}n
# User should take note of any values that are already called X{sep}n - dangerous!
  for(i in which(!duplicated(x))){
    count=0
    for(j in i:length(x)){
      if(j!=i){
        if(x[j] == x[i]){
          x[j] = paste0(x[j], sep, count)
          count = count+1
        }
      }
    }
  }
  return(x)
}

qplot3d <- function(mat3d, ...){
  # Quickly plots a 3 dimensional matrix
  # Does not add any colors, etc. by default other than those already included in plotly
  library(plotly)
  p <- plot_ly(as.data.frame(mat3d), x = ~mat3d[,1], y = ~mat3d[,2], z = ~mat3d[,3], ...) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'x'),
                        yaxis = list(title = 'y'),
                        zaxis = list(title = 'z')))
  return(p)
}

  
interpolate.matrix <- function(x){
  cm=colMeans(x, na.rm=T)
  rm=rowMeans(x, na.rm=T)
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if(is.na(x[i,j])){
        n=mean(rm[i], cm[j])
        x[i,j] = n
      }
    }
  }
  return(x)
}

reorder2 <- function(X, x.1=NULL, by){
  if(is.null(x.1) & is.null(dim(X))){
    return(X[order(match(X, by))])
  }
  
  return(X[order(match(x.1, by)),])
}

