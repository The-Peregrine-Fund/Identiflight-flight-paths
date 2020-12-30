load(file.chosse())
# functions
min.na <- function(x){
  na <- min(which(is.na(x)))
  y <- na-1
  return(y)
}
min.data <- function(x){
  y <- min(which(x==0 | x==1))
  return(y)
}
max.data <- function(x){
  y <- max(which(x==0 | x==1))
  return(y)
}
min.na.data <- function(x){
  mx <- max(which(x==0 | x==1))
  wna <- which(is.na(x))
  wna2<- min(wna[wna>mx])
  y <- wna2-1
  return(y)
}

len <- function(x){
  y<- length(which(x==0 | x==1))
  return(y)
  }
l <- apply(datl$z, 1, len)
l <- apply(datl$z, 1, min.na)
mean(l)
sd(l)
range(l)
# view weird tracks
datl$z[l==0,]