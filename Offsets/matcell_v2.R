matcell <- function(x, tilesize){
  
  sz = c(nrow(x), ncol(x))
  A = vector('list', 2)
  s = 1
  for (s in 1:2){
    dm = sz[s]
    T = min(dm, tilesize[s])
    nn = floor( dm / T ) 
    resid = dm %% tilesize[s]
    A[[s]]=c(rep(1, nn)*T,resid) 
  }

rowsizes = A[[1]];
colsizes = A[[2]];
rows = length(rowsizes)
cols = length(colsizes)
B = vector('list', rows*cols)

rowStart = 0

a = 1

for (i in 1:rows){
  colStart = 0
  for (j in 1:cols){
    B[[a]] = x[colStart+(1:colsizes[j]), rowStart+(1:rowsizes[i])]
    colStart = colStart + colsizes[j]
  a = a + 1
  }
  rowStart = rowStart + rowsizes[i]
}
return(B)
}
