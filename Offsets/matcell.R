matcell <- function(x, rowsizes, colsizes){

rows = length(rowsizes)
cols = length(colsizes)
c = vector('list', rows*cols)

rowStart = 0

a = 1

for (i in 1:rows){
  colStart = 0
  for (j in 1:cols){
    c[[a]] = x[rowStart+(1:rowsizes[i]),colStart+(1:colsizes[j])]
    colStart = colStart + colsizes[j]
  a <- a + 1
  }
  rowStart = rowStart + rowsizes[i]
}
return(c)
}
