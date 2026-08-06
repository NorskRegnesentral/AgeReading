compare_samp = function(mat)
{
  M = nrow(mat)
  p = ncol(mat)
  pmat = matrix(nrow=p,ncol=p)
  for(i in 1:p)
    for(j in i:p)
    {
      delta = mat[,i]-mat[,j]
      m = sum(delta<=0)/M
      pmat[i,j] = 2*min(m,1-m)
      pmat[j,i] = pmat[i,j]
    }
  diag(pmat) = 1
  pmat
}