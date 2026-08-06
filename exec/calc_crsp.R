rescrps01 = matrix(nrow=3,ncol=A)
for(y in 1:3)
  for(a in 1:A)
  {
    print(c(y,a))
    A_mat = matrix(res0[,y,a],ncol=1)
    rescrps01[y,a] = mean(crps_sample(y = res1[,y,a], dat = A_mat))
  }

rescrps02 = matrix(nrow=3,ncol=A)
for(y in 1:3)
  for(a in 1:A)
  {
    print(c(y,a))
    A_mat = matrix(res0[,y,a],ncol=1)
    rescrps02[y,a] = mean(crps_sample(y = res2[,y,a], dat = A_mat))
  }

rescrps12 = matrix(nrow=3,ncol=A)
for(y in 1:3)
  for(a in 1:A)
  {
    print(c(y,a))
    A_mat = matrix(res1[,y,a],ncol=1)
    rescrps12[y,a] = mean(crps_sample(y = res2[,y,a], dat = A_mat))
  }