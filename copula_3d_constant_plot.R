library(copula) 
v1 <- runif(10000)
v2 <- runif(10000)
v3 <- runif(10000)
v <- cbind(v1,v2,v3)
C <- pcopula(norm.cop, v)
j=0
for (i in 1:length(C)) {
if (C[i] >= 0.1 & C[i] < 0.11) { j = j+ 1; }
} 
vj <- cbind(c(1:j),c(1:j),c(1:j))
j=0
for (i in 1:length(C)) {
if (C[i] >= 0.1 & C[i] < 0.11) { j = j+ 1; vj[j,] <- cbind(v1[i],v2[i],v3[i]) }
} 
scatterplot3d(vj, pch=20)
