## This file creates the tables norm_table.Rdat and norm_table_squares.Rdat
## Each contains a table of P(min|X_1|,|X_2|,|X_3|>gamma) where X~N(0,Sigam)
## where Sigma has 1 on the diagonal and r_1,r_2,r_3 off-diagonal



rm(list=ls())
library(tidyverse)
library(mvtnorm)

s2 <- function(rho,gamma){
  Sigma <- matrix(c(1,rho,rho,1),2,2)
  lower1 <- c(-Inf,-Inf)
  upper1 <- c(-gamma,-gamma)
  lower2 <- c(-Inf,gamma)
  upper2 <- c(-gamma,Inf)
  return(2*pmvnorm(lower1,upper1,corr=Sigma)+2*pmvnorm(lower2,upper2,corr=Sigma))
}

s3 <- function(Sigma,gamma){
  lower1 <- c(-Inf,-Inf,-Inf)
  upper1 <- c(-gamma,-gamma,-gamma)
  lower2 <- c(-Inf,-Inf,gamma)
  upper2 <- c(-gamma,-gamma,Inf)
  lower3 <- c(-Inf,gamma,-Inf)
  upper3 <- c(-gamma,Inf,-gamma)
  lower4 <- c(gamma,-Inf,-Inf)
  upper4 <- c(Inf,-gamma,-gamma)
  return(2*pmvnorm(lower1,upper1,corr=Sigma)+
           2*pmvnorm(lower2,upper2,corr=Sigma)+
           2*pmvnorm(lower3,upper3,corr=Sigma)+
           2*pmvnorm(lower4,upper4,corr=Sigma))
}

gamma <- seq(1,10,by=0.01)
rho <- (0:9)/10
pairs <- as_tibble(expand.grid(gamma,rho))
names(pairs) <- c("gamma","rho")
triples <- as_tibble(expand.grid(gamma,rho,rho,rho))
names(triples) <- c("gamma","rho1","rho2","rho3")
triples <- filter(triples,rho1<=rho2,rho1<=rho3,rho2<=rho3)




pairs$prob <- 1
for(i in 1:nrow(pairs)){
  pairs$prob[i] <- s2(pairs$rho[i],pairs$gamma[i])
}
pairs <- mutate(pairs,alpha=2*pnorm(-gamma))
save(pairs,file="norm_table.Rdat")

triples$prob <- 1
for(i in 1:nrow(triples)){
  rho12 <- triples$rho1[i]
  rho13 <- triples$rho2[i]
  rho23 <- triples$rho3[i]
  
  Sigma <- matrix(c(1,rho12,rho13,rho12,1,rho23,rho13,rho23,1),3,3)
  triples$prob[i] <- s3(Sigma,triples$gamma[i])
}
triples <- mutate(triples,alpha=2*pnorm(-gamma))


gamma_rho_pairs <- pairs
gamma_rho_triples <- triples
save(gamma_rho_pairs,gamma_rho_triples,file="norm_table.Rdat")
