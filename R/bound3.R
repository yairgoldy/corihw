#' bound3: Bound on a uniuon of events
#'
#'
#' Compute a bound that the probability max_j |X_j|>PhiInv(1-alpha/2)
#' where X=(X_1,...,X_m)~MVN(0,Sigma)
#' and PhiInv is the inverse of the standard normal distribution
#'
#'
#'
#' @param Sigma   A sparse mxm corrolation matrix
#' @param m       Dimension of the matrix Sigma
#' @param m_block The matrix Sigma is treated as a block diagonal
#'  with m_block by m_block blocks. All elements outside the block are treated as zero.
#'
#' @return A function which for a vector of alpha return the bound
#' @seealso ??
#'
#' @examples
#' library(Matrix)
# m <- 1000
# Sigma <-  Matrix(runif(m^2),m, m,sparse=T)
# m_block <- 100
# f <- bound3(Sigma,m,m_block)
# alpha <- seq(0,0.001,by = 0.00001)
# plot(alpha,f(alpha),type="l")
#'
#' @importFrom Matrix Matrix
#' @import tibble
#' @import dplyr
#' @export
bound3 <- function(Sigma,m=NULL,m_block=100){

  if(is.null(m)) m <- nrow(Sigma)
  Sigma <- Matrix(Sigma,sparse = T)

  if(m_block>=m) m_block <- m

  # Calculate number of blocks on the diagonal
  nsubmats <- ceiling(m/m_block)

  # Create a table of the blocks (blks) and their indices (inds)
  subinds <- tibble(blks=1:nsubmats) %>% rowwise() %>%
    mutate(inds=list( seq(m_block*(blks-1)+1, min(m_block*blks,m),by = 1)))

  rho <- unique(gamma_rho_pairs$rho)
  #rho <- rho[which(rho!=0)]
  diff <- rho[2]-rho[1]

  # bound3 is based on counting the number of 2x2 matrices [1,rho; rho;1] and 3x3 matrices with off-diagonal elements (rho1,rho2,rho3)

  # Holds the 2x2 matrix-type counts as a function of rho
  rho_pairs <- tibble(rho,summands=0)

  # Holds the 3x3 matrix-type counts as a function of rho1,rho2,rho3, order does not matter
  rho_triples <- expand.grid(rho,rho,rho) %>% as_tibble() %>%
    rename(rho1=Var1,rho2=Var2,rho3=Var3) %>% dplyr::filter(rho1<=rho2, rho2<=rho3) %>% mutate(summands=0)
  rholist <- list(rho_pairs=rho_pairs,rho_triples=rho_triples)

  # Performs the counts using bound3_matcount
  for( blk in subinds$blks){
    inds <- subinds$inds[[blk]]
    Sigma_local <- Sigma[inds,inds]
    rholist <- bound3_matcount(rho_pairs = rholist$rho_pair,
                               rho_triples = rholist$rho_triples,
                               Sigma = Sigma_local,m = m_block,diff = diff)
  }

  # Add the 2x2 matrices with rho=0
  rho_pairs <- rholist$rho_pairs
  num2by2mast <- choose(m,2)
  counted2by2mats <-  dplyr::filter(rho_pairs, rho!= 0) %>% summarise( sums=sum(summands)) %>% select(sums) %>%  unlist(use.names = F)
  rho_pairs <- mutate(rho_pairs,
                      summands= ifelse(rho==0,num2by2mast-counted2by2mats,summands))

  # Add the 3x3 matrices with rho1=rho2=rho3=0
  rho_triples <- rholist$rho_triples
  num3by3mast <- choose(m,3)
  counted3by3mats <-  dplyr::filter(rho_triples,  rho1 > 0 | rho2 > 0 | rho3 >0 ) %>% summarise( sums=sum(summands)) %>% select(sums) %>% unlist(use.names = F)
  rho_triples <- mutate(rho_triples,
                        summands= ifelse(rho1 == 0 & rho2 == 0  & rho3 == 0,num3by3mast-counted3by3mats,summands))

  # Compute S1, S2 and S3 which are sum_j P(|X_j|>gamma), sum_jk P(|X_j|>gamma,|X_k|>gamma) and
  # sum_jkl P(|X_j|>gamma,|X_k|>gamma,|X_l>gamma) where gamma=PhiInv(1-alpha/2)
  gamma <- unique(gamma_rho_pairs$gamma)
  S1 <- m*2*pnorm(-gamma)
  S2 <- inner_join(gamma_rho_pairs,rho_pairs,by="rho") %>%
    group_by(gamma) %>% summarise(pr=sum(summands*prob)) %>% select(pr) %>% unlist()
  S3 <-  inner_join(gamma_rho_triples,rho_triples,by=c("rho1","rho2","rho3")) %>%
    group_by(gamma) %>% summarise(pr=sum(summands*prob))%>% select(pr) %>% unlist()



  bounds <- tibble(gamma,S1,S2,S3) %>%
    mutate(alpha=gamma2alpha(gamma),
           bound1=pmin(1,S1-S2+S3,S1-(2*m-3)*S2/choose(m,2)+(m-2)*S3/choose(m,3) ),
           beta1 = S1,
           beta2= S1+2*S2,
           beta3= S1+6*S2+6*S3,
           k=floor((beta3-beta2)/(beta2-beta1)),
           cond1= (k>=2 & k<m),
           cond2= (k*(k+1)*beta1-(2*k+1)*beta2+beta3>=0),
           cond3= (k*(k+1)-(2*k+1+k*(k+1))*beta1+(2*k+2)*beta2-beta3 >=0),
           allcond= cond1 & cond2 & cond3,
           bound2=pmin(1,1-( k*(k+1)*(1-beta1) + (2*k+1)*(beta2-beta1) -(beta3-beta2) )/(k^2+k) ),
           bound=ifelse(allcond,bound2,bound1)  )%>%
  #          bound=pmin(1,S1-S2+S3,S1-(2*m-3)*S2/choose(m,2)+(m-2)*S3/choose(m,3) )) %>%
  select(alpha,bound)
  bounds <- add_row(bounds,alpha=c(0,1),bound=c(0,1))
  return(approxfun(bounds$alpha,bounds$bound,method = "linear"))
}

#' @importFrom Matrix Matrix triu
#' @import tibble
#' @import dplyr
bound3_matcount <- function(rho_pairs,rho_triples,Sigma,m,diff){

  if(m<3){
    return(list(rho_pairs=rho_pairs,rho_triples=rho_triples))
  }

  # Prepare Sigma to be in a triplet representation
  # Note that rho is calculated with respect to the square of the correlation!!!
  Sigma <- Matrix(round(floor(Sigma^2/diff+1e-10)*diff,digits=1)) %>%  triu(k = 1) %>% slam::as.simple_triplet_matrix()
  Sigma <- tibble(i=Sigma$i,j=Sigma$j,rho=Sigma$v)


  possible_pairs <- expand.grid(1:m,1:m) %>% as_tibble() %>%
    rename(i=Var1,j=Var2) %>% dplyr::filter(i<j)

  #Sigma <- inner_join(possible_pairs,Sigma,by=c("i","j"))
  Sigma <- left_join(possible_pairs,Sigma,by=c("i","j")) %>% mutate(rho=ifelse(is.na(rho),0,rho))

  # Sum 2x2 matrices
  local_rho_pairs <- group_by(Sigma,rho) %>% summarise(summands=n()) # vector of the number of each rho
  rho_pairs <- left_join(rho_pairs,local_rho_pairs,by="rho") %>%
    mutate(summands=ifelse(is.na(summands.y),summands.x,summands.x+summands.y)) %>%
    select(-summands.x,-summands.y)


  # Sum 3x3 matrices. Note that a 3x3 matrix is defined by its 3 indices i,j,k
  Sigma_ij <- Sigma %>% rename(rho_ij=rho)
  Sigma_ik <- Sigma %>% rename(rho_ik=rho,k=j)
  Sigma_jk <- Sigma %>% rename(rho_jk=rho,k=j,j=i)

  all_triples <- inner_join(Sigma_ij,Sigma_ik,by="i") %>%  dplyr::filter(j<k) %>%
    inner_join(Sigma_jk,by=c("j","k")) %>% select(i,j,k,everything())  # Create all 3x3 submatrices

  #order the rhos such that rho1<=rho2<=rho3
  sorted_triples <- mutate(all_triples,rho1=pmin(rho_ij,rho_ik,rho_jk),
                           rho2=ifelse(rho_ij<=rho_ik & rho_ij<= rho_jk,ifelse(rho_ik<=rho_jk, rho_ik,rho_jk),
                                       ifelse(rho_jk<=rho_ik & rho_jk<= rho_ij,ifelse(rho_ij<=rho_ik, rho_ij,rho_ik),
                                              ifelse(rho_ij<=rho_jk, rho_ij,rho_jk))),
                           rho3=pmax(rho_ij,rho_ik,rho_jk)) %>%  select(rho1,rho2,rho3) %>%
    group_by(rho1,rho2,rho3) %>% summarise(summands=n())

  # Add to the original table of 3x3 matrices
  rho_triples <- left_join(rho_triples,sorted_triples,by=c("rho1","rho2","rho3")) %>%
    mutate(summands=ifelse(is.na(summands.y),summands.x,summands.x+summands.y)) %>%
    select(-summands.x,-summands.y)





  return(list(rho_pairs=rho_pairs,rho_triples=rho_triples))

}

