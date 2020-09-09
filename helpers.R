
# given list of adjusted p-values and original p-values
# calculate actual threshold!


mydiv <- function(x,y) ifelse(x == 0, 0,
                              ifelse(y==0, 1, pmin(x/y,1)))


presorted_grenander <- function(sorted_pvalues, m_total=length(sorted_pvalues),
                                quiet=TRUE){

  unique_pvalues <- unique(sorted_pvalues)
  ecdf_values <- cumsum(tabulate(match(sorted_pvalues, unique_pvalues)))/m_total

  if (min(unique_pvalues) > 0){
    # I think fdrtool neglects this borderline case and this causes returned object
    # to be discontinuous hence also not concave
    unique_pvalues <- c(0,unique_pvalues)
    ecdf_values   <- c(0, ecdf_values)
  }

  if (max(unique_pvalues) < 1){
    unique_pvalues <- c(unique_pvalues,1)
    ecdf_values    <- c(ecdf_values,1)
  }

  ll <- fdrtool::gcmlcm(unique_pvalues, ecdf_values, type="lcm")
  ll$length <- length(ll$slope.knots)
  ll$x.knots <- ll$x.knots[-ll$length]
  ll$y.knots <- ll$y.knots[-ll$length]
  # if (!quiet) message(paste("Grenander fit with", ll$length, "knots."))
  ll
}



alpha2gamma <- function(alpha){
  return(-qnorm(alpha/2))

}

gamma2alpha <- function(gamma){
  return(2*pnorm(-gamma))
}


prepare_dat <- function(){
  pvalues <- deRes$pvalue
  covariates <- deRes$baseMean
  dat_original <- tibble(id=1:length(covariates),pvalues,covariates)
  dat <-  drop_na(dat_original)
  dat <- mutate(dat,r_pvalues= ifelse(pvalues > 10^(-20), pvalues, 0),i=1:nrow(dat))
  return(dat)
}


#' @importFrom Matrix Matrix sparseMatrix
as.sparseMatrix <- function(simple_triplet_matrix_sparse) {
  retval <-  sparseMatrix(i=as.numeric(simple_triplet_matrix_sparse$i),
                          j=as.numeric(simple_triplet_matrix_sparse$j),
                          x=as.numeric(as.character(simple_triplet_matrix_sparse$v)),
                          dims=c(simple_triplet_matrix_sparse$nrow,
                                 simple_triplet_matrix_sparse$ncol),
                          dimnames = dimnames(simple_triplet_matrix_sparse),
                          giveCsparse = TRUE)
}

calc_m_effecitve <- function(Sigma,inds=1:nrow(Sigma)){
  is=unlist(inds)
  S <- Sigma[is,is]
  S <- as(S,"dgCMatrix")

  neigs <- min(100,nrow(S))
  cond <- TRUE
  while(cond){
    if(neigs==nrow(S)){
      v <- eigen(S)
    } else{
      v <- rARPACK::eigs(S,neigs)
    }
    eigvals <- v$values[v$values>1]


    if(length(eigvals)==length(v$value)) {
      neigs <- min(neigs+100,nrow(S))
    } else {
      cond <- FALSE
    }
  }
  m_effective <- nrow(S)-sum(eigvals)+length(eigvals)
  return(m_effective)
}


#' @importFrom dplyr filter
#' @importFrom slam simple_triplet_matrix
linear_constraints <- function(dat,training_indices,lambda,m,m_groups,nbins,penalty,quiet){
  if (!quiet) message("Applying Grenander estimator within each bin.")

  grenander_list  <- dplyr::filter(dat,training_indices) %>% arrange(groups,r_pvalues) %>% group_by(groups) %>%
    summarise(F=list(presorted_grenander(r_pvalues,quiet=quiet)))
  grenander_list <-  grenander_list$F
  #set up LP
  nconstraints_per_bin <- sapply(grenander_list, function(x) x$length)
  messsage(paste("nconstraints_per_bin = ",nconstraints_per_bin))

  nconstraints <- sum(nconstraints_per_bin)
  i_yi <- 1:nconstraints
  j_yi <- rep(1:nbins,  times=nconstraints_per_bin)
  v_yi <- rep(1,  nconstraints)

  i_ti <- 1:nconstraints
  j_ti <- nbins + rep(1:nbins, times=nconstraints_per_bin)
  v_ti <- unlist(lapply(grenander_list, function(x) -x$slope.knots))

  constr_matrix <- slam::simple_triplet_matrix(c(i_yi, i_ti), c(j_yi,j_ti), c(v_yi, v_ti))
  rhs <- unlist(lapply(grenander_list, function(x) x$y.knots-x$slope.knots*x$x.knots))

  obj <- c(m_groups/m*nbins*rep(1,nbins), rep(0,nbins))

  if (lambda < Inf){
    if (penalty == "total variation"){
      # -f + t_g - t_{g-1} <= 0
      i_fi <- rep(1:(nbins-1),3)
      j_fi <-	c((nbins+2):(2*nbins), (nbins+1):(2*nbins-1), (2*nbins+1):(3*nbins-1))
      v_fi <- c(rep(1,nbins-1), rep(-1,nbins-1), rep(-1, nbins-1))
      absolute_val_constr_matrix_1 <- slam::simple_triplet_matrix(i_fi,j_fi,v_fi)

      # -f - t_g + t_{g-1} <= 0
      i_fi <- rep(1:(nbins-1),3)
      j_fi <-	c((nbins+2):(2*nbins), (nbins+1):(2*nbins-1), (2*nbins+1):(3*nbins-1))
      v_fi <- c(rep(-1,nbins-1), rep(1,nbins-1), rep(-1, nbins-1))
      absolute_val_constr_matrix_2 <- slam::simple_triplet_matrix(i_fi,j_fi,v_fi)

      constr_matrix <- rbind(cbind(constr_matrix, slam::simple_triplet_zero_matrix(nconstraints, nbins-1,mode="double")),
                             absolute_val_constr_matrix_1,
                             absolute_val_constr_matrix_2)

      obj <- c(obj, rep(0, nbins-1))

      total_variation_constr <- matrix(c(rep(0,nbins), -lambda*m_groups/m, rep(1,nbins-1)),nrow=1)
      constr_matrix <- rbind(constr_matrix, total_variation_constr)

      rhs <- c(rhs,rep(0, 2*(nbins-1)),0) #add RHS for absolute differences and for total variation penalty

    } else if (penalty == "uniform deviation"){
      # -f + m t_g - sum_g m_i t_i <= 0

      absolute_val_constr_matrix_1 <- matrix(rep(-m_groups,nbins), nbins,nbins, byrow=TRUE)
      diag(absolute_val_constr_matrix_1) <- -m_groups + m
      absolute_val_constr_matrix_1 <- cbind( slam::simple_triplet_zero_matrix(nbins, nbins,mode="double"),
                                             absolute_val_constr_matrix_1,
                                             -diag(nbins))

      # -f - m t_g +  sum_g m_i t_i  <= 0

      absolute_val_constr_matrix_2 <- matrix(rep(+m_groups,nbins), nbins,nbins, byrow=TRUE)
      diag(absolute_val_constr_matrix_2) <- m_groups - m
      absolute_val_constr_matrix_2 <- cbind( slam::simple_triplet_zero_matrix(nbins, nbins,mode="double"),
                                             absolute_val_constr_matrix_2,
                                             -diag(nbins))

      constr_matrix <- rbind(cbind(constr_matrix, slam::simple_triplet_zero_matrix(nconstraints, nbins,mode="double")),
                             absolute_val_constr_matrix_1,
                             absolute_val_constr_matrix_2)

      obj <- c(obj, rep(0, nbins))

      total_variation_constr <- matrix(c(rep(0,nbins), -lambda*m_groups, rep(1,nbins)),nrow=1)
      constr_matrix <- rbind(constr_matrix, total_variation_constr)

      rhs <- c(rhs,rep(0, 2*nbins),0) #add RHS for absolute differences and for total variation penalty

    }
  }
  return(list(constr_matrix=constr_matrix,rhs=rhs,obj=obj))
}

#' @import dplyr
fwer_lin_constr <- function(funs,alpha,const=2){
  nbins <- nrow(funs)
  m <- sum(funs$ms)
  t0 <- const*alpha/m
  calc_sl <- function(f,x){ ( f(x)-f(x-1e-10) )/1e-10}
  funs <-  funs %>% rowwise() %>% mutate(y0=f(t0),sl=calc_sl(f,t0)) #%>% select(y0) %>% unlist(use.names = F)

  fwer_constr <- c(rep(0,nbins),funs$sl)
  rhs <- alpha-sum(funs$y0)+t0*sum(funs$sl)
  return(list(fwer_constr=fwer_constr,rhs=rhs))

}


# given list of adjusted p-values and original p-values
# calculate actual threshold!
mydiv <- function(x,y) ifelse(x == 0, 0,
                              ifelse(y==0, 1, pmin(x/y,1)))

#' @importFrom fdrtool gcmlcm
presorted_grenander <- function(sorted_pvalues, m_total=length(sorted_pvalues),
                                quiet=TRUE){

  unique_pvalues <- unique(sorted_pvalues)
  ecdf_values <- cumsum(tabulate(match(sorted_pvalues, unique_pvalues)))/m_total

  if (min(unique_pvalues) > 0){
    # I think fdrtool neglects this borderline case and this causes returned object
    # to be discontinuous hence also not concave
    unique_pvalues <- c(0,unique_pvalues)
    ecdf_values   <- c(0, ecdf_values)
  }

  if (max(unique_pvalues) < 1){
    unique_pvalues <- c(unique_pvalues,1)
    ecdf_values    <- c(ecdf_values,1)
  }

  ll <- fdrtool::gcmlcm(unique_pvalues, ecdf_values, type="lcm")
  ll$length <- length(ll$slope.knots)
  ll$x.knots <- ll$x.knots[-ll$length]
  ll$y.knots <- ll$y.knots[-ll$length]
  # if (!quiet) message(paste("Grenander fit with", ll$length, "knots."))
  ll
}

#' @importFrom stats qnorm
alpha2gamma <- function(alpha){
  return(-qnorm(alpha/2))

}

#' @importFrom stats pnorm
gamma2alpha <- function(gamma){
  return(2*pnorm(-gamma))
}

#' @import tibble
#' @import dplyr
prepare_dat <- function(){
  pvalues <- deRes$pvalue
  covariates <- deRes$baseMean
  dat_original <- tibble(id=1:length(covariates),pvalues,covariates)
  dat <-  drop_na(dat_original)
  dat <- mutate(dat,r_pvalues= ifelse(pvalues > 10^(-20), pvalues, 0),i=1:nrow(dat))
  return(dat)
}


#' @importFrom Matrix sparseMatrix
as.sparseMatrix <- function(simple_triplet_matrix_sparse) {
  retval <-  sparseMatrix(i=as.numeric(simple_triplet_matrix_sparse$i),
                          j=as.numeric(simple_triplet_matrix_sparse$j),
                          x=as.numeric(as.character(simple_triplet_matrix_sparse$v)),
                          dims=c(simple_triplet_matrix_sparse$nrow,
                                 simple_triplet_matrix_sparse$ncol),
                          dimnames = dimnames(simple_triplet_matrix_sparse),
                          giveCsparse = TRUE)
}

#' @importFrom  rARPACK eigs
#' @importFrom methods as
calc_m_effecitve <- function(Sigma,inds=1:nrow(Sigma)){
  is=unlist(inds)
  S <- Sigma[is,is]
  S <- as(S,"dgCMatrix")

  neigs <- min(100,nrow(S))
  cond <- TRUE
  while(cond){
    if(neigs==nrow(S)){
      v <- eigen(S)
    } else{
      v <- rARPACK::eigs(S,neigs)
    }
    eigvals <- v$values[v$values>1]


    if(length(eigvals)==length(v$value)) {
      neigs <- min(neigs+100,nrow(S))
    } else {
      cond <- FALSE
    }
  }
  m_effective <- nrow(S)-sum(eigvals)+length(eigvals)
  return(m_effective)
}

#' @importFrom slam simple_triplet_matrix
#' @importFrom dplyr filter
linear_constraints <- function(dat,training_indices,lambda,m,m_groups,nbins,penalty,quiet){
  if (!quiet) message("Applying Grenander estimator within each bin.")

  grenander_list  <- dplyr::filter(dat,training_indices) %>% arrange(groups,r_pvalues) %>% group_by(groups) %>%
    summarise(F=list(presorted_grenander(r_pvalues,quiet=quiet)))
  grenander_list <-  grenander_list$F

  #set up LP
  nconstraints_per_bin <- sapply(grenander_list, function(x) x$length)
  nconstraints <- sum(nconstraints_per_bin)
  i_yi <- 1:nconstraints
  j_yi <- rep(1:nbins,  times=nconstraints_per_bin)
  v_yi <- rep(1,  nconstraints)

  i_ti <- 1:nconstraints
  j_ti <- nbins + rep(1:nbins, times=nconstraints_per_bin)
  v_ti <- unlist(lapply(grenander_list, function(x) -x$slope.knots))

  constr_matrix <- slam::simple_triplet_matrix(c(i_yi, i_ti), c(j_yi,j_ti), c(v_yi, v_ti))
  rhs <- unlist(lapply(grenander_list, function(x) x$y.knots-x$slope.knots*x$x.knots))

  obj <- c(m_groups/m*nbins*rep(1,nbins), rep(0,nbins))

  if (lambda < Inf){
    if (penalty == "total variation"){
      # -f + t_g - t_{g-1} <= 0
      i_fi <- rep(1:(nbins-1),3)
      j_fi <-	c((nbins+2):(2*nbins), (nbins+1):(2*nbins-1), (2*nbins+1):(3*nbins-1))
      v_fi <- c(rep(1,nbins-1), rep(-1,nbins-1), rep(-1, nbins-1))
      absolute_val_constr_matrix_1 <- slam::simple_triplet_matrix(i_fi,j_fi,v_fi)

      # -f - t_g + t_{g-1} <= 0
      i_fi <- rep(1:(nbins-1),3)
      j_fi <-	c((nbins+2):(2*nbins), (nbins+1):(2*nbins-1), (2*nbins+1):(3*nbins-1))
      v_fi <- c(rep(-1,nbins-1), rep(1,nbins-1), rep(-1, nbins-1))
      absolute_val_constr_matrix_2 <- slam::simple_triplet_matrix(i_fi,j_fi,v_fi)

      constr_matrix <- rbind(cbind(constr_matrix, slam::simple_triplet_zero_matrix(nconstraints, nbins-1,mode="double")),
                             absolute_val_constr_matrix_1,
                             absolute_val_constr_matrix_2)

      obj <- c(obj, rep(0, nbins-1))

      total_variation_constr <- matrix(c(rep(0,nbins), -lambda*m_groups/m, rep(1,nbins-1)),nrow=1)
      constr_matrix <- rbind(constr_matrix, total_variation_constr)

      rhs <- c(rhs,rep(0, 2*(nbins-1)),0) #add RHS for absolute differences and for total variation penalty

    } else if (penalty == "uniform deviation"){
      # -f + m t_g - sum_g m_i t_i <= 0

      absolute_val_constr_matrix_1 <- matrix(rep(-m_groups,nbins), nbins,nbins, byrow=TRUE)
      diag(absolute_val_constr_matrix_1) <- -m_groups + m
      absolute_val_constr_matrix_1 <- cbind( slam::simple_triplet_zero_matrix(nbins, nbins,mode="double"),
                                             absolute_val_constr_matrix_1,
                                             -diag(nbins))

      # -f - m t_g +  sum_g m_i t_i  <= 0

      absolute_val_constr_matrix_2 <- matrix(rep(+m_groups,nbins), nbins,nbins, byrow=TRUE)
      diag(absolute_val_constr_matrix_2) <- m_groups - m
      absolute_val_constr_matrix_2 <- cbind( slam::simple_triplet_zero_matrix(nbins, nbins,mode="double"),
                                             absolute_val_constr_matrix_2,
                                             -diag(nbins))

      constr_matrix <- rbind(cbind(constr_matrix, slam::simple_triplet_zero_matrix(nconstraints, nbins,mode="double")),
                             absolute_val_constr_matrix_1,
                             absolute_val_constr_matrix_2)

      obj <- c(obj, rep(0, nbins))

      total_variation_constr <- matrix(c(rep(0,nbins), -lambda*m_groups, rep(1,nbins)),nrow=1)
      constr_matrix <- rbind(constr_matrix, total_variation_constr)

      rhs <- c(rhs,rep(0, 2*nbins),0) #add RHS for absolute differences and for total variation penalty

    }
  }
  return(list(constr_matrix=constr_matrix,rhs=rhs,obj=obj))
}

fwer_lin_constr_old <- function(funs,alpha,const=2){
  nbins <- nrow(funs)
  m <- sum(funs$ms)
  t0 <- const*alpha/m
  calc_sl <- function(f,x){ ( f(x)-f(x-1e-10) )/1e-10}
  funs <-  funs %>% rowwise() %>% mutate(y0=f(t0),sl=calc_sl(f,t0)) #%>% select(y0) %>% unlist(use.names = F)

  fwer_constr <- c(rep(0,nbins),funs$sl)
  rhs <- alpha-sum(funs$y0)+t0*sum(funs$sl)
  return(list(fwer_constr=fwer_constr,rhs=rhs))

}

fwer_lin_constr <- function(funs,alpha,const=2){
  nbins <- nrow(funs)
  m <- sum(funs$ms)
  x0 <- const*alpha/m
  xs <- seq(0.001,0.1,by = 1e-5)
  linear_f <- function(x0, y0, slope0, x){return(slope0*(x-x0)+y0)}
  calc_sl <- function(f,x){ ( f(x)-f(x-1e-10) )/1e-10}

  funs <-  dplyr::rowwise(data = funs) %>%
    dplyr::mutate(y0=f(x0),
                  sl=calc_sl(f,x0)) %>%
    dplyr::mutate(lf = purrr::map2(sl ,y0, function(sli, y0i){purrr::partial(.f = linear_f, x0 = x0, slope0 =sli, y0 = y0i)})) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(fdiff = min(lf(xs)-f(xs))) %>%
    dplyr::mutate( y0 = y0 -  pmin( fdiff,0)) %>%
    ungroup()


    fwer_constr <- c(rep(0,nbins),funs$sl)
  rhs <- alpha-sum(funs$y0)+x0*sum(funs$sl)
  return(list(fwer_constr=fwer_constr,rhs=rhs))

}
