#' corihw: Main function for Independent Hypothesis Weighting with Correlation
#'
#' Given a vector of p-values, a vector of covariates which are independent of the p-values under the null hypothesis,
#' a matrix Sigma of the correlation of the test statistics and a nominal significance level alpha,
#' IHW Cor learns multiple testing weights and then applies the weighted Bonferroni) procedure.
#'
#'
#' @param pvalues  Numeric vector of unadjusted p-values.
#' @param covariates  Vector which contains the one-dimensional covariates (independent under the H0 of the p-value)
#'                for each test. Assumed continuous.
#' @param alpha   Numeric, sets the nominal level for FWER control.
#' @param nbins  Integer, number of groups into which p-values will be split based on covariate. Use "auto" for
#'             automatic selection of the number of bins. Only applicable when covariates is not a factor.
#' @param quiet  Boolean, if False a lot of messages are printed during the fitting stages.
#' @param nfolds Number of folds into which the p-values will be split for the pre-validation procedure
#' @param lambda  A regularization constant
#' @param seed Integer or NULL. Split of hypotheses into folds is done randomly. To have the output of the function be reproducible,
#'	the seed of the random number generator is set to this value at the start of the function. Use NULL if you don't want to set the seed.
#' @param Sigma   A sparse mxm corrolation matrix
#' @param methods a vector that includes at least one of "IHW", CorIHW", or "M-effective"
#' @param linear_approx A flag to perform CorIHW using linear approximation. This is a much faster procedure which yields similar results.
#'
#' @return A data frame with the weights and a rejection flag per observation
#'
#' @import IHW
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @importFrom stats approxfun pnorm qnorm runif uniroot
#' @export
corihw <- function(pvalues, covariates,alpha,
                   nbins = "auto",
                   quiet = TRUE,
                   nfolds = 5,
                   lambda = 1,
                   seed = 1L,
                   Sigma= NULL,
                   methods = "CorIHW",
                   linear_approx = TRUE){

  bound3_linear_testing <- FALSE
  bound3_testing <- FALSE
  IHW_testing <- FALSE
  meffective_testing <- FALSE

  if("CorIHW"%in% methods && linear_approx == T) bound3_linear_testing <- TRUE
  if("CorIHW"%in% methods && linear_approx == F) bound3_testing <- TRUE

  if("IHW"%in% methods) IHW_testing <- TRUE
  if("M-effective"%in% methods) meffective_testing <- TRUE

  # Create the main table
  dat <- tibble::tibble(id=1:length(pvalues),pvalues,covariates)
  dat <-  tidyr::drop_na(dat)
  dat <- mutate(dat,r_pvalues= ifelse(pvalues > 10^(-20), pvalues, 0))


  nfolds <- as.integer(nfolds)


  if(is.null(Sigma)){
    message("No correlation matrix. Return IHW solution")
    sol <- IHW::ihw (dat$pvalues,dat$covariates,alpha=alpha,distrib_estimator = "grenander",adjustment_type = "bonferroni",
                     covariate_type = "ordinal", nbins = nbins, m_groups = NULL,quiet =quiet,
                     nfolds = nfolds, lambdas = lambda, seed = seed)
    dat <- as_tibble(sol@df) %>%
      mutate(id=1:n()) %>%
      select(id,pvalues=pvalue,covariates=covariate,groups=group,folds=fold) %>%
      mutate(ts=IHW::weights(sol)*alpha(sol)/nrow(sol@df),rjs=rejected_hypotheses(sol))
    return(dat)
  }



  if(all(!bound3_testing,!bound3_linear_testing,!meffective_testing,!IHW_testing)){
    stop("Need to choose at least one testing methods: IHW, CorIHW, or M-effective.")
  }

  # Split to groups
  if (nfolds==1){
    message("Using only 1 fold! Only use this if you want to learn the weights, but NEVER for testing!")
  }
  if (nbins == "auto"){
    nbins <- max(1,min(40, floor(nrow(dat)/1500))) # rule of thumb..
  }
  dat <- mutate(dat, groups = as.factor(IHW::groups_by_filter(dat$covariates, nbins, seed=seed)))

  nbins <- nlevels(dat$groups)



  if (nbins < 1) {
    stop("Cannot have less than one bin.")
  }


  groups_dat <- group_by(dat,groups) %>%  summarise(m_groups=n()) %>% arrange(groups)

  if (nbins > 1 & any(groups_dat$m_groups < 1000)){
    message("We recommend that you supply (many) more than 1000 p-values for meaningful data-driven hypothesis weighting results.")
  }
  # end splitting to groups




  if (!is.null(seed)){
    #http://stackoverflow.com/questions/14324096/setting-seed-locally-not-globally-in-r?rq=1
    tmp <- runif(1)
    old <- .Random.seed
    on.exit( { .Random.seed <<- old } )
    set.seed(as.integer(seed)) #seed
  }



  m <- nrow(dat)

  if (nbins == 1){
    message("Only 1 bin; IHW_COR reduces to bonferroni_COR")
    dat <- mutate(dat,groups=1,folds=1)
    if(bound3_testing || bound3_linear_testing){
      # find the appropriate Bonferroni bound with correlation
      fbound <- bound3(Sigma,m)
      ts <- uniroot(function(x){fbound(x)-alpha}, c(0, alpha),tol=10^(-9))$root
      if(bound3_testing)       dat <- mutate(dat,Cor=ts)
      if(bound3_linear_testing)dat <- mutate(dat,Cor=ts)
    }
    if(meffective_testing){
      meffective <- calc_m_effecitve(Sigma)
      ts <- alpha/meffective
      dat <- mutate(dat,meffective=ts)
    }
    if(IHW_testing)

      ts <- alpha/nrow(dat)
    dat <- mutate(dat,IHW=ts)
  }else if (nbins > 1) {

    #  do the k-fold strategy

    fold_lambdas <- rep(NA, nfolds)

    ##### TODO delete
    #set.seed(seed=seed)

    dat <- mutate(dat,folds=sample(1:nfolds, m, replace = TRUE))

    # Create the table dat_ts to hold the cutoff for the different gourp in each fold
    dat_ts <- expand.grid(1:nfolds,levels(dat$groups)) %>% as_tibble()
    colnames(dat_ts) <- c("folds","groups")
    if(bound3_testing) dat_ts <- mutate(dat_ts,Cor=0)
    if(bound3_linear_testing) dat_ts <- mutate(dat_ts,Cor=0)
    if(meffective_testing) dat_ts <- mutate(dat_ts,meffective=0)
    if(IHW_testing) dat_ts <- mutate(dat_ts,IHW=0)
    dat_ts <-  arrange(dat_ts,folds,groups)



    for (i in 1:nfolds){

      if (!quiet) message(paste("Estimating weights for fold", i))


      fold_lambdas[i] <- lambda

      # ok we have finally picked lambda and can proceed
      if (!quiet) message(paste( "fold=", i,"\n"))
      Sigma_i <- Sigma[dat$folds==i,dat$folds==i]
      tsi <- ihw_cor_convex(dat,dat$folds!=i,penalty=penalty, lambda=lambda,
                            alpha=alpha/nfolds,Sigma_i,quiet=quiet,
                            bound3_testing=bound3_testing,bound3_linear_testing=bound3_linear_testing,
                            meffective_testing=meffective_testing,IHW_testing=IHW_testing)
      if(bound3_testing) dat_ts$Cor[dat_ts$folds==i] <- tsi$Cor
      if(bound3_linear_testing) dat_ts$Cor[dat_ts$folds==i] <- tsi$Cor
      if(meffective_testing) dat_ts$meffective[dat_ts$folds==i] <- tsi$meffective
      if(IHW_testing)  dat_ts$IHW[dat_ts$folds==i] <- tsi$IHW



    }

    dat <- inner_join(dat,dat_ts, by=c("folds","groups"))
  }
  if(bound3_testing) dat <- mutate(dat,rjs_Cor=r_pvalues<Cor)
  if(bound3_linear_testing) dat <- mutate(dat,rjs_Cor=r_pvalues<Cor)
  if(meffective_testing) dat <- mutate(dat,rjs_meffective=r_pvalues<meffective)
  if(IHW_testing) dat <- mutate(dat,rjs_IHW=r_pvalues<IHW)
  dat <- select(dat, - r_pvalues)



  return(dat)


}


ihw_cor_convex <- function(dat,training_indices,penalty="total variation", lambda=Inf, alpha=0.1,Sigma,quiet=quiet,
                           bound3_testing=T,bound3_linear_testing=T,
                           meffective_testing=T,IHW_testing=T){


  # The training indeces are used to estimate the distribution function of the p-values

  # m_groups is vector of the size of groups in the testing data
  testing_indices <- !training_indices
  m_groups <- dplyr::filter(dat,testing_indices) %>% group_by(groups) %>%
    summarise(n=n()) %>% ungroup() %>% arrange(groups) %>% select(n) %>% unlist(use.names = F)
  m <- sum(m_groups)
  nbins <- length(unique(dat$groups))

  constraints <- linear_constraints(dat,training_indices,lambda,m,m_groups,nbins,penalty="total variation",quiet)
  constr_matrix <- constraints$constr_matrix
  rhs <- constraints$rhs
  obj <- constraints$obj


  ts_indices <- (nbins+1):(2*nbins)
  nvars <- ncol(constr_matrix)
  lb <- rep(0,ncol(constr_matrix))
  ub <- rep(1,ncol(constr_matrix))
  x0=c(rep(0,nbins),rep(alpha/m,nbins))
  objective_fun <- function(v){-sum(obj*v)}
  sol <- NULL

  if(ncol(constr_matrix)> 2*nbins){x0 <- c(x0,rep(0,nbins-1))}


  if(bound3_testing || bound3_linear_testing){
    if(!quiet) message("Calculating the probability bound.")
    create_f <- function(inds,ms){ is=unlist(inds); S <- Sigma[is,is];return(bound3(S,ms))}
    funs <- dplyr::filter(dat,testing_indices) %>%  mutate(inds=1:n()) %>% arrange(groups) %>% group_by(groups)  %>%
      summarise(ms=n(),inds=list(inds)) %>% #size and indeces of submatrices
      rowwise() %>%
      mutate(f=list(create_f(inds,ms))) # Create a list of functions that return the probablity P(union_group_g |X_i|> 1-PhiI(ts/2) )

    constr_sparse_matrix <- as.sparseMatrix(constr_matrix)


    # convert linear constraints to functions

    flc <- fwer_lin_constr(funs,alpha)
    fwer_constr <- flc$fwer_constr
    if(ncol(constr_matrix)> 2*nbins){fwer_constr <- c(fwer_constr,rep(0,nbins-1))}
    constr_matrix_lin <- rbind(constr_matrix, matrix(fwer_constr,nrow = 1))
    rhs_lin <- c(rhs,flc$rhs)

    res <- lpsymphony::lpsymphony_solve_LP(obj, constr_matrix_lin, rep("<=", nrow(constr_matrix_lin)),
                                           rhs_lin,
                                           max = TRUE, verbosity = -2, first_feasible = FALSE)

    if(bound3_linear_testing){
      sol <-cbind(sol,Cor=res$solution[ts_indices])
    }

    if(bound3_testing){
      fun_list <- list()
      g <- function(i){force(i);
        function(v){sum(v*constr_sparse_matrix[i,])-rhs[i]}
      }

      for(i in 1:nrow(constr_sparse_matrix)){

        fun_list[[i]]<- g(i)
      }

      # create the generalized Bonferroni bound

      fwer_f <- function(v){
        ts <- v[ts_indices]
        val <- 0
        for(i in 1:length(ts)){
          val <- val+funs$f[[i]](ts[i])
        }
        if(is.na(val)){val <- 1}
        return(100*(val-alpha))
      }
      fun_list <- c(fun_list,fwer_f)

      grad_obj <- function(v){return(-obj)}
      jac_constr <- function(v){
        jac <- rbind(constr_sparse_matrix,
                     nl.jacobian(v,fwer_f))
        return(as.matrix(jac))
      }

      # create one function out of all constraints that returns a vector of length nconstraints+1 (for the fwer constraint)
      constr_fun <- function(v){lapply(fun_list,function(x)do.call(x,list(v))) %>%  unlist()}

      if(!quiet) message("Start non-linear optimization")
      if(!quiet) message(paste("(LINEAR) Constraint: max (should be zero):",max(constr_fun(res$solution)),"max at constr:",
                               which.max(constr_fun(res$solution)),"out of ",nrow(constr_sparse_matrix)))
      res <- nloptr::nloptr( x0= x0,
                             eval_f=objective_fun,
                             #eval_grad_f =  grad_obj,
                             lb = lb,
                             ub = ub,
                             eval_g_ineq = constr_fun,
                             #eval_jac_g_ineq = jac_constr,
                             opts = list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-6,maxeval=100,"print_level"=0))
      sol <-cbind(sol,Cor=res$solution[ts_indices])
      if(!quiet) message(paste("Constraint: max (should be zero):",max(constr_fun(res$solution)),"max at constr:",
                               which.max(constr_fun(res$solution)),"out of ",nrow(constr_sparse_matrix)))

    }

  }



  # Get a good starting point




  if(meffective_testing){

    meffective <- dplyr::filter(dat,testing_indices) %>%  mutate(inds=1:n()) %>% arrange(groups) %>% group_by(groups)  %>%
      summarise(ms=n(),inds=list(inds)) %>% #size and indeces of submatrices
      rowwise() %>%
      mutate(m_eff=calc_m_effecitve(Sigma,inds)) %>%
      select(m_eff) %>% unlist(use.names=F)

    fwer_constr_meffective<- matrix(c(rep(0,nbins),
                                      rep(1,nbins)*meffective,
                                      rep(0,ncol(constr_matrix)-2*nbins)), nrow=1)
    constr_matrix_meffective <- rbind(constr_matrix, fwer_constr_meffective)
    rhs_meffective <- c(rhs,alpha)
    res <- lpsymphony::lpsymphony_solve_LP(obj, constr_matrix_meffective, rep("<=", nrow(constr_matrix_meffective)),
                                           rhs_meffective, #bounds= rsymphony_bounds,
                                           max = TRUE, verbosity = -2, first_feasible = FALSE)

    sol <-cbind(sol,meffective=res$solution[ts_indices])
  }
  if(IHW_testing){

    fwer_constr_IHW<- matrix(c(rep(0,nbins),
                               rep(1,nbins)*m_groups,
                               rep(0,ncol(constr_matrix)-2*nbins)), nrow=1)
    constr_matrix_IHW<- rbind(constr_matrix, fwer_constr_IHW)
    rhs_IHW <- c(rhs,alpha)
    res <- lpsymphony::lpsymphony_solve_LP(obj, constr_matrix_IHW, rep("<=", nrow(constr_matrix_IHW)),
                                           rhs_IHW, #bounds= rsymphony_bounds,
                                           max = TRUE, verbosity = -2, first_feasible = FALSE)
    sol <- cbind(sol,IHW=res$solution[ts_indices])
  }




  sol_df <- tibble(groups=unique(dat$groups)) %>% arrange(groups) %>% cbind(sol)
  sol_name <- colnames(sol)[1] # check number of rejections for cross validation
  rjs <- left_join(dplyr::filter(dat,testing_indices), sol_df,by=c("groups"))
  rjs <- sum(rjs$pvalues<rjs[,sol_name])

  sol_df <- mutate(sol_df,rjs=rjs)
  return(sol_df)
}


