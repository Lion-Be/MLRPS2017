# ----------------------------------------------------------------------------------
# Combining Evidence for Inequality Constrained Hypotheses Using Multiple Studies
# Simulation Study - Defining Parameters
# Lion Behrens
# ----------------------------------------------------------------------------------

# Loading relevant packages (my code)
library(MASS)   # to call mvrnorm
library(Matrix) # to transform vcov matrix to nearest positive definite matrix
library(peperr) # to bootstrap
library(Bain)   # to construct BFs
library(dplyr)  # to shift y to first column in std function

# Loading relevant packages (Gorica code)
library(easypackages)
list.of.packages <- c("base", "blavaan", "boot", "coda", "foreign", "FRACTION", "lavaan", "lme4",
                      "MASS", "matrixcalc", "mvtnorm", "nlme", "quadprog", "R2OpenBUGS")
libraries(list.of.packages)

# Defining relevant functions for Gorica weights analysis
gorica <-
  function(object, ..., iter=100000){
    UseMethod("gorica")
  }


gorica_penalty <-
  function(object, iter=100000, mc.cores=1){
    
    if (!(inherits(object, "ormle") )) stop("object needs to be of class ormle")
    if (all(object$constr == 0) & object$nec == 0){
      penalty <- length(object$est)
    } else {
      if (iter < 1) stop("No of iterations < 1")
      
      est<-object$est 
      K<-length(est)
      covmtrx <- object$covmtrx
      constr<-object$constr
      rhs=object$rhs
      nec=object$nec
      
      Z <- rmvnorm(n=iter, mean=rep(0, K), sigma=covmtrx)
      Dmat2=2*ginv(covmtrx)
      
      nact <- apply(Z, 1, function(z){ 
        
        dvec2=2*(z%*%ginv(covmtrx)) 
        solveQP2= solve.QP(Dmat2,dvec2,t(constr),rhs,meq =nec,factorized = FALSE)
        if (solveQP2$iact[1] == 0) return(0) else return(length(solveQP2$iact))
      })
      
      dimsol <- K - nact
      LP <- sapply(1:K, function(x) sum(x == (dimsol)))/iter
      penalty <- sum((1:K)*LP[])
      
    }
    
    return(penalty)
    
  }



gorica.ormle <-
  function(object, ..., iter=100000){
    if (!inherits(object, "ormle") & !inherits(object, "list")) stop("object needs to be of class ormle or a list of ormle objects")
    if (iter < 1) stop("No of iterations < 1")
    if (inherits(object, "ormle")) objlist <- list(object, ...) else objlist <- object
    isorlm <- sapply(objlist, function(x) inherits(x, "ormle"))
    orlmlist <- objlist[isorlm]  
    Call <- match.call()
    Call$iter <- NULL
    if (inherits(object, "ormle")) names(orlmlist) <- as.character(Call[-1L])[isorlm]
    loglik <- -2*sapply(orlmlist, function(x) x$logLik)
    penalty <- 2*sapply(orlmlist, function(x) gorica_penalty(x, iter=iter))
    gorica <- loglik + penalty
    delta <- gorica - min(gorica)
    gorica_weights <- exp(-delta/2) / sum(exp(-delta/2))
    data.frame(misfit=loglik,complexity=penalty,gorica=gorica,gorica_weights=round(gorica_weights,4))
  }


gorica.list <- function(object, ..., iter=100000){
  if (all(sapply(object, class) == "ormle")) out <- gorica.ormle(object, iter=iter)
  return(out)
}


ormle <-
  function(est,covmtrx,constr,rhs,nec){
    K=length(est)
    covmtrx = as.matrix(covmtrx)
    Dmat = 2*ginv(covmtrx)
    dvec = 2*(est%*% ginv(covmtrx))
    solveQP = solve.QP(Dmat, dvec = dvec, t(constr), rhs, meq = nec, factorized = FALSE)
    tildeQ = solveQP$solution
    restrictedest=solveQP$solution
    names(restrictedest)=names(est)
    loglik =as.numeric( ( -K/2*log(2*pi) )-( 0.5*log(det(covmtrx) ) )-( 0.5* t(est- tildeQ)%*%ginv(covmtrx)%*% (est-tildeQ)) )
    
    out <- list(est=est, covmtrx=covmtrx, constr=constr, rhs=rhs, nec=nec, logLik=loglik,restrictedest=restrictedest)
    class(out) <- "ormle"
    return(out)
    
  }

print.ormle <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\n$est\n")
  print(x$est) 
  cat("\n")
  cat("\n$restrictedest\n")
  print(x$restrictedest) 
  cat("\n")
  invisible(x)
}

# Define standardization function
std <- function(data) {
  
  # Shift y to first column
  data <- data %>%
    select(y, everything())
  
  # Standardize all other except the first/dependent variable
  for (i in 2:length(data)) {
    data[,i] <- scale(data[,i])
  }
  as.data.frame(data)
}


# ---------------------------------- #
# Specifying simulation parameters   #
# ---------------------------------- #

# -------------
# Hypotheses
# -------------

n.hypos <- 2 # simulation code generalized for up to 4 hypotheses


# Bain: Specifying hypotheses

# H1: coef1>coef2>coef3>coef4 
ERr1<-matrix(0,0,0)  
IRr1<-matrix(c(0,1,-1,0,0,0,
               0,0,1,-1,0,0,
               0,0,0,1,-1,0),nrow=3,ncol=6, byrow = TRUE)

# H2: coef2>{coef1, coef3}>coef4 
ERr2<-matrix(0,0,0)  
IRr2<-matrix(c(0,-1,1,0,0,0,
               0,0,1,-1,0,0,
               0,0,0,1,-1,0,
               0,1,0,0,-1,0),nrow=4,ncol=6, byrow = TRUE)

# H3: coef1<coef2<coef3<coef4 
ERr3<-matrix(0,0,0)  
IRr3<-matrix(c(0,-1,1,0,0,0,
               0,0,-1,1,0,0,
               0,0,0,-1,1,0),nrow=3,ncol=6, byrow = TRUE)


# Gorica: Specifying hypotheses

# H1: coef1>coef2>coef3>coef4 
constr1<-matrix(c(0,1,-1,0,0,
                  0,0,1,-1,0,
                  0,0,0,1,-1),nrow=3,ncol=5, byrow = TRUE)
rhs1 <- rep(0, 3) # right hand side values
nec1 <- 0 # number of equality constrains

# H2: coef2>{coef1, coef3}>coef4
constr2<-matrix(c(0,-1,1,0,0,
                  0,0,1,-1,0,
                  0,0,0,1,-1,
                  0,1,0,0,-1),nrow=4,ncol=5, byrow = TRUE)
rhs2 <- rep(0, 4) # right hand side values
nec2 <- 0 # number of equality constrains

# H3: coef1<coef2<coef3<coef4
constr3<-matrix(c(0,-1,1,0,0,
                  0,0,-1,1,0,
                  0,0,0,-1,1),nrow=3,ncol=5, byrow = TRUE)
rhs3 <- rep(0, 3) # right hand side values
nec3 <- 0 # number of equality constrains

# Hu: coef1, coef2, coef3, coef4, coef5 (unconstrained)
constr <- matrix(c(rep(0, 5)), nrow = 1, ncol = 5, byrow = TRUE)
rhs <- rep(0, 1)
nec <- 0



# ---------------------------------------------
# Variance/Covariance Matrices in Populations
# ---------------------------------------------

# Population 1 (columns are X1, X2, X3, X4 without intercept)  
sigma.1 <- matrix(c(1,    0.25,  0.25,  0.25,   
                    0.25,   1,     0.25,  0.25,    
                    0.25,   0.25,  1,     0.25,     
                    0.25,   0.25,  0.25,   1),nrow=4,ncol=4)
mu.1 <- rep(0,4)

# Population 2  
sigma.2 <- matrix(c(1,    0.25,  0.25,  0.25,   
                    0.25,   1,     0.25,  0.25,    
                    0.25,   0.25,  1,     0.25,     
                    0.25,   0.25,  0.25,   1),nrow=4,ncol=4)
mu.2 <- rep(0,4)

# Population 3 
sigma.3 <- matrix(c(1,    0.25,  0.25,  0.25,   
                    0.25,   1,     0.25,  0.25,    
                    0.25,   0.25,  1,     0.25,     
                    0.25,   0.25,  0.25,   1),nrow=4,ncol=4)
mu.3 <- rep(0,4)

# Population 4 
sigma.4 <- matrix(c(1,    0.25,  0.25,  0.25,   
                    0.25,   1,     0.25,  0.25,    
                    0.25,   0.25,  1,     0.25,     
                    0.25,   0.25,  0.25,   1),nrow=4,ncol=4)
mu.4 <- rep(0,4)



# -------------------
# Regression models
# -------------------

model.1 <- "lm(y~V1+V2+V3+V4, data=boot.sample1)"
model.2 <- "lm(y~V1+V2+V3+V4, data=boot.sample2)"
model.3 <- "lm(y~V1+V2+V3+V4, data=boot.sample3)"
model.4 <- "lm(y~V1+V2+V3+V4, data=boot.sample4)"



# -----------------------
# Regression parameters
# -----------------------

# Population 1
n.coef1 <- 4 # without intercept
intercept.1 <- 0
ratios.1 <- c(1.33, 1.77, 1.33, 1) # ratios between true beta coefficients
f2.1 <- 0.15 # r2 = f2/(1+f2) 
ss.1 <- 100 # sample size

# Population 2
n.coef2 <- 4
intercept.2 <- 0
ratios.2 <- c(1.33, 1.77, 1.33, 1) # ratios between true beta coefficients
f2.2 <- 0.15 # r2 = f2/(1+f2) 
ss.2 <- 100 # sample size

# Population 3
n.coef3 <- 4
intercept.3 <- 0
ratios.3 <- c(1.33, 1.77, 1.33, 1) # ratios between true beta coefficients
f2.3 <- 0.15 # r2 = f2/(1+f2) 
ss.3 <- 100 # sample size

# Population 4
n.coef4 <- 4
intercept.4 <- 0
ratios.4 <- c(1.33, 1.77, 1.33, 1) # ratios between true beta coefficients
f2.4 <- 0.15 # r2 = f2/(1+f2) 
ss.4 <- 100 # sample size




# ------------------------------------
# Simulation and Bootstap Iterations
# ------------------------------------

gorica.iterations <- 10000

# Population 1
# Simulation Iterations
n.iter.1 <- 100
# Bootstrap Iterations
n.iter1 <- 100

# Population 2
# Simulation Iterations
n.iter.2 <- 100
# Bootstrap Iterations
n.iter2 <- 100

# Population 3
# Simulation Iterations
n.iter.3 <- 100
# Bootstrap Iterations
n.iter3 <- 100

# Population 4
# Simulation Iterations
n.iter.4 <- 100
# Bootstrap Iterations
n.iter4 <- 100




# ----------------------------------------------------------------------------------
# Combining Evidence for Inequality Constrained Hypotheses Using Multiple Studies
# Simulation Study - Executing Simulation Runs
# Lion Behrens
# ----------------------------------------------------------------------------------

set.seed(301180)

# -------------------------- #
# Population 1 / Dataset 1   #
# -------------------------- #

  # Create empty matrices for storing statistics of interest
  
    # Bayes Factors against H_u
    bfmu.1 <- matrix(NA, nrow=n.iter.1, ncol=n.hypos+1)
    if (n.hypos==1){colnames(bfmu.1) <- c("BF(H1u)", "BF(Huc)")} else if (n.hypos==2){
                    colnames(bfmu.1) <- c("BF(H1u)", "BF(H2u)", "BF(Huc)")} else if (n.hypos==3){
                    colnames(bfmu.1) <- c("BF(H1u)", "BF(H2u)", "BF(H3u)", "BF(Huc)")} else if (n.hypos==4){
                    colnames(bfmu.1) <- c("BF(H1u)", "BF(H2u)", "BF(H3u)", "BF(H4u)", "BF(Huc)")}
    
    # Posterior Model Probabilities
    pmp.1 <- matrix(NA, nrow=n.iter.1, ncol=n.hypos+1)
    if (n.hypos==1){colnames(pmp.1) <- c("P(H1)", "P(Hu)")} else if (n.hypos==2){
                    colnames(pmp.1) <- c("P(H1)", "P(H2)", "P(Hu)")} else if (n.hypos==3){
                    colnames(pmp.1) <- c("P(H1)", "P(H2)", "P(H3)", "P(Hu)")} else if (n.hypos==4){
                    colnames(pmp.1) <- c("P(H1)", "P(H2)", "P(H3)", "P(H4)", "P(Hu)")}
    
    # Gorica values
    gorica.1 <- matrix(NA, nrow=n.iter.1, ncol=n.hypos+1)
    if (n.hypos==1){colnames(gorica.1) <- c("Value(H1)", "Value(Hu)")} else if (n.hypos==2){
                    colnames(gorica.1) <- c("Value(H1)", "Value(H2)", "Value(Hu)")} else if (n.hypos==3){
                    colnames(gorica.1) <- c("Value(H1)", "Value(H2)", "Value(H3)", "Value(Hu)")} else if (n.hypos==4){
                    colnames(gorica.1) <- c("Value(H1)", "Value(H2)", "Value(H3)", "Value(H4)", "Value(Hu)")}
    
    # Gorica weights
    weights.1 <- matrix(NA, nrow=n.iter.1, ncol=n.hypos+1)
    if (n.hypos==1){colnames(weights.1) <- c("Weight(H1)", "Weight(Hu)")} else if (n.hypos==2){
                    colnames(weights.1) <- c("Weight(H1)", "Weight(H2)", "Weight(Hu)")} else if (n.hypos==3){
                    colnames(weights.1) <- c("Weight(H1)", "Weight(H2)", "Weight(H3)", "Weight(Hu)")} else if (n.hypos==4){
                    colnames(weights.1) <- c("Weight(H1)", "Weight(H2)", "Weight(H3)", "Weight(H4)", "Weight(Hu)")}
    
    
    
  # Simulation
  for (i in 1:n.iter.1){
    
    # Create sample of predictor variables
    sample1 <- as.data.frame(mvrnorm(n = ss.1, mu.1, sigma.1, empirical=FALSE))
    
    # Determine x and true beta coefficients
    r2.1 <- f2.1/(1+f2.1)
    fun.1 <- function (x) sum((ratios.1[1]*x)^2,(ratios.1[2]*x)^2,(ratios.1[3]*x)^2,(ratios.1[4]*x)^2) +
      2*(sum((ratios.1[1]*x)*(ratios.1[2]*x)*cor(sample1$V1, sample1$V2), 
             (ratios.1[1]*x)*(ratios.1[3]*x)*cor(sample1$V1, sample1$V3), 
             (ratios.1[1]*x)*(ratios.1[4]*x)*cor(sample1$V1, sample1$V4), 
             (ratios.1[2]*x)*(ratios.1[3]*x)*cor(sample1$V2, sample1$V3), 
             (ratios.1[2]*x)*(ratios.1[4]*x)*cor(sample1$V2, sample1$V4), 
             (ratios.1[3]*x)*(ratios.1[4]*x)*cor(sample1$V3, sample1$V4)))-r2.1 
    
    x.1 <- uniroot(fun.1, lower=0, upper=100)$root
    betas.1 <- ratios.1*x.1
    
    # Add error variable based on desired r2 (because all vars are standardized, var(resid)=1-R2)
    var.epsilon.1 <- 1-r2.1
    sample1$epsilon <- rnorm(ss.1, sd=sqrt(var.epsilon.1))
    
    # Add y variable based on true coefficients and desired r2
    sample1$y <- intercept.1 + betas.1[1]*sample1$V1 + betas.1[2]*sample1$V2 + 
      betas.1[3]*sample1$V3 + betas.1[4]*sample1$V4 + sample1$epsilon
    
    
    # Obtain Variance/Covariance Matrix and Parameter Estimates via Bootstrapping
    
      # Create empty matrix with n.iter as rows and n.coef as columns
      betas1 <- matrix(data=NA, nrow=n.iter1, ncol=n.coef1+1)
      
      # Bootstrap observations, store them in matrix
      for (j in 1:n.iter1){
        
        # Bootstrap sample
        indices1 <- resample.indices(nrow(sample1), n=nrow(sample1), method="boot")
        boot.sample1 <- sample1[indices1$sample.index[[1]],] 
        
        # Standardizing independent variables
        boot.sample1 <- std(boot.sample1)
        
        # Regression Model 
        betas1[j,] <- coef(eval(parse(text=model.1)))
      }
    
      # Obtain Vcov and beta-hat
      ones1       <- rep(1, n.iter1)
      meanvector1 <- (1 / n.iter1) * t(betas1) %*% ones1
      meanmatrix1 <- ones1 %*% t(meanvector1)
      VCovMatrix1  <- (1 / (n.iter1-1)) * t(betas1 - meanmatrix1) %*% (betas1 - meanmatrix1)
      betahat1 <- colMeans(betas1)
    
    # Analysis in Bain
      if (n.hypos==1){
        result1 <- Bain(estimate=betahat1,
                        covariance=VCovMatrix1, 
                        samp=ss.1,
                        nspec=0,
                        njoint=n.coef1+1,
                        ERr1, IRr1)}
      else if (n.hypos==2){
        result1 <- Bain(estimate=betahat1,
                        covariance=VCovMatrix1, 
                        samp=ss.1,
                        nspec=0,
                        njoint=n.coef1+1,
                        ERr1, IRr1, ERr2, IRr2)}
      else if (n.hypos==3){
        result1 <- Bain(estimate=betahat1,
                          covariance=VCovMatrix1, 
                          samp=ss.1,
                          nspec=0,
                          njoint=n.coef1+1,
                          ERr1, IRr1, ERr2, IRr2, ERr3, IRr3)}
      else if (n.hypos==4){
        result1 <- Bain(estimate=betahat1,
                          covariance=VCovMatrix1, 
                          samp=ss.1,
                          nspec=0,
                          njoint=n.coef1+1,
                          ERr1, IRr1, ERr2, IRr2, ERr3, IRr3, ERr4, IRr4)}
        
      # Store statistics of interest
      if (n.hypos==1){    
        bfmu.1[i,1] <- result1$testResult$fit[1]/result1$testResult$complexity[1]
        bfmu.1[i,2] <- 1} 
      else if (n.hypos==2){
        bfmu.1[i,1] <- result1$testResult$fit[1]/result1$testResult$complexity[1]
        bfmu.1[i,2] <- result1$testResult$fit[2]/result1$testResult$complexity[2]
        bfmu.1[i,3] <- 1}
      else if (n.hypos==3){
        bfmu.1[i,1] <- result1$testResult$fit[1]/result1$testResult$complexity[1]
        bfmu.1[i,2] <- result1$testResult$fit[2]/result1$testResult$complexity[2]
        bfmu.1[i,3] <- result1$testResult$fit[3]/result1$testResult$complexity[3]
        bfmu.1[i,4] <- 1}
      else if (n.hypos==4){
        bfmu.1[i,1] <- result1$testResult$fit[1]/result1$testResult$complexity[1]
        bfmu.1[i,2] <- result1$testResult$fit[2]/result1$testResult$complexity[2]
        bfmu.1[i,3] <- result1$testResult$fit[3]/result1$testResult$complexity[3]
        bfmu.1[i,4] <- result1$testResult$fit[4]/result1$testResult$complexity[4]
        bfmu.1[i,5] <- 1}
      
      pmp.1[i,] <- c(result1$testResult$PMPb, 1-sum(result1$testResult$PMPb))
    
    
    # Gorica Analysis
      if (n.hypos==1){
        # Obtain order-restricted estimates
        H1 <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr1,
                    nec = nec1, rhs = rhs1)
        Hu <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr,
                    nec = nec, rhs = rhs)
        # Construct Gorica values and weights
        gorica.output1 <- gorica(H1, Hu, iter = gorica.iterations)
        gorica.output1}
      else if (n.hypos==2){
        # Obtain order-restricted estimates
        H1 <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr1,
                    nec = nec1, rhs = rhs1)
        H2 <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr2,
                    nec = nec2, rhs = rhs2)
        Hu <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr,
                    nec = nec, rhs = rhs)
        # Construct Gorica values and weights
        gorica.output1 <- gorica(H1, H2, Hu, iter = gorica.iterations)
        gorica.output1}
      else if (n.hypos==3){
        # Obtain order-restricted estimates
        H1 <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr1,
                    nec = nec1, rhs = rhs1)
        H2 <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr2,
                    nec = nec2, rhs = rhs2)
        H3 <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr3,
                    nec = nec3, rhs = rhs3)
        Hu <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr,
                    nec = nec, rhs = rhs)
        # Construct Gorica values and weights
        gorica.output1 <- gorica(H1, H2, H3, Hu, iter = gorica.iterations)
        gorica.output1}
      else if (n.hypos==4){
        # Obtain order-restricted estimates
        H1 <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr1,
                    nec = nec1, rhs = rhs1)
        H2 <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr2,
                    nec = nec2, rhs = rhs2)
        H3 <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr3,
                    nec = nec3, rhs = rhs3)
        H4 <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr4,
                    nec = nec4, rhs = rhs4)
        Hu <- ormle(est = betahat1, covmtrx = VCovMatrix1, const = constr,
                    nec = nec, rhs = rhs)
        # Construct Gorica values and weights
        gorica.output1 <- gorica(H1, H2, H3, H4, Hu, iter = gorica.iterations)
        gorica.output1} 
      
        
      # Store statistics of interest 
      gorica.1[i,] <- gorica.output1$gorica
      weights.1[i,] <- gorica.output1$gorica_weights
      
    
    # Print out simulation process
    cat(paste0("Drawing from Population 1, "))
    cat(paste0("Simulation iteration:",i, "/", n.iter.1 ))
   
  }
  
    
    
    
# -------------------------- #
# Population 2 / Dataset 2   #
# -------------------------- #
    
  # Create empty matrices for storing statistics of interest
  
    # Bayes Factors against H_u
    bfmu.2 <- matrix(NA, nrow=n.iter.2, ncol=n.hypos+1)
    if (n.hypos==1){colnames(bfmu.2) <- c("BF(H1u)", "BF(Huc)")} else if (n.hypos==2){
                    colnames(bfmu.2) <- c("BF(H1u)", "BF(H2u)", "BF(Huc)")} else if (n.hypos==3){
                    colnames(bfmu.2) <- c("BF(H1u)", "BF(H2u)", "BF(H3u)", "BF(Huc)")} else if (n.hypos==4){
                    colnames(bfmu.2) <- c("BF(H1u)", "BF(H2u)", "BF(H3u)", "BF(H4u)", "BF(Huc)")}
    
    # Posterior Model Probabilities
    pmp.2 <- matrix(NA, nrow=n.iter.2, ncol=n.hypos+1)
    if (n.hypos==1){colnames(pmp.2) <- c("P(H1)", "P(Hu)")} else if (n.hypos==2){
                    colnames(pmp.2) <- c("P(H1)", "P(H2)", "P(Hu)")} else if (n.hypos==3){
                    colnames(pmp.2) <- c("P(H1)", "P(H2)", "P(H3)", "P(Hu)")} else if (n.hypos==4){
                    colnames(pmp.2) <- c("P(H1)", "P(H2)", "P(H3)", "P(H4)", "P(Hu)")}
  
    # Gorica values
    gorica.2 <- matrix(NA, nrow=n.iter.2, ncol=n.hypos+1)
    if (n.hypos==1){colnames(gorica.2) <- c("Value(H1)", "Value(Hu)")} else if (n.hypos==2){
                    colnames(gorica.2) <- c("Value(H1)", "Value(H2)", "Value(Hu)")} else if (n.hypos==3){
                    colnames(gorica.2) <- c("Value(H1)", "Value(H2)", "Value(H3)", "Value(Hu)")} else if (n.hypos==4){
                    colnames(gorica.2) <- c("Value(H1)", "Value(H2)", "Value(H3)", "Value(H4)", "Value(Hu)")}
    
    # Gorica weights
    weights.2 <- matrix(NA, nrow=n.iter.2, ncol=n.hypos+1)
    if (n.hypos==1){colnames(weights.2) <- c("Weight(H1)", "Weight(Hu)")} else if (n.hypos==2){
                    colnames(weights.2) <- c("Weight(H1)", "Weight(H2)", "Weight(Hu)")} else if (n.hypos==3){
                    colnames(weights.2) <- c("Weight(H1)", "Weight(H2)", "Weight(H3)", "Weight(Hu)")} else if (n.hypos==4){
                    colnames(weights.2) <- c("Weight(H1)", "Weight(H2)", "Weight(H3)", "Weight(H4)", "Weight(Hu)")}
    
  # Simulation
  for (i in 1:n.iter.2){
    
    # Create sample of predictor variables
    sample2 <- as.data.frame(mvrnorm(n = ss.2, mu.2, sigma.2, empirical=FALSE))
    
    # Determine x and true beta coefficients
    r2.2 <- f2.2/(1+f2.2)
    fun.2 <- function (x) sum((ratios.2[1]*x)^2,(ratios.2[2]*x)^2,(ratios.2[3]*x)^2,(ratios.2[4]*x)^2) +
      2*(sum((ratios.2[1]*x)*(ratios.2[2]*x)*cor(sample2$V1, sample2$V2), 
             (ratios.2[1]*x)*(ratios.2[3]*x)*cor(sample2$V1, sample2$V3), 
             (ratios.2[1]*x)*(ratios.2[4]*x)*cor(sample2$V1, sample2$V4), 
             (ratios.2[2]*x)*(ratios.2[3]*x)*cor(sample2$V2, sample2$V3), 
             (ratios.2[2]*x)*(ratios.2[4]*x)*cor(sample2$V2, sample2$V4), 
             (ratios.2[3]*x)*(ratios.2[4]*x)*cor(sample2$V3, sample2$V4)))-r2.2 
    
    x.2 <- uniroot(fun.2, lower=0, upper=100)$root
    betas.2 <- ratios.2*x.2
    
    # Add error variable based on desired r2 (because all vars are standardized, var(resid)=1-R2)
    var.epsilon.2 <- 1-r2.2
    sample2$epsilon <- rnorm(ss.2, sd=sqrt(var.epsilon.2))
    
    # Add y variable based on true coefficients and desired r2
    sample2$y <- intercept.2 + betas.2[1]*sample2$V1 + betas.2[2]*sample2$V2 + 
      betas.2[3]*sample2$V3 + betas.2[4]*sample2$V4 + sample2$epsilon
    
    # Obtain Variance/Covariance Matrix and Parameter Estimates via Bootstrapping
    
      # Create empty matrix with n.iter as rows and n.coef as columns
      betas2 <- matrix(data=NA, nrow=n.iter2, ncol=n.coef2+1)
      
      # Bootstrap observations, store them in matrix
      for (j in 1:n.iter2){
        
        # Bootstrap sample
        indices2 <- resample.indices(nrow(sample2), n=nrow(sample2), method="boot")
        boot.sample2 <- sample2[indices2$sample.index[[1]],] 
        
        # Standardizing independent variables
        boot.sample2 <- std(boot.sample2)
        
        # Regression Model 
        betas2[j,] <- coef(eval(parse(text=model.2)))
      }
      
      # Obtain Vcov and beta-hat
      ones2       <- rep(1, n.iter2)
      meanvector2 <- (1 / n.iter2) * t(betas2) %*% ones2
      meanmatrix2 <- ones2 %*% t(meanvector2)
      VCovMatrix2  <- (1 / (n.iter2-1)) * t(betas2 - meanmatrix2) %*% (betas2 - meanmatrix2)
      betahat2 <- colMeans(betas2)
    
    # Analysis in Bain
      if (n.hypos==1){
        result2 <- Bain(estimate=betahat2,
                        covariance=VCovMatrix2, 
                        samp=ss.2,
                        nspec=0,
                        njoint=n.coef2+1,
                        ERr1, IRr1)}
      else if (n.hypos==2){
        result2 <- Bain(estimate=betahat2,
                        covariance=VCovMatrix2, 
                        samp=ss.2,
                        nspec=0,
                        njoint=n.coef2+1,
                        ERr1, IRr1, ERr2, IRr2)}
      else if (n.hypos==3){
        result2 <- Bain(estimate=betahat2,
                        covariance=VCovMatrix2, 
                        samp=ss.2,
                        nspec=0,
                        njoint=n.coef2+1,
                        ERr1, IRr1, ERr2, IRr2, ERr3, IRr3)}
      else if (n.hypos==4){
        result2 <- Bain(estimate=betahat2,
                        covariance=VCovMatrix2, 
                        samp=ss.2,
                        nspec=0,
                        njoint=n.coef2+1,
                        ERr1, IRr1, ERr2, IRr2, ERr3, IRr3, ERr4, IRr4)}
      
      # Store statistics of interest
      if (n.hypos==1){    
        bfmu.2[i,1] <- result2$testResult$fit[1]/result2$testResult$complexity[1]
        bfmu.2[i,2] <- 1} 
      else if (n.hypos==2){
        bfmu.2[i,1] <- result2$testResult$fit[1]/result2$testResult$complexity[1]
        bfmu.2[i,2] <- result2$testResult$fit[2]/result2$testResult$complexity[2]
        bfmu.2[i,3] <- 1}
      else if (n.hypos==3){
        bfmu.2[i,1] <- result2$testResult$fit[1]/result2$testResult$complexity[1]
        bfmu.2[i,2] <- result2$testResult$fit[2]/result2$testResult$complexity[2]
        bfmu.2[i,3] <- result2$testResult$fit[3]/result2$testResult$complexity[3]
        bfmu.2[i,4] <- 1}
      else if (n.hypos==4){
        bfmu.2[i,1] <- result2$testResult$fit[1]/result2$testResult$complexity[1]
        bfmu.2[i,2] <- result2$testResult$fit[2]/result2$testResult$complexity[2]
        bfmu.2[i,3] <- result2$testResult$fit[3]/result2$testResult$complexity[3]
        bfmu.2[i,4] <- result2$testResult$fit[4]/result2$testResult$complexity[4]
        bfmu.2[i,5] <- 1}
      
      pmp.2[i,] <- c(result2$testResult$PMPb, 1-sum(result2$testResult$PMPb))
    
    
    # Gorica Analysis
    if (n.hypos==1){
      # Obtain order-restricted estimates
      H1 <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr1,
                  nec = nec1, rhs = rhs1)
      Hu <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr,
                  nec = nec, rhs = rhs)
      # Construct Gorica values and weights
      gorica.output2 <- gorica(H1, Hu, iter = gorica.iterations)
      gorica.output2}
    else if (n.hypos==2){
      # Obtain order-restricted estimates
      H1 <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr1,
                  nec = nec1, rhs = rhs1)
      H2 <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr2,
                  nec = nec2, rhs = rhs2)
      Hu <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr,
                  nec = nec, rhs = rhs)
      # Construct Gorica values and weights
      gorica.output2 <- gorica(H1, H2, Hu, iter = gorica.iterations)
      gorica.output2}
    else if (n.hypos==3){
      # Obtain order-restricted estimates
      H1 <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr1,
                  nec = nec1, rhs = rhs1)
      H2 <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr2,
                  nec = nec2, rhs = rhs2)
      H3 <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr3,
                  nec = nec3, rhs = rhs3)
      Hu <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr,
                  nec = nec, rhs = rhs)
      # Construct Gorica values and weights
      gorica.output2 <- gorica(H1, H2, H3, Hu, iter = gorica.iterations)
      gorica.output2}
    else if (n.hypos==4){
      # Obtain order-restricted estimates
      H1 <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr1,
                  nec = nec1, rhs = rhs1)
      H2 <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr2,
                  nec = nec2, rhs = rhs2)
      H3 <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr3,
                  nec = nec3, rhs = rhs3)
      H4 <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr4,
                  nec = nec4, rhs = rhs4)
      Hu <- ormle(est = betahat2, covmtrx = VCovMatrix2, const = constr,
                  nec = nec, rhs = rhs)
      # Construct Gorica values and weights
      gorica.output2 <- gorica(H1, H2, H3, H4, Hu, iter = gorica.iterations)
      gorica.output2} 
    
    
    # Store statistics of interest 
    gorica.2[i,] <- gorica.output2$gorica
    weights.2[i,] <- gorica.output2$gorica_weights
    
    
    
    # Print out simulation process
    cat(paste0("Drawing from Population 2, "))
    cat(paste0("Simulation iteration:",i, "/", n.iter.2 ))
    
  }
  

    
    
    
# -------------------------- #
# Population 3 / Dataset 3   #
# -------------------------- #
  
  # Create empty matrices for storing statistics of interest
  
    # Bayes Factors against H_u
    bfmu.3 <- matrix(NA, nrow=n.iter.3, ncol=n.hypos+1)
    if (n.hypos==1){colnames(bfmu.3) <- c("BF(H1u)", "BF(Huc)")} else if (n.hypos==2){
                    colnames(bfmu.3) <- c("BF(H1u)", "BF(H2u)", "BF(Huc)")} else if (n.hypos==3){
                    colnames(bfmu.3) <- c("BF(H1u)", "BF(H2u)", "BF(H3u)", "BF(Huc)")} else if (n.hypos==4){
                    colnames(bfmu.3) <- c("BF(H1u)", "BF(H2u)", "BF(H3u)", "BF(H4u)", "BF(Huc)")}
    
    # Posterior Model Probabilities
    pmp.3 <- matrix(NA, nrow=n.iter.3, ncol=n.hypos+1)
    if (n.hypos==1){colnames(pmp.3) <- c("P(H1)", "P(Hu)")} else if (n.hypos==2){
                    colnames(pmp.3) <- c("P(H1)", "P(H2)", "P(Hu)")} else if (n.hypos==3){
                    colnames(pmp.3) <- c("P(H1)", "P(H2)", "P(H3)", "P(Hu)")} else if (n.hypos==4){
                    colnames(pmp.3) <- c("P(H1)", "P(H2)", "P(H3)", "P(H4)", "P(Hu)")}
    
    # Gorica values
    gorica.3 <- matrix(NA, nrow=n.iter.3, ncol=n.hypos+1)
    if (n.hypos==1){colnames(gorica.3) <- c("Value(H1)", "Value(Hu)")} else if (n.hypos==2){
                    colnames(gorica.3) <- c("Value(H1)", "Value(H2)", "Value(Hu)")} else if (n.hypos==3){
                    colnames(gorica.3) <- c("Value(H1)", "Value(H2)", "Value(H3)", "Value(Hu)")} else if (n.hypos==4){
                    colnames(gorica.3) <- c("Value(H1)", "Value(H2)", "Value(H3)", "Value(H4)", "Value(Hu)")}
    
    # Gorica weights
    weights.3 <- matrix(NA, nrow=n.iter.3, ncol=n.hypos+1)
    if (n.hypos==1){colnames(weights.3) <- c("Weight(H1)", "Weight(Hu)")} else if (n.hypos==2){
                    colnames(weights.3) <- c("Weight(H1)", "Weight(H2)", "Weight(Hu)")} else if (n.hypos==3){
                    colnames(weights.3) <- c("Weight(H1)", "Weight(H2)", "Weight(H3)", "Weight(Hu)")} else if (n.hypos==4){
                    colnames(weights.3) <- c("Weight(H1)", "Weight(H2)", "Weight(H3)", "Weight(H4)", "Weight(Hu)")}
  
  # Simulation
  for (i in 1:n.iter.3){
    
    # Create sample of predictor variables
    sample3 <- as.data.frame(mvrnorm(n = ss.3, mu.3, sigma.3, empirical=FALSE))
    
    # Determine x and true beta coefficients
    r2.3 <- f2.3/(1+f2.3)
    fun.3 <- function (x) sum((ratios.3[1]*x)^2,(ratios.3[2]*x)^2,(ratios.3[3]*x)^2,(ratios.3[4]*x)^2) +
      2*(sum((ratios.3[1]*x)*(ratios.3[2]*x)*cor(sample3$V1, sample3$V2), 
             (ratios.3[1]*x)*(ratios.3[3]*x)*cor(sample3$V1, sample3$V3), 
             (ratios.3[1]*x)*(ratios.3[4]*x)*cor(sample3$V1, sample3$V4), 
             (ratios.3[2]*x)*(ratios.3[3]*x)*cor(sample3$V2, sample3$V3), 
             (ratios.3[2]*x)*(ratios.3[4]*x)*cor(sample3$V2, sample3$V4), 
             (ratios.3[3]*x)*(ratios.3[4]*x)*cor(sample3$V3, sample3$V4)))-r2.3 
    
    x.3 <- uniroot(fun.3, lower=0, upper=100)$root
    betas.3 <- ratios.3*x.3
    
    # Add error variable based on desired r2 (because all vars are standardized, var(resid)=1-R2)
    var.epsilon.3 <- 1-r2.3
    sample3$epsilon <- rnorm(ss.3, sd=sqrt(var.epsilon.3))
    
    # Add y variable based on true coefficients and desired r2
    sample3$y <- intercept.3 + betas.3[1]*sample3$V1 + betas.3[2]*sample3$V2 + 
      betas.3[3]*sample3$V3 + betas.3[4]*sample3$V4 + sample3$epsilon
    
    # Obtain Variance/Covariance Matrix and Parameter Estimates via Bootstrapping
    
      # Create empty matrix with n.iter as rows and n.coef as columns
      betas3 <- matrix(data=NA, nrow=n.iter3, ncol=n.coef3+1)
      
      # Bootstrap observations, store them in matrix
      for (j in 1:n.iter3){
        
        # Bootstrap sample
        indices3 <- resample.indices(nrow(sample3), n=nrow(sample3), method="boot")
        boot.sample3 <- sample3[indices3$sample.index[[1]],] 
        
        # Standardizing independent variables
        boot.sample3 <- std(boot.sample3)
        
        # Regression Model 
        betas3[j,] <- coef(eval(parse(text=model.3)))
      }
      
      # Obtain Vcov and beta-hat
      ones3       <- rep(1, n.iter3)
      meanvector3 <- (1 / n.iter3) * t(betas3) %*% ones3
      meanmatrix3 <- ones3 %*% t(meanvector3)
      VCovMatrix3  <- (1 / (n.iter3-1)) * t(betas3 - meanmatrix3) %*% (betas3 - meanmatrix3)
      betahat3 <- colMeans(betas3)
    
    # Analysis in Bain
      if (n.hypos==1){
        result3 <- Bain(estimate=betahat3,
                        covariance=VCovMatrix3, 
                        samp=ss.3,
                        nspec=0,
                        njoint=n.coef3+1,
                        ERr1, IRr1)}
      else if (n.hypos==2){
        result3 <- Bain(estimate=betahat3,
                        covariance=VCovMatrix3, 
                        samp=ss.3,
                        nspec=0,
                        njoint=n.coef3+1,
                        ERr1, IRr1, ERr2, IRr2)}
      else if (n.hypos==3){
        result3 <- Bain(estimate=betahat3,
                        covariance=VCovMatrix3, 
                        samp=ss.3,
                        nspec=0,
                        njoint=n.coef3+1,
                        ERr1, IRr1, ERr2, IRr2, ERr3, IRr3)}
      else if (n.hypos==4){
        result3 <- Bain(estimate=betahat3,
                        covariance=VCovMatrix3, 
                        samp=ss.3,
                        nspec=0,
                        njoint=n.coef3+1,
                        ERr1, IRr1, ERr2, IRr2, ERr3, IRr3, ERr4, IRr4)}
      
      
      
      # Store statistics of interest
      if (n.hypos==1){    
        bfmu.3[i,1] <- result3$testResult$fit[1]/result3$testResult$complexity[1]
        bfmu.3[i,2] <- 1} 
      else if (n.hypos==2){
        bfmu.3[i,1] <- result3$testResult$fit[1]/result3$testResult$complexity[1]
        bfmu.3[i,2] <- result3$testResult$fit[2]/result3$testResult$complexity[2]
        bfmu.3[i,3] <- 1}
      else if (n.hypos==3){
        bfmu.3[i,1] <- result3$testResult$fit[1]/result3$testResult$complexity[1]
        bfmu.3[i,2] <- result3$testResult$fit[2]/result3$testResult$complexity[2]
        bfmu.3[i,3] <- result3$testResult$fit[3]/result3$testResult$complexity[3]
        bfmu.3[i,4] <- 1}
      else if (n.hypos==4){
        bfmu.3[i,1] <- result3$testResult$fit[1]/result3$testResult$complexity[1]
        bfmu.3[i,2] <- result3$testResult$fit[2]/result3$testResult$complexity[2]
        bfmu.3[i,3] <- result3$testResult$fit[3]/result3$testResult$complexity[3]
        bfmu.3[i,4] <- result3$testResult$fit[4]/result3$testResult$complexity[4]
        bfmu.3[i,5] <- 1}
      
      pmp.3[i,] <- c(result3$testResult$PMPb, 1-sum(result3$testResult$PMPb))
    
    
    # Gorica Analysis
    if (n.hypos==1){
      # Obtain order-restricted estimates
      H1 <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr1,
                  nec = nec1, rhs = rhs1)
      Hu <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr,
                  nec = nec, rhs = rhs)
      # Construct Gorica values and weights
      gorica.output3 <- gorica(H1, Hu, iter = gorica.iterations)
      gorica.output3}
    else if (n.hypos==2){
      # Obtain order-restricted estimates
      H1 <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr1,
                  nec = nec1, rhs = rhs1)
      H2 <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr2,
                  nec = nec2, rhs = rhs2)
      Hu <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr,
                  nec = nec, rhs = rhs)
      # Construct Gorica values and weights
      gorica.output3 <- gorica(H1, H2, Hu, iter = gorica.iterations)
      gorica.output3}
    else if (n.hypos==3){
      # Obtain order-restricted estimates
      H1 <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr1,
                  nec = nec1, rhs = rhs1)
      H2 <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr2,
                  nec = nec2, rhs = rhs2)
      H3 <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr3,
                  nec = nec3, rhs = rhs3)
      Hu <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr,
                  nec = nec, rhs = rhs)
      # Construct Gorica values and weights
      gorica.output3 <- gorica(H1, H2, H3, Hu, iter = gorica.iterations)
      gorica.output3}
    else if (n.hypos==4){
      # Obtain order-restricted estimates
      H1 <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr1,
                  nec = nec1, rhs = rhs1)
      H2 <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr2,
                  nec = nec2, rhs = rhs2)
      H3 <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr3,
                  nec = nec3, rhs = rhs3)
      H4 <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr4,
                  nec = nec4, rhs = rhs4)
      Hu <- ormle(est = betahat3, covmtrx = VCovMatrix3, const = constr,
                  nec = nec, rhs = rhs)
      # Construct Gorica values and weights
      gorica.output3 <- gorica(H1, H2, H3, H4, Hu, iter = gorica.iterations)
      gorica.output3} 
    
    
    # Store statistics of interest 
    gorica.3[i,] <- gorica.output3$gorica
    weights.3[i,] <- gorica.output3$gorica_weights
    
    
    
    # Print out simulation process
    cat(paste0("Drawing from Population 3, "))
    cat(paste0("Simulation iteration:",i, "/", n.iter.3 ))
    
  }

    
    

# -------------------------- #
# Population 4 / Dataset 4   #
# -------------------------- #

# Create empty matrices for storing statistics of interest

# Bayes Factors against H_u
bfmu.4 <- matrix(NA, nrow=n.iter.4, ncol=n.hypos+1)
    if (n.hypos==1){colnames(bfmu.4) <- c("BF(H1u)", "BF(Huc)")} else if (n.hypos==2){
                    colnames(bfmu.4) <- c("BF(H1u)", "BF(H2u)", "BF(Huc)")} else if (n.hypos==3){
                    colnames(bfmu.4) <- c("BF(H1u)", "BF(H2u)", "BF(H3u)", "BF(Huc)")} else if (n.hypos==4){
                    colnames(bfmu.4) <- c("BF(H1u)", "BF(H2u)", "BF(H3u)", "BF(H4u)", "BF(Huc)")}

# Posterior Model Probabilities
pmp.4 <- matrix(NA, nrow=n.iter.4, ncol=n.hypos+1)
if (n.hypos==1){colnames(pmp.4) <- c("P(H1)", "P(Hu)")} else if (n.hypos==2){
                colnames(pmp.4) <- c("P(H1)", "P(H2)", "P(Hu)")} else if (n.hypos==3){
                colnames(pmp.4) <- c("P(H1)", "P(H2)", "P(H3)", "P(Hu)")} else if (n.hypos==4){
                colnames(pmp.4) <- c("P(H1)", "P(H2)", "P(H3)", "P(H4)", "P(Hu)")}

# Gorica values
gorica.4 <- matrix(NA, nrow=n.iter.4, ncol=n.hypos+1)
if (n.hypos==1){colnames(gorica.4) <- c("Value(H1)", "Value(Hu)")} else if (n.hypos==2){
                colnames(gorica.4) <- c("Value(H1)", "Value(H2)", "Value(Hu)")} else if (n.hypos==3){
                colnames(gorica.4) <- c("Value(H1)", "Value(H2)", "Value(H3)", "Value(Hu)")} else if (n.hypos==4){
                colnames(gorica.4) <- c("Value(H1)", "Value(H2)", "Value(H3)", "Value(H4)", "Value(Hu)")}

# Gorica weights
weights.4 <- matrix(NA, nrow=n.iter.4, ncol=n.hypos+1)
if (n.hypos==1){colnames(weights.4) <- c("Weight(H1)", "Weight(Hu)")} else if (n.hypos==2){
                colnames(weights.4) <- c("Weight(H1)", "Weight(H2)", "Weight(Hu)")} else if (n.hypos==3){
                colnames(weights.4) <- c("Weight(H1)", "Weight(H2)", "Weight(H3)", "Weight(Hu)")} else if (n.hypos==4){
                colnames(weights.4) <- c("Weight(H1)", "Weight(H2)", "Weight(H3)", "Weight(H4)", "Weight(Hu)")}


# Simulation
for (i in 1:n.iter.4){
  
  # Create sample of predictor variables
  sample4 <- as.data.frame(mvrnorm(n = ss.4, mu.4, sigma.4, empirical=FALSE))
  
  # Determine x and true beta coefficients
  r2.4 <- f2.4/(1+f2.4)
  fun.4 <- function (x) sum((ratios.4[1]*x)^2,(ratios.4[2]*x)^2,(ratios.4[3]*x)^2,(ratios.4[4]*x)^2) +
    2*(sum((ratios.4[1]*x)*(ratios.4[2]*x)*cor(sample4$V1, sample4$V2), 
           (ratios.4[1]*x)*(ratios.4[3]*x)*cor(sample4$V1, sample4$V3), 
           (ratios.4[1]*x)*(ratios.4[4]*x)*cor(sample4$V1, sample4$V4), 
           (ratios.4[2]*x)*(ratios.4[3]*x)*cor(sample4$V2, sample4$V3), 
           (ratios.4[2]*x)*(ratios.4[4]*x)*cor(sample4$V2, sample4$V4), 
           (ratios.4[3]*x)*(ratios.4[4]*x)*cor(sample4$V3, sample4$V4)))-r2.4 
  
  x.4 <- uniroot(fun.4, lower=0, upper=100)$root
  betas.4 <- ratios.4*x.4
  
  # Add error variable based on desired r2 (because all vars are standardized, var(resid)=1-R2)
  var.epsilon.4 <- 1-r2.4
  sample4$epsilon <- rnorm(ss.4, sd=sqrt(var.epsilon.4))
  
  # Add y variable based on true coefficients and desired r2
  sample4$y <- intercept.4 + betas.4[1]*sample4$V1 + betas.4[2]*sample4$V2 + 
    betas.4[3]*sample4$V3 + betas.4[4]*sample4$V4 + sample4$epsilon
  
  # Obtain Variance/Covariance Matrix and Parameter Estimates via Bootstrapping
  
    # Create empty matrix with n.iter as rows and n.coef as columns
    betas4 <- matrix(data=NA, nrow=n.iter4, ncol=n.coef4+1)
    
    # Bootstrap observations, store them in matrix
    for (j in 1:n.iter4){
      
      # Bootstrap sample
      indices4 <- resample.indices(nrow(sample4), n=nrow(sample4), method="boot")
      boot.sample4 <- sample4[indices4$sample.index[[1]],] 
      
      # Standardizing independent variables
      boot.sample4 <- std(boot.sample4)
      
      # Regression Model 
      betas4[j,] <- coef(eval(parse(text=model.4)))
    }
    
    # Obtain Vcov and beta-hat
    ones4       <- rep(1, n.iter4)
    meanvector4 <- (1 / n.iter4) * t(betas4) %*% ones4
    meanmatrix4 <- ones4 %*% t(meanvector4)
    VCovMatrix4  <- (1 / (n.iter4-1)) * t(betas4 - meanmatrix4) %*% (betas4 - meanmatrix4)
    betahat4 <- colMeans(betas4)
  
  # Analysis in Bain
    if (n.hypos==1){
      result4 <- Bain(estimate=betahat4,
                      covariance=VCovMatrix4, 
                      samp=ss.4,
                      nspec=0,
                      njoint=n.coef4+1,
                      ERr1, IRr1)}
    else if (n.hypos==2){
      result4 <- Bain(estimate=betahat4,
                      covariance=VCovMatrix4, 
                      samp=ss.4,
                      nspec=0,
                      njoint=n.coef4+1,
                      ERr1, IRr1, ERr2, IRr2)}
    else if (n.hypos==3){
      result4 <- Bain(estimate=betahat4,
                      covariance=VCovMatrix4, 
                      samp=ss.4,
                      nspec=0,
                      njoint=n.coef4+1,
                      ERr1, IRr1, ERr2, IRr2, ERr3, IRr3)}
    else if (n.hypos==4){
      result4 <- Bain(estimate=betahat4,
                      covariance=VCovMatrix4, 
                      samp=ss.4,
                      nspec=0,
                      njoint=n.coef4+1,
                      ERr1, IRr1, ERr2, IRr2, ERr3, IRr3, ERr4, IRr4)}
    
    
    
    # Store statistics of interest
    if (n.hypos==1){    
      bfmu.4[i,1] <- result4$testResult$fit[1]/result4$testResult$complexity[1]
      bfmu.4[i,2] <- 1} 
    else if (n.hypos==2){
      bfmu.4[i,1] <- result4$testResult$fit[1]/result4$testResult$complexity[1]
      bfmu.4[i,2] <- result4$testResult$fit[2]/result4$testResult$complexity[2]
      bfmu.4[i,3] <- 1}
    else if (n.hypos==3){
      bfmu.4[i,1] <- result4$testResult$fit[1]/result4$testResult$complexity[1]
      bfmu.4[i,2] <- result4$testResult$fit[2]/result4$testResult$complexity[2]
      bfmu.4[i,3] <- result4$testResult$fit[3]/result4$testResult$complexity[3]
      bfmu.4[i,4] <- 1}
    else if (n.hypos==4){
      bfmu.4[i,1] <- result4$testResult$fit[1]/result4$testResult$complexity[1]
      bfmu.4[i,2] <- result4$testResult$fit[2]/result4$testResult$complexity[2]
      bfmu.4[i,3] <- result4$testResult$fit[3]/result4$testResult$complexity[3]
      bfmu.4[i,4] <- result4$testResult$fit[4]/result4$testResult$complexity[4]
      bfmu.4[i,5] <- 1}
    
    pmp.4[i,] <- c(result4$testResult$PMPb, 1-sum(result4$testResult$PMPb))
  
  
  # Gorica Analysis
  if (n.hypos==1){
    # Obtain order-restricted estimates
    H1 <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr1,
                nec = nec1, rhs = rhs1)
    Hu <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr,
                nec = nec, rhs = rhs)
    # Construct Gorica values and weights
    gorica.output4 <- gorica(H1, Hu, iter = gorica.iterations)
    gorica.output4}
  else if (n.hypos==2){
    # Obtain order-restricted estimates
    H1 <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr1,
                nec = nec1, rhs = rhs1)
    H2 <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr2,
                nec = nec2, rhs = rhs2)
    Hu <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr,
                nec = nec, rhs = rhs)
    # Construct Gorica values and weights
    gorica.output4 <- gorica(H1, H2, Hu, iter = gorica.iterations)
    gorica.output4}
  else if (n.hypos==3){
    # Obtain order-restricted estimates
    H1 <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr1,
                nec = nec1, rhs = rhs1)
    H2 <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr2,
                nec = nec2, rhs = rhs2)
    H3 <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr3,
                nec = nec3, rhs = rhs3)
    Hu <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr,
                nec = nec, rhs = rhs)
    # Construct Gorica values and weights
    gorica.output4 <- gorica(H1, H2, H3, Hu, iter = gorica.iterations)
    gorica.output4}
  else if (n.hypos==4){
    # Obtain order-restricted estimates
    H1 <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr1,
                nec = nec1, rhs = rhs1)
    H2 <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr2,
                nec = nec2, rhs = rhs2)
    H3 <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr3,
                nec = nec3, rhs = rhs3)
    H4 <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr4,
                nec = nec4, rhs = rhs4)
    Hu <- ormle(est = betahat4, covmtrx = VCovMatrix4, const = constr,
                nec = nec, rhs = rhs)
    # Construct Gorica values and weights
    gorica.output4 <- gorica(H1, H2, H3, H4, Hu, iter = gorica.iterations)
    gorica.output4} 
  
  
  # Store statistics of interest 
  gorica.4[i,] <- gorica.output4$gorica
  weights.4[i,] <- gorica.output4$gorica_weights
  
  # Print out simulation process
  cat(paste0("Drawing from Population 4, "))
  cat(paste0("Simulation iteration:",i, "/", n.iter.1 ))
  
}

    
    
    
    
    
    
# ---------------------------------------------------------------------- #
# Running the Aggregation Method Using Bayes Factors and Gorica Weights  #
# ---------------------------------------------------------------------- #

  # ---------------
  # Preparation
  # ---------------

  # Specify number of populations/datasets
  n.studies <- 4
  
  # Define three-dimensional arrays to store updated PMPs and Gorica Weights
  probs.arr <- array(0,dim=c(n.studies,n.hypos+1,n.iter.1))
  weights.arr <- array(0,dim=c(n.studies,n.hypos+1,n.iter.1))
    
  # Define matrix to store final PMPs of upcoming updating process
  probs.final <- matrix(NA, nrow=n.iter.1, ncol=n.hypos+1)
  if (n.hypos==1){colnames(probs.final) <- c("P(H1)", "P(Hu)")} else if (n.hypos==2){
                  colnames(probs.final) <- c("P(H1)", "P(H2)", "P(Hu)")} else if (n.hypos==3){
                  colnames(probs.final) <- c("P(H1)", "P(H2)", "P(H3)", "P(Hu)")} else if (n.hypos==4){
                  colnames(probs.final) <- c("P(H1)", "P(H2)", "P(H3)", "P(H4)", "P(Hu)")}
  
  # Define matrix to store final Gorica Weights of upcoming updating process
  weights.final <- matrix(NA, nrow=n.iter.1, ncol=n.hypos+1)
  if (n.hypos==1){colnames(weights.final) <- c("Weight(H1)", "Weight(Hu)")} else if (n.hypos==2){
                  colnames(weights.final) <- c("Weight(H1)", "Weight(H2)", "Weight(Hu)")} else if (n.hypos==3){
                  colnames(weights.final) <- c("Weight(H1)", "Weight(H2)", "Weight(H3)", "Weight(Hu)")} else if (n.hypos==4){
                  colnames(weights.final) <- c("Weight(H1)", "Weight(H2)", "Weight(H3)", "Weight(H4)", "Weight(Hu)")}
  
  # --------------
  # Application
  # --------------
  for (i in 1:n.iter.1){
    
    # Bayesian Analysis using Bayes Factors
    
      # Fill all matrix elements with NAs
      probs.arr[,,i] <- matrix(NA, nrow=n.studies, ncol=n.hypos+1)
      if (n.hypos==1){colnames(probs.arr) <- c("P(H1)", "P(Hu)")} else if (n.hypos==2){
                      colnames(probs.arr) <- c("P(H1)", "P(H2)", "P(Hu)")} else if (n.hypos==3){
                      colnames(probs.arr) <- c("P(H1)", "P(H2)", "P(H3)", "P(Hu)")} else if (n.hypos==4){
                      colnames(probs.arr) <- c("P(H1)", "P(H2)", "P(H3)", "P(H4)", "P(Hu)")}
      
      # Fill first matrix rows with 1/3
      probs.arr[1,,i] <- 1/3
      
      # Combine probabilities (this array stores the complete updating process)
      probs.arr[1,,i] <- (probs.arr[1,,i] * bfmu.1[i,]) / (sum(probs.arr[1,,i] * bfmu.1[i,]))
      probs.arr[2,,i] <- (probs.arr[1,,i] * bfmu.2[i,]) / (sum(probs.arr[1,,i] * bfmu.2[i,]))
      probs.arr[3,,i] <- (probs.arr[2,,i] * bfmu.3[i,]) / (sum(probs.arr[2,,i] * bfmu.3[i,]))
      probs.arr[4,,i] <- (probs.arr[3,,i] * bfmu.4[i,]) / (sum(probs.arr[3,,i] * bfmu.4[i,]))
      
      # Extract solely final probabilities
      probs.final[i,] <- c(probs.arr[4,,i])
      
    
    # Frequentist Analysis using Gorica Weights
      
      # Fill all matrix elements with NAs
      weights.arr[,,i] <- matrix(NA, nrow=n.studies, ncol=n.hypos+1)
      if (n.hypos==1){colnames(weights.arr) <- c("Weight(H1)", "Weight(Hu)")} else if (n.hypos==2){
                      colnames(weights.arr) <- c("Weight(H1)", "Weight(H2)", "Weight(Hu)")} else if (n.hypos==3){
                      colnames(weights.arr) <- c("Weight(H1)", "Weight(H2)", "Weight(H3)", "Weight(Hu)")} else if (n.hypos==4){
                      colnames(weights.arr) <- c("Weight(H1)", "Weight(H2)", "Weight(H3)", "Weight(H4)", "Weight(Hu)")}
      
      # Combine weights (this array stores the complete updating process)
      weights.arr[1,,i] <- weights.1[i,]
      weights.arr[2,,i] <- (weights.arr[1,,i] * weights.2[i,]) / (sum(weights.arr[1,,i] * weights.2[i,]))
      weights.arr[3,,i] <- (weights.arr[2,,i] * weights.3[i,]) / (sum(weights.arr[2,,i] * weights.3[i,]))
      weights.arr[4,,i] <- (weights.arr[3,,i] * weights.4[i,]) / (sum(weights.arr[3,,i] * weights.4[i,]))
      
      # Extract solely final probabilities
      weights.final[i,] <- weights.arr[4,,i]  
   
  }
  
  
  # -------------
  # Evaluation
  # -------------
  
  # Inspect final PMPs and Weights of all iterations
  probs.final
  weights.final
  
  # Plot final PMPs and Weights
  par(mfrow=c(1,ncol(probs.final)))
  for (i in 1:ncol(probs.final)){
    hist(probs.final[,i], main=colnames(probs.final)[i], xlab="", breaks=n.iter.1)
  }
  
  par(mfrow=c(1,(n.hypos+1)))
  for (i in 1:(n.hypos+1)){
    hist(weights.final[,i], main=colnames(weights.final)[i], xlab="", breaks=n.iter.1)
  }
  
  
  
  # Compute True Hypothesis Rate (given that H1 is true!)
  
    # Based on PMPs
    true.h.rate.pmp <- rep(NA, n.iter.1)
    for (i in 1:n.iter.1){
      if (n.hypos==1){
        if (probs.arr[4,1,i] > probs.arr[4,2,i])
            true.h.rate.pmp[i] <- 1 else {
            true.h.rate.pmp[i] <- 0 }}
      if (n.hypos==2){
        if (probs.arr[4,1,i] > probs.arr[4,2,i] & probs.arr[4,3,i])
          true.h.rate.pmp[i] <- 1 else {
          true.h.rate.pmp[i] <- 0 }}
      if (n.hypos==3){
        if (probs.arr[4,1,i] > probs.arr[4,2,i] & probs.arr[4,3,i] & probs.arr[4,4,i])
          true.h.rate.pmp[i] <- 1 else {
          true.h.rate.pmp[i] <- 0 }}
      if (n.hypos==4){
        if (probs.arr[4,1,i] > probs.arr[4,2,i] & probs.arr[4,3,i] & probs.arr[4,4,i] & probs.arr[4,5,i])
          true.h.rate.pmp[i] <- 1 else {
            true.h.rate.pmp[i] <- 0 }}
    }
    
    # Based on Gorica Weights
    true.h.rate.weights <- rep(NA, n.iter.1)
    for (i in 1:n.iter.1){
        if (n.hypos==2){
          if (weights.arr[4,1,i] > weights.arr[4,2,i] & weights.arr[4,3,i] )
              true.h.rate.weights[i] <- 1 else {
              true.h.rate.weights[i] <- 0 }}
      if (n.hypos==3){
        if (weights.arr[4,1,i] > weights.arr[4,2,i] & weights.arr[4,3,i] & weights.arr[4,4,i])
            true.h.rate.weights[i] <- 1 else {
            true.h.rate.weights[i] <- 0 }}
      if (n.hypos==4){
        if (weights.arr[4,1,i] > weights.arr[4,2,i] & weights.arr[4,3,i] & weights.arr[4,4,i] & weights.arr[4,5,i])
            true.h.rate.weights[i] <- 1 else {
            true.h.rate.weights[i] <- 0 }}
    }
    
    
  
  # Compute Ratios 
    
    # Based on PMPs
    pmp.ratios <- matrix(NA, nrow=n.iter.1, ncol=((n.hypos+1)*((n.hypos+1)-1))/2)
    if (n.hypos==2){
      colnames(pmp.ratios) <- c("P(H1)/P(H2)", "P(H1)/P(Hu)", "P(H2)/P(Hu)")
      for (i in 1:n.iter.1){
        pmp.ratios[i,1] <- probs.arr[4,1,i]/probs.arr[4,2,i]
        pmp.ratios[i,2] <- probs.arr[4,1,i]/probs.arr[4,3,i]
        pmp.ratios[i,3] <- probs.arr[4,2,i]/probs.arr[4,3,i]
      }}
    if (n.hypos==3){
      colnames(pmp.ratios) <- c("P(H1)/P(H2)", "P(H1)/P(H3)", "P(H1)/P(Hu)", "P(H2)/P(H3)", "P(H2)/P(Hu)", "P(H3)/P(Hu)")
      for (i in 1:n.iter.1){
        pmp.ratios[i,1] <- probs.arr[4,1,i]/probs.arr[4,2,i]
        pmp.ratios[i,2] <- probs.arr[4,1,i]/probs.arr[4,3,i]
        pmp.ratios[i,3] <- probs.arr[4,1,i]/probs.arr[4,4,i]
        pmp.ratios[i,4] <- probs.arr[4,2,i]/probs.arr[4,3,i]
        pmp.ratios[i,5] <- probs.arr[4,2,i]/probs.arr[4,4,i]
        pmp.ratios[i,6] <- probs.arr[4,3,i]/probs.arr[4,4,i]
      }}
    if (n.hypos==4){
      colnames(pmp.ratios) <- c("P(H1)/P(H2)", "P(H1)/P(H3)", "P(H1)/P(H4)", "P(H1)/P(Hu)", "P(H2)/P(H3)", "P(H2)/P(H4)", "P(H2)/P(Hu)", "P(H3)/P(H4)", "P(H3)/P(Hu)", "P(H4)/P(Hu)" )
      for (i in 1:n.iter.1){
        pmp.ratios[i,1] <- probs.arr[4,1,i]/probs.arr[4,2,i]
        pmp.ratios[i,2] <- probs.arr[4,1,i]/probs.arr[4,3,i]
        pmp.ratios[i,3] <- probs.arr[4,1,i]/probs.arr[4,4,i]
        pmp.ratios[i,4] <- probs.arr[4,1,i]/probs.arr[4,5,i]
        pmp.ratios[i,5] <- probs.arr[4,2,i]/probs.arr[4,3,i]
        pmp.ratios[i,6] <- probs.arr[4,2,i]/probs.arr[4,4,i]
        pmp.ratios[i,7] <- probs.arr[4,2,i]/probs.arr[4,5,i]
        pmp.ratios[i,8] <- probs.arr[4,3,i]/probs.arr[4,4,i]
        pmp.ratios[i,9] <- probs.arr[4,3,i]/probs.arr[4,5,i]
        pmp.ratios[i,10] <- probs.arr[4,4,i]/probs.arr[4,5,i]
        
      }}
  
    # Based on Gorica Weights
      
      weights.ratios <- matrix(NA, nrow=n.iter.1, ncol=((n.hypos+1)*((n.hypos+1)-1))/2)
      if (n.hypos==2){
        colnames(weights.ratios) <- c("Weight(H1)/Weight(H2)", "Weight(H1)/Weight(Hu)", "Weight(H2)/Weight(Hu)")
        for (i in 1:n.iter.1){
          weights.ratios[i,1] <- weights.arr[4,1,i]/weights.arr[4,2,i]
          weights.ratios[i,2] <- weights.arr[4,1,i]/weights.arr[4,3,i]
          weights.ratios[i,3] <- weights.arr[4,2,i]/weights.arr[4,3,i]
        }}
      if (n.hypos==3){
        colnames(weights.ratios) <- c("Weight(H1)/Weight(H2)", "Weight(H1)/Weight(H3)", "Weight(H1)/Weight(Hu)", "Weight(H2)/Weight(H3)", "Weight(H2)/Weight(Hu)", "Weight(H3)/Weight(Hu)")
        for (i in 1:n.iter.1){
          weights.ratios[i,1] <- weights.arr[4,1,i]/weights.arr[4,2,i]
          weights.ratios[i,2] <- weights.arr[4,1,i]/weights.arr[4,3,i]
          weights.ratios[i,3] <- weights.arr[4,1,i]/weights.arr[4,4,i]
          weights.ratios[i,4] <- weights.arr[4,2,i]/weights.arr[4,3,i]
          weights.ratios[i,5] <- weights.arr[4,2,i]/weights.arr[4,4,i]
          weights.ratios[i,6] <- weights.arr[4,3,i]/weights.arr[4,4,i]
        }}
      if (n.hypos==4){
        colnames(weights.ratios) <- c("Weight(H1)/Weight(H2)", "Weight(H1)/Weight(H3)", "Weight(H1)/Weight(H4)", "Weight(H1)/Weight(Hu)", "Weight(H2)/Weight(H3)", "Weight(H2)/Weight(H4)", "Weight(H2)/Weight(Hu)", "Weight(H3)/Weight(H4)", "Weight(H3)/Weight(Hu)", "Weight(H4)/Weight(Hu)" )
        for (i in 1:n.iter.1){
          weights.ratios[i,1] <- weights.arr[4,1,i]/weights.arr[4,2,i]
          weights.ratios[i,2] <- weights.arr[4,1,i]/weights.arr[4,3,i]
          weights.ratios[i,3] <- weights.arr[4,1,i]/weights.arr[4,4,i]
          weights.ratios[i,4] <- weights.arr[4,1,i]/weights.arr[4,5,i]
          weights.ratios[i,5] <- weights.arr[4,2,i]/weights.arr[4,3,i]
          weights.ratios[i,6] <- weights.arr[4,2,i]/weights.arr[4,4,i]
          weights.ratios[i,7] <- weights.arr[4,2,i]/weights.arr[4,5,i]
          weights.ratios[i,8] <- weights.arr[4,3,i]/weights.arr[4,4,i]
          weights.ratios[i,9] <- weights.arr[4,3,i]/weights.arr[4,5,i]
          weights.ratios[i,10] <- weights.arr[4,4,i]/weights.arr[4,5,i]
         
        }}
  
  
  par(mfrow=c(1,ncol(pmp.ratios)))
  for (i in 1:ncol(pmp.ratios)){
    hist(pmp.ratios[,i], main=colnames(pmp.ratios)[i], xlab="", breaks=n.iter.1)
  }
  
  par(mfrow=c(1,ncol(weights.ratios)))
  for (i in 1:ncol(weights.ratios)){
    hist(weights.ratios[,i], main=colnames(weights.ratios)[i], xlab="", breaks=n.iter.1)
  }
  
   
  # How to inspect results?
  
    # Histogram 1: Distribution of PMPs 
    # Histogram 2: Distribution of Gorica Weights
    # Histogram 3: Distribution of PMP Ratios
    # Histogram 4: Distribution of Gorica Weight Ratios
    mean(true.h.rate.pmp) # THR based on PMPs
    mean(true.h.rate.weights) # THR based on Gorica 
    





















