##=======================================================
## Confidence Intervals for
##     (a) Ellipse: P-dimensional normal distribution
##     (b) Difference between two means
##-------------------------------------------------------
## Simultaneous Confidence intervals
##    (1) Single mean
##    (2) Bonferroni
##    (3) Large sample
##    (4) Large sample ellipses
##-------------------------------------------------------
##=======================================================

common.ftable <- function(n, p, alpha=0.05)
  {
    warning("alpha=",alpha)
    (((n - 1) * p)/(n - p))*qf(1 - alpha, p, n - p)
  }

common.chisq <- function(df, alpha=0.05)
  {
    warning("alpha=",alpha)
    qchisq(1-alpha, df)
  }

common.t <- function(n, p, alpha=0.05)
  {
    warning("alpha=",alpha)
    qt(1-(alpha)/(2*p),n-1)
  }

common.pooled_varcov <- function(data1, data2)
  {
    sp1 <- cov(data1)
    n1 <- nrow(data1)
    
    sp2 <- cov(data2)
    n2 <- nrow(data2)

    sp <- ((n1 - 1)*sp1 + (n2 - 1)*sp2)/(n1 + n2 - 2)
  }

common.Fhotelling <- function(q, n1, n2, alpha=0.05)
  {
    warning("alpha=", alpha)
    qf(1-(alpha),q, (n1 + n2 - q - 1))
  }

common.hotelling.means <- function(means1, gmeans, varcov,
                                   n, alpha=0.05, rnd=3)
  {
    
    ## Mahalanobi distance
    D2 <- t(means1 - gmeans)%*%solve(varcov)%*%(means1 - gmeans)
    ## T-squared
    T2 <- n*D2
    ## Number of parameters
    p <- ncol(varcov)

    ## The calculated F-statistic
    fcalc <- common.ftable(n, p, alpha)

    fcalc <- round(fcalc, digits=3)
    T2 <- round(T2, digits=3)
    message("\tHotelling T^2: ", T2,
            " compare with F(p,n-p,alpha=",alpha,
            ") = ", fcalc)


  }

common.hotelling.data <- function(data1, data2, alpha=0.05, rnd=3)
  {
    sp <- common.pooled_varcov(data1, data2)
    
    n1 <- nrow(data1)
    n2 <- nrow(data2)
    
    means1 <- apply(data1, 2, mean)
    means2 <- apply(data2, 2, mean)

    ## Mahalanobi distance
    D <- t(means1 - means2)%*%solve(sp)%*%(means1 - means2)
    ## T-squared
    T <- (n1*n2)*D/(n1 + n2)
    ## Number of parameters
    q <- ncol(sp)

    ## The calculated F-statistic
    fcalc <- (n1 + n2 - q - 1)*T/((n1 + n2 - 2)*q)
    fcalc <- round(fcalc, digits=3)
    
    fval <- common.Fhotelling(q, n1, n2, alpha)
    fval <- round(fval, digits=3)
    
    message("\tHotelling T^2: ", fcalc,
            " compare with F(q,n1+n2-q-1) = ",
            fval)
  }

common.paired <- function(data1, data2, hypomeans,
                          alpha=0.05, rnd=3)
  {
    ## Paired comparison
    data.diff <- data1 - data2
    data.mean <- apply(data.diff, 2, mean)

    data.cov <- cov(data.diff)

    # Conduct hypothesis with U_0 = [X, Y, Z,...]
    hmean <- matrix(hypomeans,ncol=1)

    # Number of observations
    n <- nrow(data.diff)
    p <- ncol(data.diff)
    # The T-Squared test statistic
    data.tsq <- n*t(data.mean - hmean)%*%solve(data.cov)%*%(data.mean - hmean)
    data.tsq <- round(data.tsq, digits=rnd)
    
    #Compare with the following F-statistic
    data.fcalc <- common.ftable(n, p, alpha)
    data.fcalc <- round(data.fcalc, digits=rnd)

    message("\tDependent Comparison of Means\n",
            "Compare: T^2=",data.tsq," with ",
            "F-table=",data.fcalc)
  }

common.dependent <- function(...)
  {
    common.paired(...)
  }
common.repeated <- function(...)
  {
    common.paired(...)
  }

ci.ellipse.pdim <- function(means, hmeans, varcov,
                            n, alpha=0.05, rnd=3)
  {
    
    mean.diff <- means - hmeans
    mean.diff <- matrix(mean.diff, ncol=1)

    p <- nrow(varcov)

    ftab <- common.ftable(n, p, alpha)
    ftab <- round(ftab, digits=rnd)

    cie <- n*t(mean.diff)%*%solve(varcov)%*%(mean.diff)
    cie <- round(cie, digits=rnd)

    means <- round(means, digits=3)
    hmeans <- round(hmeans, digits=3)
    
    message("\n\t100(1-",alpha,")% Confidence ellipse of",
            " a P-Dimensional normal distribution.\n",
            " for _X_ = (", toString(means),")",
            " == MU = (", toString(hmeans),")",
            "\n : Ellipse is { ", cie,
            " <= [(((n - 1)*p)/(n - p))F(df1=",p,"df2=",n - p,
            ",alpha=",alpha,") = ",
            ftab,"] }\n\n")
  }

ci.ellipse.twomeans <- function(means1, means2,
                                varcov, n, alpha=0.05, rnd=3)
  {
    p <- nrow(varcov)
    
    fval <- common.ftable(n, p, alpha)
    fval <- round(fval, digits=rnd)
    
    for(i in 1:(p-1)){
      
      xbar1 <- means1[i]
      mu1 <- means2[i]
      
      for(j in (i+1):p){

        xbar2 <- means1[j]
        mu2 <- means2[j]

        dbar <- c(xbar1, xbar2) - c(mu1, mu2)
        dbar <- matrix(dbar, ncol=1)
        
        sii <- varcov[i,i]
        sij <- varcov[i,j]
        sjj <- varcov[j,j]
        smat <- matrix(c(sii, sij, sij, sjj),ncol=2)
        
        cie <- n*t(dbar)%*%solve(smat)%*%dbar
        cie <- round(cie, digits=rnd)

        mu1 <- round(mu1, digits=rnd)
        mu2 <- round(mu2, digits=rnd)
        
        message("\t\t100(1-",alpha,")% Confidence ellipse",
                " of two means.\n",
                "(U_", i,",U_", j,
                ") = (", toString(c(mu1, mu2)),")",
                " : ", cie,
                " <= [(((n - 1)*p)/(n - p))F(df1=",p,",df2=",n - p,
                ",alpha=",alpha,") = ",
                fval,"]")        
      }
    }
  }



ci.ellipse.largesample <- function(means1, means2,
                                   varcov, n, alpha=0.05, rnd=3)
  {
    p <- nrow(varcov)

    Xsq <- common.chisq(p, alpha)
    Xsq <- round(Xsq, digits=rnd)
    
    for(i in 1:(p-1)){
      xb1 <- means1[i]
      mu1 <- means2[i]
      
      for(j in (i+1):p){
        xb2 <- means1[j]
        mu2 <- means2[j]

        diffb <- c(xb1, xb2) - c(mu1, mu2)
        diffb <- matrix(diffb, ncol=1)

        sii <- varcov[i,i]
        sij <- varcov[i,j]
        sjj <- varcov[j,j]
        smat <- matrix(c(sii, sij, sij, sjj),ncol=2)
        
        cie <- n*t(diffb)%*%solve(smat)%*%diffb
        cie <- round(cie, digits=rnd)

        mu1 <- round(mu1, digits=rnd)
        mu2 <- round(mu2, digits=rnd)
        
        message("\t100(1-",alpha,")% Large-sample Confidence",
                " ellipse.\n",
                "(U_", i,",U_", j,
                ") = (", toString(c(mu1,mu2)),")",
                " : ", cie,
                " <= [Chisq(df=",p, ",alpha=",alpha,") = ",
                Xsq,"]")
      }
    }
  }


##-------------------------------------------------------
## varcov: the variance-covariance matrix
## n:      number of observations.
##-------------------------------------------------------
ci.simult.singlemu <- function(means, varcov, n, alpha=0.05, rnd=3)
  {
    
    p <- nrow(varcov)
    fval <- common.ftable(n, p, alpha)
    fval <- sqrt(fval)

    for(i in 1:p){
      sii <- varcov[i,i]
      sii <- sqrt(sii/n)

      xbi <- means[i]
      
      interval <- fval*sii
      interval <- xbi + c((-1)*interval,interval)
      interval <- round(interval, digits=rnd)

      message("\tT^2 Simultaneous CI - Mu = ",
              "U_",i," : (",toString(interval),")")

    }
  }

##-------------------[CORRECT]---------------------------
## means: the means of the p-factors
## varcov: the variance-covariance matrix
## n:      number of observations.
##-------------------------------------------------------
ci.simult.bonferroni <- function(means, varcov,
                                 n, alpha=0.05, rnd=3)
  {
    p <- nrow(varcov)
    tval <- common.t(n, p, alpha)
    for(i in 1:p){
      xbar <- means[i]
      
      sii <- (varcov[i,i])/n
      sii <- sqrt(sii)

      interval <- tval*sii
      
      interval <- c((-1)*interval, interval)
      interval <- xbar + interval
      interval <- round(interval, digits=rnd)
      
      message("\tBonferroni CI - Mu=U_",i,
              " : (",toString(interval),")")
    }
  }

##-------------------------------------------------------
## means1: the means of the p-factors
## means2: the given means
## varcov: the variance-covariance matrix
## n:      number of observations.
##-------------------------------------------------------
ci.simult.largesample <- function(means, varcov,
                                  n, alpha=0.05, rnd=3)
  {
    p <- nrow(varcov)
    Xsq <- common.chisq(p, alpha)
    Xsq <- sqrt(Xsq)
    
    for(i in 1:p){
      xbar <- means[i]
      
      sii <- (varcov[i,i])/n
      sii <- sqrt(sii)

      interval <- Xsq*sii
      
      interval <- c((-1)*interval, interval)
      interval <- xbar + interval
      interval <- round(interval, digits=rnd)

      message("\tSimultaneous Large-sample C.I - Mu=U_",i,
              " : (",toString(interval),")")
    }
  }

ci.simult.diff <- function(means, varcov, n, alpha=0.05, rnd=3)
  {
    p <- nrow(varcov)

    Fval <- common.ftable(n, p, alpha)
    Fval <- sqrt(Fval)
    
    for(i in 1:(p-1)){
      
      xbar_i<- means[i]
      
      for(j in (i+1):p){
          
        xbar_j<- means[j]
        diff_bar <- xbar_i - xbar_j

        # Variance-Covariance components
        sii <- varcov[i,i]
        sij <- varcov[i,j]
        sjj <- varcov[j,j]
        ss <- (sii + sjj - (2 * sij))/n
        ss <- sqrt(ss)

        interval <- c((-1)*Fval*ss, Fval*ss)
        interval <- diff_bar + interval
        # Round off to specified significant digits
        interval <- round(interval, digits=rnd)
        
        message("\tC.I Difference of 2 means",
                " (U_", i," - U_", j,
                ") = (",toString(interval),")")
        
      }
    }
  }
