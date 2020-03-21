#' Calculate stepwise PIC
#'
#' @param x       A vector of response.
#' @param py      A matrix containing possible predictors of x.
#' @param alpha   The significance level used to judge whether the sample estimate in Equation \deqn{\hat{PIC} = sqrt(1-exp(-2\hat{PI})} is significant or not. A default alpha value is 0.1.
#' @param method  A logic value. Use Fortran (T) or R (F). Default T.
#'
#' @return A list of 2 elements: the column numbers of the meaningful predictors (cpy), and partial informational correlation (cpyPIC).
#' @export
#'
#' @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
#' @examples
#' data(data1) # AR9 model   x(i)=0.3*x(i-1)-0.6*x(i-4)-0.5*x(i-9)+eps
#' x<-data1[,1]   # response
#' py<-data1[,-1]  # possible predictors
#' stepwise.PIC(x,py)
#' 
#' data(data2) # AR4 model:  x(i)=0.6*x(i-1)-0.4*x(i-4)+eps
#' x<-data2[,1]   # response
#' py<-data2[,-1]  # possible predictors
#' stepwise.PIC(x,py)
#' 
#' data(data3) # AR1 model  x(i)=0.9*x(i-1)+0.866*eps
#' x<-data3[,1]   # response
#' py<-data3[,-1]  # possible predictors
#' stepwise.PIC(x,py)

stepwise.PIC <- function (x, py, alpha = 0.1) 
{
  x = as.matrix(x)
  n = nrow(x)
  npy = ncol(py)
  cpy = cpyPIC = NULL
  icpy = 0
  z = NULL
  isig = T
  icoloutz = 1:npy
  #cat("calc.PIC-----------","\n")
  while (isig) {
    npicmax = npy - icpy
    pictemp = rep(0, npicmax)
    y = py[, icoloutz]
    pictemp = pic.calc(x,as.matrix(y),z)
    ctmp = order(-pictemp)[1]
    cpytmp = icoloutz[ctmp]
    picmaxtmp = pictemp[ctmp]
    if (!is.null(z)) {
      df = n - ncol(z)
    } else {
      df = n
    }
    
    t <-  qt(1-alpha, df=df)
    picthres <- sqrt(t^2/(t^2+df)) 
    #cat("picthres",picthres,"\n")
    
    if (picmaxtmp > picthres) {
      cpy = c(cpy, cpytmp)
      cpyPIC = c(cpyPIC, picmaxtmp)
      #cat(cpyPIC,"\n")
      z = cbind(z, py[, cpytmp])
      icoloutz = icoloutz[-ctmp]
      icpy = icpy + 1
      if ((npy - icpy) == 0) isig = F
    } else {
      isig = F
    }
  }
  #cat("calc.PW------------","\n")
  if (!is.null(z)) {
    out = pw.calc(x, py, cpy, cpyPIC)
    outwt = out$pw
    lstwt = abs(lsfit(z, x)$coef)
    return(list(cpy = cpy, cpyPIC = cpyPIC, wt = outwt, 
                lstwet = lstwt))
  } else {
    message("None of the provided predictors is related to the response variable")
  }
}

#' Calculate the ratio of conditional error standard deviations
#'
#' @param x     A vector of response.
#' @param zin   A matrix containing the meaningful predictors selected from a large set of possible predictors (z).
#' @param zout  A matrix containing the remaining possible predictors after taking out the meaningful predictors (zin).
#'
#' @return The STD ratio.
#' @export
#'
#' @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
calc.scaleSTDratio <- function (x, zin, zout) 
{
  if(!missing(zout)){
    xhat = knnregl1cv(x, zout)
    stdratxzout = sqrt(var(x - xhat)/var(x))
    zinhat = knnregl1cv(zin, zout)
    stdratzinzout = sqrt(var(zin - zinhat)/var(zin))
    return(0.5 * (stdratxzout + stdratzinzout))
  } else {
    return(1) 
  }

}

#' Leave one out cross validation.
#'
#' @param x A vector of response.
#' @param z A matrix of predictors.
#' @param k The number of nearest neighbours used. The default is 0, indicating Lall and Sharma default is used.
#' @param pw A vector of partial weights of the same length of z.
#'
#' @return A vector of L1CV estimates of the response.
#' @export
#'
#' @references Lall, U., Sharma, A., 1996. A Nearest Neighbor Bootstrap For Resampling Hydrologic Time Series. Water Resources Research, 32(3): 679-693.
#' @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
knnregl1cv <- function (x, z, k = 0, pw) 
{
  x = as.matrix(x)
  n = nrow(x)
  if (k == 0) 
    k = floor(0.5 + 3 * (sqrt(n)))
  z = as.matrix(z)
  nz = ncol(z)
  sd = sqrt(diag(var(z)))
  if (missing(pw)) {
    pw = rep(1, nz)
  }
  for (j in 1:nz) z[, j] = z[, j]/(sd[j]/pw[j])
  d = as.matrix(dist(z))
  ord1 = apply(d, 2, order)
  ord = ord1[2:(k + 1), ]
  kern = 1/(1:k)/sum(1/(1:k))
  xhat = rep(0, n)
  for (j in 1:k) xhat = xhat + x[ord[j, ]] * kern[j]
  return(xhat)
}

#-------------------------------------------------------------------------------
kernel.est.uvn <- function(Z) {
  
  N <- length(Z)
  d <- 1
  # compute sigma & constant
  sigma <- 1.5*bw.nrd0(Z)
  constant <- sqrt(2*pi) * sigma * N
  
  # Commence main loop
  dens <- vector()
  for(h in 1:N) {
    dis.Z <- (Z - Z[h])^2
    exp.dis <- exp(-dis.Z / (2*sigma^2))
    dens[h] <- sum(exp.dis) / constant
  }
  return(dens)
}
#-------------------------------------------------------------------------------
kernel.est.mvn <- function(Z) {
  
  # Compute covariance and determinant of cov
  N <- nrow(Z)
  d <- ncol(Z)
  Cov <- cov(Z)
  det.Cov <- det(Cov)

  # Compute sigma & constant
  sigma <- 1.5 * (4/(d + 2))^(1/(d + 4)) * N^(-1/(d + 4))
  constant <- (sqrt(2*pi)*sigma)^d * sqrt(det.Cov) * N
  
  # Commence main loop
  dens <- vector()
  for(h in 1:N) {
    dist.val <- mahalanobis(Z, center = Z[h,], cov = Cov)
    exp.dis <- exp(-dist.val / (2*sigma^2))
    dens[h] <- sum(exp.dis) / constant
  }
  
  return(dens)
}
#-------------------------------------------------------------------------------
pmi.calc <- function(Y, X) {
  # Use Gaussian kernel to compute marginal densities F(Y) & F(X) and
  # joint density F(X,Y)
  N <- length(Y)
  pdf.X <- kernel.est.uvn(X)
  
  pdf.Y <- kernel.est.uvn(Y)
  pdf.XY <- kernel.est.mvn(cbind(Y, X))
  
  calc <- log(pdf.XY / (pdf.Y * pdf.X))
  return(sum(calc) / N)
  
}
#-------------------------------------------------------------------------------
pic.calc <- function(Y, X, Z) {
  
  if(is.null(Z)){
    x.in <- X
    y.in <- Y
  } else {
    x.in <- apply(X, 2, function(x) knnregl1cv(x, Z)-x)
    y.in <- knnregl1cv(Y, Z)-Y
  }
  
  pmi <- apply(x.in, 2, function(z) pmi.calc(y.in,z))
  pmi[pmi<0] <- 0
  pic <- sqrt(1-exp(-2*pmi))
  
  return(as.numeric(pic))
}
#-------------------------------------------------------------------------------
pw.calc <- function(Y, X, cpy, cpyPIC){
  
  wt <- NA
  Z = as.matrix(X[,cpy])
  if(ncol(Z)==1) {
    wt <- calc.scaleSTDratio(Y, Z)*cpyPIC
  } else {
    for(i in 1:length(cpy)) wt[i] <- calc.scaleSTDratio(Y, Z[,i], Z[,-i])*cpyPIC[i]
  }
  return(list(pw=wt))
}


