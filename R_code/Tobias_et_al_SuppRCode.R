# 0. Front Matter ####
# Purpose: 
#     This is the code for the analysis in 
#     Tobias, V.D., E. Chen, J. Hobbs, M. Eakin, and S. Detwiler. in press. 
#          Biological Conservation. Preprint DOI: 
#          https://doi.org/10.1101/2023.04.17.537200 
# Author: 
#     Vanessa Tobias <vanessa_tobias@fws.gov>
# USFWS Disclaimer:
#     Code written by staff at the United States Fish and Wildlife Service (FWS) 
#     is provided on an “as is” basis and the user assumes responsibility for 
#     its use. FWS has relinquished control of the information and no longer 
#     has responsibility to protect the integrity, confidentiality, or 
#     availability of the information. Any reference to specific commercial 
#     products, processes, or services by service mark, trademark, manufacturer, 
#     or otherwise, does not constitute or imply their endorsement, 
#     recommendation or favoring by FWS. The FWS seal and logo shall not be 
#     used in any manner to imply endorsement of any commercial product or 
#     activity by FWS or the United States Government.
# License: 
#     This project is licensed under the terms of the GNU General public
#     license, Version 3.0.
# Format:
#     When this code is opened with RStudio, it will have a table of contents.
#     The user should set the root folder for figure output in section 1.



# 1. Setup ####
## 1.1 Load Packages ####
library(popbio)
library(lmtest)
library(plotrix)
library(dplyr)
library(meta)
library(MASS)
library(MuMIn)
library(trend)

## 1.2 Set folders ####
root_figs <- "./Figures/Revised_20231019/"

## 1.3 Load Data ####
dat <- read.csv("./DataProcessed/LFS_Indices_oldnames.csv",
                header = TRUE,
                stringsAsFactors = FALSE)

surveyList <- c("mwt.0", "mwt.1", "mwt.2", # Bay Study MWT
                "ot.0", "ot.1", "ot.2",    # Bay Study OT
                "Total",                   # FMWT
                "index")                   # 20mm

## 1.4 Define Functions ####
negloglike.ml = function(theta, yt, tt)  
{
  # Code originally from Humbert et al. 2009
  # A better way to estimate population trends
  # Oikos 118: 1940-1946, 2009
  # doi: 10.1111/j.1600-0706.2009.17839.x
  
  #  ML objective function "negloglike.ml" is negative of log-likelihood;
  #  the Nelder-Mead optimization routine in R, "optim", is a minimization
  #  routine.  The ML objective function uses equations 24-26 from Dennis et
  #  al. (2006).  The three function arguments are:  theta, vector of
  #  parameters (transformed to the real line), yt, vector of time series
  #  observations, and tt, vector of observation times.
  muu = theta[1];
  sigmasq = exp(theta[2]);      #  Constrains ssq > 0. 
  tausq = exp(theta[3]);        #  Constrains tsq > 0.
  xzero = theta[4];
  q = length(yt)-1;
  qp1 = q+1;
  yt = matrix(yt, nrow = qp1, ncol = 1);
  vx = matrix(0, qp1, qp1);
  for (ti in 1:q)
  {
    vx[(ti+1):qp1,(ti+1):qp1] = matrix(1, 1, (qp1-ti)) * tt[ti+1];
  }
  Sigma.mat=sigmasq*vx;
  Itausq=matrix(rep(0,(qp1*qp1)),nrow=q+1,ncol=q+1);
  diag(Itausq)=rep(tausq,q+1);
  V=Sigma.mat+Itausq;
  mu=matrix((xzero+muu*tt),nrow=qp1,ncol=1);
  ofn=((qp1)/2)*log(2*pi)+(0.5*log(det(V)))+
    (0.5*(t(yt-mu)%*%ginv(V)%*%(yt-mu)));
  return(ofn);
}



negloglike.reml=function(theta,yt,tt)
{
  # Code originally from Humbert et al. 2009
  # A better way to estimate population trends
  # Oikos 118: 1940-1946, 2009
  # doi: 10.1111/j.1600-0706.2009.17839.x
  
  #  REML objective function "negloglike.reml" is negative of log-likelihood
  #  for second differences of the log-scale observations.  The REML objective
  #  function uses equations A18-A25 from Humbert et al. (2009).  The three
  #  function arguments are:  theta, vector of parameters (transformed to the
  #  real line), yt, vector of time series observations (log scale), and
  #  tt, vector of observation times.  Function performs the differencing.
  sigsq=exp(theta[1]);         #  Constrains ssq > 0.
  tausq=exp(theta[2]);         #  Constrains tsq > 0.
  q=length(yt)-1;
  qp1=q+1;
  vx=matrix(0,qp1,qp1);
  for (ti in 1:q)
  {
    vx[(ti+1):qp1,(ti+1):qp1]=matrix(1,1,(qp1-ti))*tt[ti+1];
  }
  Sigma.mat=sigsq*vx;
  Itausq=matrix(rep(0,(qp1*qp1)),nrow=q+1,ncol=q+1);
  diag(Itausq)=rep(tausq,q+1);
  V=Sigma.mat+Itausq;
  ss=tt[2:qp1]-tt[1:q];
  D1mat=cbind(-diag(1/ss),matrix(0,q,1))+cbind(matrix(0,q,1),diag(1/ss));
  D2mat=cbind(-diag(1,q-1),matrix(0,q-1,1))+
    cbind(matrix(0,q-1,1),diag(1,q-1));
  Phi.mat=D2mat%*%D1mat%*%V%*%t(D1mat)%*%t(D2mat);
  wt=(yt[2:qp1]-yt[1:q])/ss;
  ut=wt[2:q]-wt[1:q-1];
  ofn=(q/2)*log(2*pi)+(0.5*log(det(Phi.mat)))+
    (0.5*(ut%*%ginv(Phi.mat)%*%ut));
  return(ofn);
}

extCDFAC <- function (mu, 
                      sig2, 
                      Nc, 
                      Ne,
                      dd = "none", #density dependence: "none" or "ricker"
                      # "ricker" needs troubleshooting- shouldn't be able to go infinite
                      K = NA, #equilibrium population size
                      rho, #correlation coefficient
                      numReps, #number of sims to run or pops to simulate
                      tmax = 50) {
  # modified version of popbio::countCDFxt
  # that can account for autocorrelation in annual growth rates
  # using Morris & Doak box 4.6 for guidance
  
  sizeMat <- matrix(0, nrow = numReps, ncol = tmax)
  
  sig <- sqrt(sig2)
  beta1 <- sqrt(1-rho^2) #fixes total env var = sig2
  numExt <- 0 #start with no extinct populations
  for(i in 1:numReps){
    N <- Nc #start at current population size
    eold <- rnorm(n = 1) #random deviation drawn from standard normal distribution
    for(j in 1:tmax){
      # make a new autocorrelated deviation
      
      enew <- (rho * eold) + (sig * beta1 * rnorm(n = 1))
      # project one year ahead with Ricker model
      # maybe change to BH later?
      if(dd == "none"){
        N <- N * exp(mu + enew)
      } else if(dd == "ricker"){
        N <- N * exp(mu*(1-N/K) + enew)
      }
      
      eold <- enew #save this e for next iteration
      N <- if(N < Ne) 0 else N # if N is less than Ne, then crash
      sizeMat[i, j] <- N
    }
  }
  return(sizeMat)
}

rsq2 <- function(data, # matrix result of extCDFAC
                 n, # size of sample to take from data
                 B){ # number of bootstrap samples
  # Bootstrap confidence intervals for extCDFAC
  x <- matrix(NA, B, ncol(data))
  for(i in 1:B){
    # randomly sample row numbers
    s.vec <- sample(x = 1:nrow(data), n, replace = TRUE)
    # use vector of row numbers to select rows from the data
    # --> matrix (n X ncol)
    d <- data[s.vec,]
    # calculate the proportion of 0s in each column
    # --> vector of length ncol
    m.prop <- apply(d, 2, FUN = function(x) 1-(sum(x>0)/n))
    x[i,] <- m.prop
  }
  m.x <- apply(x, 2, mean)
  se.m <- apply(x, 2, sd)
  quant <- apply(x, 2, function(x) quantile(x,c(0.05, 0.5, 0.95)))
  ret <- data.frame(cbind(t(quant), m.x, se.m))
  names(ret) <- c("lower", "median", "upper", "mean", "SE")
  return(ret)
}



calcLam2 <- function(dataset, abCol, yearCol){
  # Function to calculate lambda with the diffusion approximation
  # skips over NA values and still calculate lambda correctly
  # -- original calcLam just didn't use pairs with NAs
  # -- TL;DR calcLam2 uses more of the existing data
  
  # Set up data
  dat <- dataset                             # rename input dataset
  
  year <- dat[, which(names(dat)==yearCol)]  # identify year data
  ab <- dat[, which(names(dat)==abCol)]      # identify abundance data
  
  dat <- data.frame("year" = year, "ab" = ab)                # combine year and abundance into a dataframe
  dat <- na.omit(dat)                        # remove records for missing abundance values
  
  dat$lam <- c(dat$ab[2:length(dat$ab)]/dat$ab[1:length(dat$ab)-1], NA)  # calculate aparent lambda from existing ab ests
  dat$loglam <- log(dat$lam)                     # calculate laglambda
  
  # Run regression analysis
  regdat <- data.frame(x = c(sqrt(dat$year[2:length(dat$year)]-dat$year[1:length(dat$year)-1]), NA))
  regdat$y <- dat$loglam/regdat$x
  regdat$y[which(is.infinite(regdat$y))] <- NA
  
  # reg <- glm(x ~ y + 0, data = na.omit(regdat[which(dat$year > 2000),]))
  reg <- glm(y ~ x + 0, data = na.omit(regdat))
  
  df <- summary(reg)$df[2]
  tcrit <- qt(0.05/2, df, lower.tail = FALSE)
  
  mu <- summary(reg)$coefficients[1]
  mu_se <- summary(reg)$coefficients[2]
  mu_lcl <- mu - tcrit*mu_se
  mu_ucl <- mu + tcrit*mu_se
  
  lambda <- exp(mu)
  lambda_lcl <- exp(mu_lcl)
  lambda_ucl <- exp(mu_ucl)
  
  chi025 <- qchisq(1-0.025, df=df)
  chi975 <- qchisq(1-0.975, df=df)  
  
  sigmasq <- mean(reg$residuals^2)
  sigmasq_lcl <- df*sigmasq/chi025
  sigmasq_ucl <- df*sigmasq/chi975
  
  out <- list(reg, mu, mu_se, mu_lcl, mu_ucl, sigmasq, sigmasq_lcl, sigmasq_ucl,
              lambda, lambda_lcl, lambda_ucl, dat[, c("year", "lam", "loglam", "ab")])
  names(out) <- c("reg", "mu", "mu_se", "mu_lcl", "mu_ucl", "sigmasq","sigmasq_lcl", "sigmasq_ucl",
                  "lambda", "lambda_lcl", "lambda_ucl", "datDF")
  return(out)
}


# 2. Individual Surveys ####
## 2.1 Set up data ####
# Make a table to hold results:
Tab2_EGSS_REML <- data.frame(Survey = surveyList,
                             mu.reml = NA,
                             mu_lo.reml = NA,
                             mu_hi.reml = NA)

# Make a table for summaries of mu & variability:
muDF_EGSS <- data.frame("survey" = surveyList,
                        "mu" = rep(NA,length(surveyList)),
                        "sd" = rep(NA,length(surveyList)),
                        "n"  = rep(NA,length(surveyList)))


# Start calculations #
# Much of the code in this section comes from Humbert et al. 2009,
#   with some modifications

for(i in surveyList){ 
  # Loop through all of the surveys
  # This is a very big chunk of code to loop through.
  
  # Prep the data for individual surveys
  my_data <- dat[c("Year", i)]
  # remove lines with NAs
  my_data <- my_data[complete.cases(my_data),]
  names(my_data) <- c("Time.t", "Observed.t")
  # zeros break the log part of the code. Replace them with 1's
  my_data$Observed.t[which(my_data$Observed.t == 0)] <- 1
  # pull out vectors for the years and observations
  Observed.t <- my_data$Observed.t
  Time.t <- my_data$Time.t

  # Set values 

  T.t = Time.t-Time.t[1];     #  Time starts at zero.
  Y.t =log(Observed.t);      #  Log-transform the observations.
  q = length(Y.t)-1;          #  Number of time series transitions, q.
  qp1 = q+1;                  #  q+1 gets used a lot, too.
  S.t = T.t[2:qp1]-T.t[1:q];  #  Time intervals.
  m = rep(1,qp1);             #  Will contain Kalman means for Kalman calculations.
  v = rep(1,qp1);             #  Will contain variances for Kalman calculations.

  # The EGOE estimates
  Ybar=mean(Y.t);
  Tbar=mean(T.t);
  mu.egoe=sum((T.t-Tbar)*(Y.t-Ybar))/sum((T.t-Tbar)*(T.t-Tbar));
  x0.egoe=Ybar-mu.egoe*Tbar;
  ssq.egoe=0;
  Yhat.egoe=x0.egoe+mu.egoe*T.t;
  tsq.egoe=sum((Y.t-Yhat.egoe)*(Y.t-Yhat.egoe))/(q-1);
  
  # The EGPN estimates
  Ttr=sqrt(S.t);
  Ytr=(Y.t[2:qp1]-Y.t[1:q])/Ttr;
  mu.egpn=sum(Ttr*Ytr)/sum(Ttr*Ttr);
  Ytrhat=mu.egpn*Ttr;
  ssq.egpn=sum((Ytr-Ytrhat)*(Ytr-Ytrhat))/(q-1);
  tsq.egpn=0;
  x0.egpn=Y.t[1];
  
  # Initial values for EGSS are averages of EGOE and EGPN values 
  mu0=(mu.egoe+mu.egpn)/2;    #  For ML only 
  ssq0=ssq.egpn/2;            #  For ML and REML
  tsq0=tsq.egoe/2;            #  For ML and REML
  x00=x0.egoe;                #  For ML only     
  
  
## 2.2 Calculate ML and REML Estimates ####
  
  # The ML estimates.
  EGSSml=optim(par=c(mu0,log(ssq0),log(tsq0),x00),
               negloglike.ml,NULL,method="Nelder-Mead",yt=Y.t,tt=T.t);
  params.ml=c(EGSSml$par[1],exp(EGSSml$par[2]),exp(EGSSml$par[3]),
              EGSSml$par[4]);
  lnlike.ml=-EGSSml$value[1];
  AIC.egss=-2*lnlike.ml+2*length(params.ml);
  
  mu.ml=params.ml[1];           # These are the ML estimates.
  ssq.ml=params.ml[2];          #          --
  tsq.ml=params.ml[3];          #          --
  x0.ml=params.ml[4];           #          --
  
  # The REML estimates.
  EGSSreml=optim(par=c(log(ssq0),log(tsq0)),
                 negloglike.reml,NULL,method="Nelder-Mead",yt=Y.t,tt=T.t);
  params.reml=c(exp(EGSSreml$par[1]),exp(EGSSreml$par[2]))
  
  ssq.reml=params.reml[1];   	#  These are the REML estimates.
  tsq.reml=params.reml[2];   	#           --
  
  vx=matrix(0,qp1,qp1);
  for (ti in 1:q)
  {
    vx[(ti+1):qp1,(ti+1):qp1]=matrix(1,1,(qp1-ti))*T.t[ti+1];
  }
  Sigma.mat=ssq.reml*vx;
  Itausq=matrix(rep(0,(qp1*qp1)),nrow=q+1,ncol=q+1);
  diag(Itausq)=rep(tsq.reml,q+1);
  V=Sigma.mat+Itausq;
  D1mat=cbind(-diag(1/S.t),matrix(0,q,1))+cbind(matrix(0,q,1),diag(1/S.t));
  V1mat=D1mat%*%V%*%t(D1mat);
  W.t=(Y.t[2:qp1]-Y.t[1:q])/S.t;
  j1=matrix(1,q,1);
  V1inv=ginv(V1mat);
  mu.reml=(t(j1)%*%V1inv%*%W.t)/(t(j1)%*%V1inv%*%j1);
  j=matrix(1,qp1,1);
  Vinv=ginv(V);
  x0.reml=(t(j)%*%Vinv%*%(Y.t-mu.reml*T.t))/(t(j)%*%Vinv%*%j);
  Var_mu.reml=1/(t(j1)%*%V1inv%*%j1);         #  Variance of mu
  mu_hi.reml=mu.reml+1.96*sqrt(Var_mu.reml);  #  95% CI for mu
  mu_lo.reml=mu.reml-1.96*sqrt(Var_mu.reml);  #       --
  
  #  Calculate estimated population sizes for EGSS model
  #    with Kalman filter, for plotting.
  #
  #  Choose ML or REML estimates here for calculating model values
  #  for plotting (by commenting out the unwanted, default is REML).
  #  mu=mu.ml;  ssq=ssq.ml;  tsq=tsq.ml;  x0=x0.ml;
  mu=mu.reml;  ssq=ssq.reml;  tsq=tsq.reml;  x0=x0.reml;
  
  m[1]=x0;       	#  Initial mean of Y(t).
  v[1]=tsq;      	#  Initial variance of Y(t).
  
  for (ti in 1:q)   #  Loop to generate estimated population abundances
  {                 #    using Kalman filter (see equations 6 & 7,
    #    Dennis et al. (2006)).
    m[ti+1]=mu+(m[ti]+((v[ti]-tsq)/v[ti])*(Y.t[ti]-m[ti]));
    v[ti+1]=tsq*((v[ti]-tsq)/v[ti])+ssq+tsq;
  }
  
  #  The following statement calculates exp{E[X(t) | Y(t), Y(t-1),...,Y(0)]};
  #    see equation 54 in Dennis et al. (2006).  
  Predict.t=exp(m+((v-tsq)/v)*(Y.t-m));
  
  #  Plot the data & model-fitted values
  plot(Observed.t ~ Time.t,xlab="time",ylab="population abundance",
       type="o",lty="solid",pch=1,cex=1);
  #  Population data are circles.
  points(Predict.t ~ Time.t,type="l",lty="dashed",lwd=1);
  #  Estimated abundances are dashed line.
  legend("top", c("Observed.t","Predict.t"),lty=c(1,2),pch=c("o",""),bty="n")
  #  Graph legend
  
  #  Print the parameter estimates
  parms.egoe=c(mu.egoe,ssq.egoe,tsq.egoe,x0.egoe); #  Collect for printing
  parms.egpn=c(mu.egpn,ssq.egpn,tsq.egpn,x0.egpn); #          --
  parms.reml=c(mu.reml,ssq.reml,tsq.reml,x0.reml); #          --
  parms.ml=c(mu.ml,ssq.ml,tsq.ml,x0.ml);           #          --
  names=c("mu","ssq","tsq","x0");                  #          --
  types=c("EGOE","EGPN","EGSS-ML","EGSS-REML");    #          --
  
  #  Print stuff
  matrix(cbind(parms.egoe,parms.egpn,parms.ml,parms.reml),
         nrow=4,ncol=4,byrow=TRUE,dimnames=list(types,names));	
  
  # Add values to Table 3 for the manuscript #
  Tab2_EGSS_REML[which(Tab2_EGSS_REML$Survey == i), c(2:4)] <- c(mu.reml,
                                                                 mu_lo.reml,
                                                                 mu_hi.reml)

  # Add values to a table to use in meta-analysis
  muDF_EGSS[which(muDF_EGSS$survey == i), 2:4] <- c(mu.reml, 
                                                    sqrt(ssq.reml),
                                                    length(my_data[,1]))
   
}

muDF_EGSS$lam <- exp(muDF_EGSS$mu)
muDF_EGSS$lam_lo <- exp(Tab2_EGSS_REML$mu_lo.reml[1:8])
muDF_EGSS$lam_hi <- exp(Tab2_EGSS_REML$mu_hi.reml[1:8])



## 2.2 Quasi-Extinction Thresholds - Surveys ####
# Calculate the mean value for the past 10 years:

Ne_vals <- dat %>% dplyr::select(Year, 
                                 c("mwt.0", "mwt.1", "mwt.2", "ot.0",  "ot.1",  
                                   "ot.2",  "Total", "index")) %>% 
  filter(Year > 2008) %>% 
  summarize_all(mean, na.rm = TRUE) %>% 
  dplyr::select(c("mwt.0", "mwt.1", "mwt.2", "ot.0",  "ot.1",  "ot.2",  
                  "Total", "index"))
# Set a scaling parameter for a fraction of the mean:
Ne_scale <- 0.01 
Ne_vals <- Ne_vals*Ne_scale
# make the data frame tall for convenience:
Ne_vals <- t(Ne_vals)
Ne_vals <- round(Ne_vals, 0)
# set the minimum value to 1:
Ne_vals[which(Ne_vals < 1)] <- 1

# lamDF$Ne <- Ne_vals[,1]
muDF_EGSS$Ne <- Ne_vals[,1]
# set starting size
#  here: the most recent non-NA value:
#  setting it to the lowest value didn't work because some has 0s 
#  in the time series
muDF_EGSS$Nc <- apply(dat[,surveyList], 
                  2, 
                  FUN = function(x) tail(na.omit(x), 1))
muDF_EGSS$Nc <- round(muDF_EGSS$Nc, 0)

# Investigate alternative Ne values - reviewer suggested using 2005:2015 #
# We didn't end up using these, but they are similar to what we did use
Ne_vals_alt <- dat %>% dplyr::select(Year, 
                                     c("mwt.0", "mwt.1", "mwt.2", "ot.0",  
                                       "ot.1",  "ot.2",  "Total", "index")) %>% 
  filter(Year > 2004) %>% 
  filter(Year < 2016) %>% 
  summarize_all(mean, na.rm = TRUE) %>% 
  dplyr::select(c("mwt.0", "mwt.1", "mwt.2", "ot.0",  "ot.1",  "ot.2",  "Total", 
                  "index"))
round(Ne_vals_alt * Ne_scale, 0)
mean(as.numeric(round(Ne_vals_alt * Ne_scale, 0)))

## 2.3 Cumulative Extinction Probabilities - surveys ####
# using countCDFxt to get CDF with CIs 
#    based on suggestion from SSA peer review

for(i in 1:length(surveyList)){
  r = which(muDF_EGSS$survey == surveyList[i])
  assign(paste0("lam_", surveyList[i], "_QEP_EGSS"), 
         countCDFxt(mu = muDF_EGSS$mu[r], 
                    sig2 = (muDF_EGSS$sd[r])^2, 
                    Nc = muDF_EGSS$Nc[i],
                    Ne = muDF_EGSS$Ne[i], 
                    nt = muDF_EGSS$n[r],
                    tmax = 50,
                    plot = FALSE))
}



# 3. Check Assumptions ####

## 3.1 Calculate year-specific lambda values ####
#     & the mean using the diffusion approximation 

# Note: we used the diffusion approximation to calculate the mean 
#    lambda values in an earlier version of the manuscript. 
#    The following code was originally developed for that analysis. 
#    It also calculates the year-specific values so we're including 
#    it here to assist with assumption checking.


for(i in surveyList){
  assign(paste0("lam_", i), 
         calcLam2(dataset = dat, 
                  abCol = i, 
                  yearCol = "Year"))
}

lamDF <- data.frame(matrix(unlist(lapply(paste0("lam_", surveyList), 
                                         function(x) data.frame(get(x)[9:11]))),
                           ncol = 3, byrow = TRUE))
names(lamDF) <- c("lambda", "lambda_lcl", "lambda_ucl")
lamDF$survey <- surveyList

lamDF$Ne <- Ne_vals[,1]
lamDF$Nc <- apply(dat[,surveyList], 
                  2, 
                  FUN = function(x) tail(na.omit(x), 1))

## 3.2 No change in growth rate over time ####
# Graph Lambda over time 
png(paste0(root_figs, "Assumptions_LambdaTime.png"),
    width = 8.5, height = 11, units = "in", res = 300)
par(mar = c(4.1, 4.1, 2.1, 1.1),
    cex = 1.25,
    mfrow = c(4, 2))
for(i in 2:length(names(dat))){
  datplot <- eval(parse(text = paste("lam", 
                                     names(dat)[i],
                                     sep = "_")))
  datplot <- datplot$datDF[, c("year", "loglam")]
  datplot$loglam[which(is.infinite(datplot$loglam))] <- NA
  dlm <- lm(datplot$loglam ~ datplot$year)
  plot(datplot$year, 
       datplot$loglam,
       pch = 16,
       xlim = c(1967, 2021),
       ylim = c(-6, 5),
       xlab = "Year",
       ylab = "log(Lambda)",
       main = names(dat)[i])
  text(x = 1975, y =4, 
       labels = paste("R2 = ", round(summary(dlm)$r.squared, 2),
                      ", p = ", round(summary(dlm)$coefficients[2, 4], 2),
                      sep = ""), 
       cex = 0.85)
}
dev.off()

## 3.3 Density Dependence ####
# see Morris & Doak p. 92
# see p. 94 for quantitative methods for checking assumptions
# lag() numbers are for the previous year = t-1, or t in M&D
# unlagged values = year t, or t+1 in M&D

png(paste0(root_figs, "Assumptions_DDep.png"),
    width = 8.5, height = 11, units = "in", res = 300)
par(mar = c(4.1, 4.1, 2.1, 1.1),
    cex = 1.25,
    mfrow = c(4, 2))
for(i in 2:length(names(dat))){
  x <- lag(dat[, i])
  y <- log(dat[, i]/lag(dat[, i]))
  y[which(is.infinite(y))] <- NA
  dlm <- lm(y ~ x)
  plot(x,
       y, 
       pch = 16,
       xlab = paste(names(dat)[i], "previous year (Nt)"),
       ylab = "log(Nt+1/Nt)",
       main = names(dat)[i])
  text(x = 0.8*max(x, na.rm = TRUE), 
       y = 0.8*max(y, na.rm = TRUE), 
       labels = paste("R2 = ", round(summary(dlm)$r.squared, 2), 
                      ", Slope = ", round(summary(dlm)$coefficients[2, 1], 4), 
                      ", p = ", round(summary(dlm)$coefficients[2, 4], 2),
                      sep = ""), 
       cex = 0.85)
  
}
dev.off()


# "linear" models = no density dependence
ret <- data.frame("Survey" = NA,
                  "DensityDep" = NA, 
                  "df" = NA, 
                  "AICc" = NA, 
                  "Converged" = NA,
                  "r"   = NA,
                  "R_p" = NA,
                  "K"   = NA,
                  "K_p" = NA)
for(i in 2:length(names(dat))){
  print(names(dat)[i])
  df <- eval(parse(text = paste("lam_",
                                names(dat)[i],
                                "$datDF",
                                sep = "")))
  nt <- df$ab        #N_t
  nt_1 <- lag(nt)    #N_tminus1
  loglam <- df$loglam
  loglam[is.infinite(loglam)] <- NA
  
  tryCatch({
    # Beverton Holt
    regBH <- nls(loglam ~ R/(1 + (nt_1*R/K)),
                 start = list(R =1,
                              K = 1000),
                 control = list(warnOnly = TRUE))
    # Ricker
    regR <- nls(loglam ~ R*(1 - nt_1/K),
                start = list(R = 1,
                             K = 8000),
                control = list(warnOnly = TRUE))
    # Density Independent
    regL <- glm(loglam ~ nt_1)    
    # compile results
    
  },
  error = function(e) print(e),
  finally = {
    ret2 <- cbind(
      rep(names(dat)[i], 3),
      c("Beverton-Holt", "Ricker", "Linear"),
      AICc(regBH, regR, regL),
      c(regBH$convInfo$isConv,
        regR$convInfo$isConv,
        regL$converged),
      c(coefficients(regBH)[1], #R
        if (regR$convInfo$stopMessage == "singular gradient") NA else coefficients(regR)[1],
        coefficients(regL)[2]),
      c(summary(regBH)$coefficients[1,4], #p-value for R
        if (regR$convInfo$stopMessage == "singular gradient") NA else summary(regR)$coefficients[1,4],
        summary(regL)$coefficients[2,4]),
      c(coefficients(regBH)[2], #K
        if (regR$convInfo$stopMessage == "singular gradient") NA else coefficients(regR)[2],
        NA),
      c(summary(regBH)$coefficients[2,4], #p-value for k
        if (regR$convInfo$stopMessage == "singular gradient") NA else summary(regR)$coefficients[2,4],
        NA))
    names(ret2) <- c("Survey", "DensityDep", "df", 
                     "AICc", "Converged", "r",  
                     "R_p",
                     "K", 
                     "K_p")
    ret <- rbind(ret, ret2)
    rm(df, nt, nt_1, regBH, regR, regL, ret2)
    next
  } 
  )
}
ret <- ret[-1,]
ret 


## 3.4 Check for autocorrelation ####
png(paste0(root_figs, "Assumptions_Autocorr.png"),
    width = 8.5, height = 11, units = "in", res = 300)
par(mar = c(4.1, 4.1, 2.1, 1.1),
    cex = 1.25,
    mfrow = c(4, 2))
for(i in 2:length(names(dat))){
  datplot <- eval(parse(text = paste("lam", 
                                     names(dat)[i],
                                     sep = "_")))
  datplot2 <- datplot$datDF[, c("year", "loglam")]
  datplot2$loglam[which(is.infinite(datplot2$loglam))] <- NA
  plot(
    datplot$reg$residuals,
    pch = 16,
    xlab = "Year",
    ylab = "Residuals",
    main = names(dat)[i])
  abline(h = 0, col = "grey", lty = 2)
  text(x = 5, y = 0.9*max(datplot$reg$residuals), 
       labels = paste("p =", 
                      round(dwtest(datplot$reg, alternative = "two.sided")$p.value, 2)), 
       cex = 0.85)
}
dev.off()


for(i in 2:length(names(dat))){
  datplot <- eval(parse(text = paste("lam", 
                                     names(dat)[i],
                                     sep = "_")))
  datplot2 <- datplot$datDF[, c("year", "loglam")]
  datplot2$loglam[which(is.infinite(datplot2$loglam))] <- NA
  
  plot(lag(datplot$reg$residuals),
       datplot$reg$residuals)
  
  print(names(dat)[i])
  print("Slope:")
  print(summary(lm(datplot$reg$residuals ~ lag(datplot$reg$residuals)))$coefficients[2, 4])
  print("DW Test:")
  print(dwtest(datplot$reg, alternative = "two.sided")$p.value)
  
}

# ACF
acf.vec <- numeric(length(names(dat))-1)
for(i in 2:length(names(dat))){
  datplot <- eval(parse(text = paste("lam", 
                                     names(dat)[i],
                                     sep = "_")))
  datplot2 <- datplot$datDF[, c("year", "loglam")]
  datplot2$loglam[which(is.infinite(datplot2$loglam))] <- NA
  
  print(names(dat)[i])
  acf.vec[i-1] <- pacf(datplot$reg$residuals, main = names(dat)[i])[1]$acf[1,1,1] #corr coeffs
}

acf.vec <- unlist(acf.vec)
# Not sure where this came from: 
# acf.vec[c(1, 8)] <- 0 #FMWT have doesn't significant autocorr
acf.vec[c(4, 9)] <- 0 #mwt.2 & 20-mm don't have significant autocorr
m.corr <- mean(acf.vec) #didn't change much

for(i in 2:length(names(dat))){
  datplot <- eval(parse(text = paste("lam", 
                                     names(dat)[i],
                                     sep = "_")))
  datplot2 <- datplot$datDF[, c("year", "loglam")]
  datplot2$loglam[which(is.infinite(datplot2$loglam))] <- NA
  
  print(pacf(datplot$reg$residuals, main = names(dat)[i])[1])
  
  print(names(dat)[i])
  # p < 0.05 ==> autocorr
  print(dwtest(datplot$reg, alternative = "greater"))
  print(dwtest(datplot$reg, alternative = "two.sided"))
  print(dwtest(datplot$reg, alternative = "less"))
  print("-----")
}


# 4. Meta-Analysis ####
## 4.1 Calculate Meta-analysis mean and CIs ####
metaMu_EGSS <- metamean(n = muDF_EGSS$n,
                        mean = muDF_EGSS$mu,
                        sd = muDF_EGSS$sd,
                        studlab = surveyList)
summary(metaMu_EGSS)
forest(metaMu_EGSS)

# Save estimates to the summary table
Tab2_EGSS_REML[9, 1] <- "meta-analysis mean"
Tab2_EGSS_REML[9, 2:4] <- c(metaMu_EGSS$TE.random,
                            metaMu_EGSS$lower.random,
                            metaMu_EGSS$upper.random)

muDF_EGSS[9, 1] <- "meta-analysis mean"
muDF_EGSS[9, c(2, 5:7)] <- c(metaMu_EGSS$TE.random,
                             exp(metaMu_EGSS$TE.random),
                             exp(metaMu_EGSS$lower.random),
                             exp(metaMu_EGSS$upper.random))
mean(muDF_EGSS$Ne[1:8]/muDF_EGSS$Nc[1:8])
mean(muDF_EGSS$Ne[1:7]/muDF_EGSS$Nc[1:7])
muDF_EGSS[9, 8] <- 50 #Ne value for meta-analysis CDF = 2.6% of Nc = 1854*0.026
muDF_EGSS[9, 9] <- round(mean(muDF_EGSS$Nc[1:8]), 0)



## 4.2 Quasi-Extinction Probabilities - Meta-Analysis ####
# Bootstrapped QE Prob from Random model
# using mean "current abundance" and mean QE threshold:
meanCDF_EGSS <- countCDFxt(mu = metaMu_EGSS$TE.random, 
                           sig2 = metaMu_EGSS$seTE.fixed*sqrt(sum(metaMu_EGSS$n)), 
                           Nc = muDF_EGSS$Nc[9],  
                           Ne = muDF_EGSS$Ne[9], 
                           nt = 40,
                           tmax = 30,
                           plot = FALSE)
test_EGSS <- extCDFAC(mu = metaMu_EGSS$TE.random ,
                      sig2 = metaMu_EGSS$seTE.fixed*sqrt(sum(metaMu_EGSS$n)),
                      Nc = muDF_EGSS$Nc[9],
                      Ne = muDF_EGSS$Ne[9], 
                      dd = "none",
                      rho = m.corr, #correlation coefficient
                      numReps = 10000,
                      tmax = 50)
test.g_EGSS <- rsq2(test_EGSS, n = 100, B = 5000)


## 4.3 Consider stable population growth ####
# like in Schultz et al --> What if mu = 0?
# These are the means and CIs in the text of the results section

testStable <- extCDFAC(mu = 0 ,
                       sig2 = metaMu_EGSS$seTE.fixed*sqrt(sum(metaMu_EGSS$n)), #matches original inputs to meta-analysis cdf
                       Nc = muDF_EGSS$Nc[9],  
                       Ne = muDF_EGSS$Ne[9],
                       dd = "none",
                       rho = m.corr, #correlation coefficient
                       numReps = 10000,
                       tmax = 50)
testStable.g <- rsq2(testStable, n = 100, B = 5000)
testStable.g





# 5. Percent Decline ####
decline <- data.frame(Survey = Tab2_EGSS_REML$Survey, # Names of surveys
                      # Percent decline over 15 years:
                      dec_10y = round(100 *(1 - (exp(Tab2_EGSS_REML$mu.reml) ^ 10)), 1),
                      # Percent decline over 15 years:
                      dec_15y = round(100 *(1 - (exp(Tab2_EGSS_REML$mu.reml) ^ 15)), 1),
                      # Percent decline over 20 years:
                      dec_20y = round(100 *(1 - (exp(Tab2_EGSS_REML$mu.reml) ^ 20)), 1))



# 6. Manuscript Figures ####
## Fig. 2 Meta-analysis plot ####
png(paste0(root_figs, "Fig2_EGSSmeta.png"),
    width = 8, height = 5, units = "in", res = 300)
par(mar = c(5.1, 10.5, 4.1, 2.1), #5.1 4.1 4.1 2.1
    cex.axis = 1.25, cex.lab = 1.25)
plotCI(exp(c(metaMu_EGSS$mean, metaMu_EGSS$TE.random)), 
       1:(length(surveyList)+1), 
       li = exp(c(metaMu_EGSS$lower, metaMu_EGSS$lower.random)),
       ui = exp(c(metaMu_EGSS$upper, metaMu_EGSS$upper.random)),
       pch = 16, lwd = 2,
       err = "x",
       xlab = "Annual Population Growth Rate",
       ylab = "",
       yaxt = "n",
       xlim = c(0.5, 1.7),
       cex = 1.5,
       bty = "n")
axis(side = 2, 
     at = 1:(length(surveyList)+1), 
     las = 2,
     labels =
       c( "SFBS MWT Age-0",
          "SFBS MWT Age-1",
          "SFBS MWT Age-2",
          "SFBS OT Age-0",
          "SFBS OT Age-1",
          "SFBS OT Age-2",
          "FMWT",
          "20-mm",
          "Meta-Analysis"),
     tick = FALSE,
     pos = 0.6)
abline(v = 1, lty = 2)
dev.off()



## Fig. 3 Composite plot of extinction probabilities ####
# Figure 3 in manuscript now

# create a crosswalk for names of surveys so the labels stay in the right order
surveyCross <- data.frame(short = Tab2_EGSS_REML$Survey,
                          long = c("SFBS MWT Age-0",
                                   "SFBS MWT Age-1",
                                   "SFBS MWT Age-2",
                                   "SFBS OT Age-0",
                                   "SFBS OT Age-1",
                                   "SFBS OT Age-2",
                                   "FMWT",
                                   "20-mm",
                                   "Meta-Analysis"))

# plotting parameters in the order of surveyList
ltyList <- c(1, 2, 3,
             1, 2, 3,
             1, 2)
colList <- c("black", "black", "black",          #original sfbs MWT
             "grey", "grey", "grey",             #original sfbs OT
             "cyan4", "cyan4")                   #FMWT, 20-mm
legList <- NA
for (i in 1:length(surveyList)){
  legList[i] <- surveyCross$long[which(surveyCross$short == surveyList[i])]
}


legList2 <- NA
for (i in 1:length(Tab2_EGSS_REML$Survey)){
  legList2[i] <- surveyCross$long[which(surveyCross$short == Tab2_EGSS_REML$Survey[i])]
}

png(paste0(root_figs, "Fig3_CompCDF_EGSS.png"),
    width = 8, height = 5, units = "in", res = 300)
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(0, 0, type = "n",
     xlim = c(1, 50),
     ylim = c(0, 1),
     xlab = "",
     ylab = "Probability of Quasi-Exinction",
     xaxt = "n",
     xaxs = "i",
     yaxs = "i")
axis(side = 1, at = seq(2, 120, by = 5), labels = seq(2020, (2020+119), by = 5))
# reviewer suggested horizontal line at 0.2
abline(h = 0.2, lwd = 4, col = "honeydew3")
for(i in 1:length(surveyList)){
  lines(1:50, get(paste0("lam_", surveyList[i], "_QEP_EGSS"))$Gbest, 
        lty = ltyList[i],
        col = colList[i],
        lwd = 2)
}
legend("topleft", 
       # reverse the order so it matches the order on the lambda summary plot
       lty = rev(ltyList),
       col = rev(colList),
       lwd = 2,
       # use the survey name crosswalk
       legend = rev(legList2[-9]),
       title = "Survey",
       bty = "n",
       cex = 0.85)
# Add reference lines for IUCN Redlist Criteria E:
# segments(x0 = c(0, 0, 0), x1 = c(10, 20, 100), y0 = c(0.5, 0.2, 0.1))
# segments(x0 = c(10, 20, 100), y0 = c(0, 0, 0), y1 = c(0.5, 0.2, 0.1))
dev.off()


## Fig. 4 Meta Analysis CDF with bootstrapped CIs ####
png(paste0(root_figs, "Fig4_Meta_CDF_boot_EGSS.png"),
    width = 8, height = 5, units = "in", res = 300)
par(mar = c(5.1, 4.1, 4.1, 2.1), #5.1 4.1 4.1 2.1
    mfrow = c(1, 1),
    cex.axis = 1.25, cex.lab = 1.25)
plot(test.g_EGSS$mean, 
     type = "l", 
     lwd = 2,
     ylim = c(0.0001, 1),
     xaxt = "n",
     xlab = "",
     ylab = "Probability of Quasi-Exinction",
     #log = "y",
     xaxs = "i",
     yaxs = "i")
# add the horizontal line at 20% suggested by reviewer
abline(h = 0.2, lwd = 4, col = "honeydew3")
# plot the mean line back over the top of the 20% line
lines(test.g_EGSS$mean, lwd = 2)
lines(test.g_EGSS$lower, lty = 2, lwd = 2)
lines(test.g_EGSS$upper, lty = 2, lwd = 2)
axis(side = 1, at = seq(1, 120, by = 5), labels = seq(2020, (2020+119), by = 5))
dev.off()



## Fig. 5 Sensitivity Graphs ####
# color palette for graphs:
colFun <- colorRampPalette(c("#cccccc", "#252525"))

# Vary starting abundance (Nc) 
png(paste0(root_figs, "Fig5_Sens_Nc_Ne_EGSS.png"),
    width = 8, height = 10, units = "in", res = 300)
par(cex = 1.25, 
    mfrow = c(2, 1))
plot(0, 0, type = "n",
     ylim = c(0, 1),
     xlim = c(1, 30),
     xaxt = "n",
     main = "Sensitivity to Starting Abundance Value",
     xlab = "Year",
     ylab = "Probability of Quasi-Extinction")
NcSens <- seq(1, 3500, length.out = 50)
for(i in 1:50){ 
  yvals <- countCDFxt(mu = metaMu_EGSS$TE.random, 
                      sig2 = metaMu_EGSS$seTE.fixed*sqrt(sum(metaMu_EGSS$n)), 
                      Nc = NcSens[i], 
                      Ne = round(mean(Ne_vals), 0), 
                      nt = 40,
                      tmax = 30,
                      plot = FALSE)
  lines(1:30,
        yvals$Gbest,
        col = colFun(50)[i])
}
# add a line to highlight a particular value
# lines(1:30,
#       countCDFxt(mu = metaMu_EGSS$TE.random,
#                  sig2 = metaMu_EGSS$seTE.fixed*sqrt(sum(metaMu_EGSS$n)),
#                  Nc = round(mean(lamDF$Nc), 0),
#                  Ne = round(mean(Ne_vals), 0),
#                  nt = 40,
#                  tmax = 30,
#                  plot = FALSE)$Gbest,
#       lwd = 3, col = "red")
axis(side = 1, at = seq(1, 35, by = 5), labels = seq(2020, 2050, by = 5))
legend("top",
       title = "Color Ramp for Starting Values",
       horiz = TRUE,
       lwd = 2,
       bty = "n", cex = 0.8,
       col = colFun(50)[c(1, 13, 25, 38, 50)],
       legend = round(seq(1, 3500, length.out = 50)[c(1, 13, 25, 38, 50)], 0))

# Vary QE Thresholds (Ne)
plot(0, 0, type = "n",
     ylim = c(0, 1),
     xlim = c(1, 30),
     xaxt = "n",
     main = "Sensitivity to Quasi-Extinction Threshold Value",
     xlab = "Year",
     ylab = "Probability of Quasi-Extinction")
for(i in 1:50){
  yvals <- countCDFxt(mu = metaMu_EGSS$TE.random, 
                      sig2 = metaMu_EGSS$seTE.fixed*sqrt(sum(metaMu_EGSS$n)),
                      Nc = round(mean(lamDF$Nc), 0),
                      Ne = i,
                      nt = 40,
                      tmax = 30,
                      plot = FALSE)
  lines(1:30,
        yvals$Gbest,
        col = colFun(50)[i])
}
# add a line to highlight a particular value
# lines(1:30,
#       countCDFxt(mu = metaMu_EGSS$TE.random, 
#                  sig2 = metaMu_EGSS$seTE.fixed*sqrt(sum(metaMu_EGSS$n)),
#                  Nc = round(mean(lamDF$Nc), 0), 
#                  Ne = 50, #round(mean(Ne_vals), 0), 
#                  nt = 40,
#                  tmax = 30,
#                  plot = FALSE)$Gbest,
#       lwd = 3, col = "red")
axis(side = 1, at = seq(1, 35, by = 5), labels = seq(2020, 2050, by = 5))
legend("top",
       horiz = TRUE,
       title = "Color Ramp for Extinction Thresholds",
       lwd = 2,
       bty = "n", cex = 0.8,
       col = colFun(50)[c(1, 13, 25, 38, 50)],
       legend = c(1, 13, 25, 38, 50))
dev.off()
par(mfrow = c(1, 1))

# 7. Mann-Kendall Trends ####
trend.table <- data.frame(survey = names(dat)[-1],
                          z = NA,  # statistic
                          p = NA,  # p.value
                          n = NA,  # parameter
                          S = NA,  # estimates[1]
                          varS = NA, # estimates[2]
                          tau = NA)  # estimates[3]
for(i in 2:ncol(dat)){
  
  # label the output with the colum name
  print(names(dat)[i])
  
  # make a vector out of the column of interest
  vect <- dat[,i]
  
  # trim the NAs off the ends of the vector
  vect <- vect[min(which(!is.na(vect))):max(which(!is.na(vect)))]
  nopes <- which(is.na(vect))
  
  # replace the internal NAs with mean of two neighboring values
  
  for(j in nopes){
    x1 <- if(is.na(vect[j-1])) vect[j-2] else vect[j-1]
    x2 <- if(is.na(vect[j+1])) vect[j+2] else vect[j+1]
    vect[j] <- mean(x1, x2)
  }
  
  # run a Mann-Kendall test for a trend
  # print(mk.test(vect))
  mk <- mk.test(vect)
  
  # put a place marker at the end for readability
  # print ("------")
  
  # put values into the table on the row that corresponds to the survey
  trend.table[i-1, 2:7] <- c(mk[[3]], #"statistic",
                             mk[[2]], #"p.value",
                             mk[[5]], #"parameter",
                             mk[[6]])   #, #"estimates")
}
# clean up format for including in papers:
trend.table$z <- round(trend.table$z, 2)
trend.table$p <- round(trend.table$p, 2)
trend.table$varS <- round(trend.table$varS, 2)
trend.table$tau <- round(trend.table$tau, 2)