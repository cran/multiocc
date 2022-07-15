## We assume the raw data are collected in space and time, for S species.
## S = number of species
## T = number of unique time periods
## L = number of unique locations

## Dependencies
#library(MASS)
#library(truncnorm)
#library(tmvtnorm)
#library(fields)
#library(akima)
#library(MCMCpack)

#' This function runs the MCMC.
#'
#' @param M.iter The total number of iterations in MCMC
#' @param M.burn The length of the burn in
#' @param M.thin The number to thin the chain.  Thinning by 10 only stores every 10th run.
#' @param model.input A list of output created by running the create.data.R function
#' @param q Desired number of Moran's I basis functions in the restricted spatial regression model
#' @param sv A TRUE/FALSE on whether or not the MCMC output should be saved as 'MCMC.Rdata' and overwritten every 1000 iterations.  Defaults to false.
#' @return A list with all standard MCMC output
#' @export
#' @examples
#'
#' data(Example)
#' head(detection)
#' head(occupancy)
#' head(coords)
#' names = list("species"=colnames(detection)[4:9],
#' "detection"=c("duration"),"occupancy"=c("forest","elev"))
#' model.input = create.data(detection, occupancy, coords, names, threshold = 15000)
#' out = run.mcmc(M.iter=3, M.burn=1, M.thin=1, model.input, q=10, sv=FALSE)

run.mcmc = function(M.iter, M.burn=NULL, M.thin=NULL, model.input, q=NULL, sv=FALSE){

  if(q>min(summary(as.factor(model.input$occupancy.info$season)))){
    stop("Number of basis functions cannot be larger than the smallest number of sites in any season")
  }

  if(!is.null(q)){
    basis = make.basis(q,model.input)
  }else{
    q = ceiling(.1*min(summary(as.factor(model.input$occupancy.info$season))))
    basis = make.basis(q,model.input)
  }

  if(is.null(M.thin)){
    M.thin=1
  }

  if(is.null(M.burn)){
    M.burn= floor(M.iter/2)
  }

  y = as.matrix(model.input$y)
  seasonInd = model.input$detection.info$season
  surveyInd = model.input$detection.info$survey
  site = model.input$detection.info$siteID
  X = as.matrix(model.input$X)
  W = as.matrix(model.input$W)

  ## Indicator if a site/season/survey combo led to observations
  onI = model.input$detection.info$observations

  #onI.temp = aggregate(x = model.input$detection.info$visited,
  #          by = list(model.input$detection.info$siteID, model.input$detection.info$season),
  #          FUN = sum)

  #onI.occupancy = 1*(onI.temp$x>0) ### indicates if we have an observation or not
  ### note if we do not have an observation, our model simulates from the posterior predictive distribution for z


  KtQK=basis[[1]]
  q=dim(KtQK)[1]
  K=basis[[2]]
  ### note seasonInd is about y which may be longer than z...need another seasonInd for z
  #zseasonInd = seasonInd[which(surveyInd==1)]
  #zsite = site[which(surveyInd==1)]
  zseasonInd = model.input$occupancy.info$season
  zsite = model.input$occupancy.info$siteID

  ### number of seasons
  T = max(zseasonInd)
  nT = nrow(X)

  ### number of site visits for each location in each season
  ### for each location each season, get an index of all surveys - all rows in y that correspond to each site/season
  #Ktot = rep(0,nT)
  #Kind = list()
  #for(j in 1:nT){
  #  Kind[[j]] = which(seasonInd==zseasonInd[j] & site==zsite[j])
  #  Ktot[j] = length(Kind[[j]])
  #}

  S = dim(y)[2]
  # y = species ### column for each species;
  ### make y=0 where we do not have an observation. These aren't used for the posterior since we subset them out,
  ### but just using to help with knowing where to predict z
  y[onI==0,]=rep(0,S)
  ysum = matrix(0,nT,S)
  for(i in 1:nT){
    II = which(site==zsite[i] & seasonInd == zseasonInd[i])
    ysum[i,] = colSums(matrix(y[II,],length(II),S))
  }
  #ysum = as.matrix(aggregate(y,by=list(site,seasonInd),FUN=sum)[,3:5]) ### number of detections in the Ktot site visits
  z = 1*(ysum>0) ### if y>0 for any site visit then z=1  ap = dim(X)[2]

  ap = dim(X)[2]
  bp = dim(W)[2]

  alpha = matrix(0,ap,S)
  beta = matrix(0,bp,S)
  Sig = diag(S)
  rho.num=0.9
  rho = rho.num*diag(S)

  ### gg = gamma
  gg = matrix(0,q*T,S)
  gvec = as.vector(gg)
  ztilde = z-0.5 ### need initial value negative if z<0 and positive if z>0.
  ytilde = y-0.5

  ### hyperparameters
  Signu = S+2
  SigPsi = diag(S)

  Zout = array(0,c(dim(z)[1],(M.iter-M.burn)/M.thin,S))
  betaout = array(0,c(bp,(M.iter-M.burn)/M.thin,S))
  alphaout = array(0,c(ap,(M.iter-M.burn)/M.thin,S))
  gamout = array(0,c(q*T,(M.iter-M.burn)/M.thin,S))
  sigout = array(0,c(S^2,(M.iter-M.burn)/M.thin))
  rhoout = array(0,c(S,(M.iter-M.burn)/M.thin))
  betaburn = array(0,c(bp,M.burn/M.thin,S))
  alphaburn = array(0,c(ap,M.burn/M.thin,S))
  sigburn = array(0,c(S^2,M.burn/M.thin))
  gamburn = array(0,c(q*T,(M.burn)/M.thin,S))
  psiout = array(0,c(nrow(X),(M.iter-M.burn)/M.thin,S))
  pout = array(0,c(nrow(W),(M.iter-M.burn)/M.thin,S))
  loglik = array(0,c(nrow(W),(M.iter-M.burn)/M.thin,S))

  #### Note we only need to update the z's where the corresponding y=0
  #### If y=1, then z has to be 1. This is already true from the initial
  #### condition, so those z's never need to be updated...

  nyz = dim(ysum)[1]-colSums(1*(ysum)>0) ### how many of the y's are not 1
  zInd = list()
  for(s in 1:S){
    zInd[[s]] = which(ysum[,s]==0)
  }

  st = Sys.time()

  for(m in 1:M.iter){

    ####### update the z's where the corresponding y are 0 (or missing) for each site visit
    ## for each species
    for (s in 1:S){
      probpsi0 = X%*%alpha[,s]+K%*%gg[,s]
      psi0 = pnorm(probpsi0)
      probp0 = W%*%beta[,s]
      p0 = pnorm(probp0)
      ### need the product of the 1-p0 for each site/season combo
      temp  = p0*onI
      temp[onI==0]=0
      p0cit = aggregate(1-temp,by=list(site,seasonInd),FUN=prod)[,3]
      psibar0 = (psi0*p0cit)/(1-psi0+psi0*p0cit)
      z[zInd[[s]],s] = (log(runif(nyz[s]))<log(psibar0[zInd[[s]]]))
    }

    ######## update each ztilde
    for (s in 1:S){
      Ind0 = which(z[,s]==0)
      Ind1 = which(z[,s]==1)
      r0=rtruncnorm(n=length(Ind0), a=-Inf, b=0, mean = X[Ind0,]%*%alpha[,s]+K[Ind0,]%*%gg[,s], sd = 1)
      ztilde[Ind0,s] = r0
      r1=rtruncnorm(n=length(Ind1), a=0, b=Inf, mean = X[Ind1,]%*%alpha[,s]+K[Ind1,]%*%gg[,s], sd = 1)
      ztilde[Ind1,s] = r1
    }

    ######## update each ytilde
    for (s in 1:S){
      Ind0 = which(y[,s]==0)
      Ind1 = which(y[,s]==1)
      r0=rtruncnorm(n=length(Ind0), a=-Inf, b=0, mean = W[Ind0,]%*%beta[,s], sd = 1)
      ytilde[Ind0,s] = r0
      r1=rtruncnorm(n=length(Ind1), a=0, b=Inf, mean = W[Ind1,]%*%beta[,s], sd = 1)
      ytilde[Ind1,s] = r1
    }

    ######## update gvec and gg
    ### t=1
    t=1
    TID = c() ### observations from same time period
    TIDg = c() ### rows of g from same time period
    TIDg1 = c() ### rows of g from next time period
    Tind = which(zseasonInd==t)
    for(s in 1:S){
      TID = c(TID, (s-1)*nT+Tind)
      TIDg = c(TIDg, ((s-1)*q*T+(t-1)*q+1):((s-1)*q*T+t*q))
      TIDg1 = c(TIDg1, ((s-1)*q*T+t*q+1):((s-1)*q*T+(t+1)*q))
    }
    ### compute precision matrix for full conditional of gamma_t
    gp = diag(S*q)+kronecker(solve(Sig),KtQK[,,t])+t(kronecker(rho,diag(q)))%*%kronecker(solve(Sig),KtQK[,,t+1])%*%kronecker(rho,diag(q))
    gcov=chol2inv(chol(gp))
    gmu = gcov%*%(t(kronecker(diag(S),K[Tind,((t-1)*q+1):(t*q)]))%*%as.vector(ztilde[Tind,]-X[Tind,]%*%alpha)+t(kronecker(rho,diag(q)))%*%kronecker(solve(Sig),KtQK[,,(t+1)])%*%gvec[TIDg1])
    gstar = mvrnorm(1,gmu,gcov)
    for (s in 1:S){
      gvec[TIDg[((s-1)*q+1):(s*q)]] =  gstar[((s-1)*q+1):(s*q)]-mean(gstar[((s-1)*q+1):(s*q)])
    }

    if (T>2) {
      ### update for t=2,...,T-1
      for (t in 2:(T-1)){
        TID = c() ### observations from same time period
        TIDg = c() ### rows of g from same time period
        TIDg0 = c() ### rows of g from previous time period
        TIDg1 = c() ### rows of g from next time period
        Tind = which(zseasonInd==t)
        for (s in 1:S){
          TID = c(TID, (s-1)*nT+Tind)
          TIDg=c(TIDg, ((s-1)*q*T+(t-1)*q+1):((s-1)*q*T+t*q))
          TIDg0 = c(TIDg0,((s-1)*q*T+(t-2)*q+1):((s-1)*q*T+(t-1)*q))
          TIDg1=c(TIDg1, ((s-1)*q*T+t*q+1):((s-1)*q*T+(t+1)*q))
        }
        gp = diag(S*q)+kronecker(solve(Sig),KtQK[,,t])+t(kronecker(rho,diag(q)))%*%kronecker(solve(Sig),KtQK[,,t+1])%*%kronecker(rho,diag(q))
        gcov=chol2inv(chol(gp))
        gmu = gcov%*%(t(kronecker(diag(S),K[Tind,((t-1)*q+1):(t*q)]))%*%as.vector(ztilde[Tind,]-X[Tind,]%*%alpha)+t(kronecker(rho,diag(q)))%*%kronecker(solve(Sig),KtQK[,,t+1])%*%gvec[TIDg1]+kronecker(solve(Sig),KtQK[,,t])%*%kronecker(rho,diag(q))%*%gvec[TIDg0])
        gstar = mvrnorm(1,gmu,gcov)
        for (s in 1:S){
          gvec[TIDg[((s-1)*q+1):(s*q)]] =  gstar[((s-1)*q+1):(s*q)]-mean(gstar[((s-1)*q+1):(s*q)])
        }
      }
    }
    if (T>1) {
      #### update for t=T
      t=T
      TID = c() ### observations from same time period
      TIDg = c() ### rows of g from same time period
      TIDg0 = c() ### rows of g from previous time period
      Tind = which(zseasonInd==t)
      for (s in 1:S){
        TID = c(TID, (s-1)*nT+Tind)
        TIDg= c(TIDg, ((s-1)*q*T+(t-1)*q+1):((s-1)*q*T+t*q))
        TIDg0 = c(TIDg0,((s-1)*q*T+(t-2)*q+1):((s-1)*q*T+(t-1)*q))
      }
      gp = diag(S*q)+kronecker(solve(Sig),KtQK[,,t])
      gcov=chol2inv(chol(gp))
      gmu = gcov%*%(t(kronecker(diag(S),K[Tind,((t-1)*q+1):(t*q)]))%*%as.vector(ztilde[Tind,]-X[Tind,]%*%alpha)+kronecker(solve(Sig),KtQK[,,t])%*%kronecker(rho,diag(q))%*%gvec[TIDg0])
      gstar = mvrnorm(1,gmu,gcov)
      for (s in 1:S){
        gvec[TIDg[((s-1)*q+1):(s*q)]] =  gstar[((s-1)*q+1):(s*q)]-mean(gstar[((s-1)*q+1):(s*q)])
      }
    }
    gg = matrix(gvec,q*T,S)

    #### update beta
    ## This requires expanding z
    z.temp = z
    colnames(z.temp) = paste("z",1:dim(z.temp)[2],sep=".")
    hold = data.frame(zsite, zseasonInd, z.temp)
    colnames(hold) = c("siteID", "season",colnames(z.temp))
    expanded.z = merge(model.input$detection.info, hold, by=c("siteID","season"), all.x=TRUE)
    expanded.z = expanded.z[,colnames(z.temp)]

    for (s in 1:S){
      #### pull out the rows of y for which the corresponding z and onI.occupancy is 1

      zstar = which(expanded.z[,s]*onI==1)
      # star = unlist(Kind[zstar])
      bcov = solve(t(W[zstar,])%*%W[zstar,])
      bmu = bcov%*%t(W[zstar,])%*%ytilde[zstar,s]
      beta[,s] = mvrnorm(1,bmu,bcov)
    }

    ### update alpha
    for (s in 1:S){
      acov = solve(t(X)%*%X)
      amu = acov%*%(t(X)%*%(ztilde[,s]-K%*%gg[,s]))
      alpha[,s]=mvrnorm(1,amu,acov)
    }

    ### update rho
    if(T>1){
      tempv = matrix(0,S,S)
      tempm = rep(0,S)
      for(t in 2:T){
        gs0 = matrix(0,q*S,S)
        TIDg = c()
        for(s in 1:S){
          gs0[((s-1)*q+1):(s*q),s] = gg[((t-2)*q+1):((t-1)*q),s]
          TIDg= c(TIDg, ((s-1)*q*T+(t-1)*q+1):((s-1)*q*T+t*q))
        }
        tempv = tempv + t(gs0)%*%kronecker(solve(Sig),KtQK[,,t])%*%gs0
        tempm = tempm + t(gs0)%*%kronecker(solve(Sig),KtQK[,,t])%*%gvec[TIDg]
      }
      vv = solve(tempv)
      mm = vv%*%tempm
      ### simulate from truncated normal
      rhovec = as.vector(rtmvnorm(n=1, mean = as.vector(mm), sigma=vv,lower=rep(0,S), upper=rep(1,S)))
      rho = diag(rhovec)
    }

    ### update Sig
    At = matrix(0,S,S)
    t=1
    gt = as.vector(gg[((t-1)*q+1):(t*q),])
    C = gt%*%t(gt)
    for(i in 1:S){
      for(j in 1:S){
        At[i,j] = sum(diag(C[((i-1)*q+1):(i*q),((j-1)*q+1):(j*q)]%*%KtQK[,,t]))
      }
    }
    for(t in 2:T){
      gt = as.vector(gg[((t-1)*q+1):(t*q),])
      gt0 = as.vector(gg[((t-2)*q+1):((t-1)*q),])
      C = (gt-kronecker(rho,diag(q))%*%gt0)%*%t((gt-kronecker(rho,diag(q))%*%gt0))
      for(i in 1:S){
        for(j in 1:S){
          At[i,j] = At[i,j]+sum(diag(C[((i-1)*q+1):(i*q),((j-1)*q+1):(j*q)]%*%KtQK[,,t]))
        }
      }
    }
    nunew = Signu+q*T
    psinew = SigPsi+At
    Sig = riwish(nunew,psinew)

    if(m %% 1000 == 0){message(m)}

    if (m %% M.thin == 0){
      if (m<=M.burn) {
        betaburn[,m/M.thin,]=beta
        alphaburn[,m/M.thin,] = alpha
        sigburn[,m/M.thin] = matrix(Sig,nrow=S^2,byrow=T)
      } ## closes if
      else {
        Zout[,(m-M.burn)/M.thin,] = z
        alphaout[,(m-M.burn)/M.thin,] = alpha
        betaout[,(m-M.burn)/M.thin,]=beta
        gamout[,(m-M.burn)/M.thin,] = gg
        sigout[,(m-M.burn)/M.thin] = matrix(Sig,nrow=S^2,byrow=T)
        rhoout[,(m-M.burn)/M.thin] = diag(rho)
        for(s in 1:S){
          psiout[,(m-M.burn)/M.thin,s] = pnorm(X%*%alpha[,s]+K%*%gg[,s])
          pout[,(m-M.burn)/M.thin,s] = pnorm(W%*%beta[,s])

          ## This requires expanding psi
          psi.temp = pnorm(X%*%alpha[,s]+K%*%gg[,s])
          hold = data.frame(zsite, zseasonInd, psi.temp)
          colnames(hold) = c("siteID", "season", "psi.temp")
          expanded.psi = merge(model.input$detection.info, hold, by=c("siteID","season"), all.x=TRUE)
          expanded.psi = expanded.psi$psi.temp

          loglik[,(m-M.burn)/M.thin,s] = y[,s]*log(onI*expanded.psi*pout[,(m-M.burn)/M.thin,s])+(rep(1,nrow(W))-y[,s])*log(1-onI*expanded.psi*pout[,(m-M.burn)/M.thin,s])
          II = which(is.na(loglik[,(m-M.burn)/M.thin,s]))
          loglik[II,(m-M.burn)/M.thin,s] = 0
        }
      }
    }

    ## Save output every 1000 iterations.  Write over file.
    if (sv == "TRUE"){
      if (m %% 1000 == 0) {
        save(file="MCMCoutput.Rdata", list=c("alphaout","betaout","Zout","sigout","gamout","psiout","pout"))
      }
    }
  }

  run.time = Sys.time()-st

  WAIC = rep(0,S)
  for(s in 1:S){
    lppd = sum(log(rowMeans(exp(loglik[,,s]))))
    pWAIC2 = sum(apply(loglik[,,s],2,FUN=var ))
    WAIC[s] = -2*(lppd-pWAIC2)
  }

  #out=list(alphaout,betaout,Zout,alphaburn,betaburn,sigburn,sigout,gamout,gamburn,run.time)
  out = list("alpha"=alphaout,"beta"=betaout,"Sigma"=sigout,"rho"=rhoout, "gamma" = gamout, "psi"=psiout,"p"=pout,"run.time"=run.time,"WAIC"=WAIC,"occupancy.info"=model.input$occupancy.info,"detection.info"=model.input$detection.info[,c("siteID","site","season","survey")],"X"=X,"W"=W,"y"=y, "basis.K" = basis[[2]])

  ##.GlobalEnv$mcmc.output = out
  return(out)
}

