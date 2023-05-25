#' This function runs the MCMC.
#'
#' @param M.iter The total number of iterations in MCMC
#' @param M.burn The length of the burn in
#' @param M.thin The number to thin the chain.  Thinning by 10 only stores every 10th run.
#' @param model.input A list of output created by running the create.data.R function
#' @param q Desired number of Moran's I basis functions in the restricted spatial regression model
#' @param sv A TRUE/FALSE on whether or not the MCMC output should be saved as 'MCMC.Rdata' and overwritten every 1000 iterations.  Defaults to false.
#' @param every A number to determine how frequently MCMC output is saved along the chain.  Defaults to 1000.
#' @param WAIC A TRUE/FALSE on whether or not the MCMC should compute and save WAIC.  Defaults to false.
#' @param param2keep A character vector that governs which outputs are saved.  Permissible entries are "alpha", "beta", "gamma", "rho", "sigma", "psi", "z", "p", and "loglik"

#' @return A list with all standard MCMC output
#' @import MASS tmvtnorm truncnorm
#' @export
#' @examples
#' head(detection)
#' head(occupancy)
#' head(coords)
#' DataNames = list("species"=colnames(detection)[4:9],
#' "detection"=c("duration"),"occupancy"=c("forest","elev"))
#' model.input = multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000)
#' out = GibbsSampler(M.iter=3, M.burn=1, M.thin=1, model.input, q=10, sv=FALSE)



GibbsSampler = function(M.iter, M.burn=NULL, M.thin=NULL, model.input, q=NULL, sv=FALSE, every=1000,WAIC=FALSE,param2keep=c("alpha","beta","gamma","rho","sigma","psi")){

  ### check that param2keep only contains variables that can be saved
  if( min(as.numeric(param2keep %in% c("alpha","beta","gamma","rho","sigma","psi","z","p","loglik")))==0){
    stop("param2keep must be a character vector, whose only permissible entries are `alpha`, `beta`, `gamma`, `rho`, `sigma`, `psi`, `z`, `p`, and `loglik`")
  }


  if(q>min(summary(as.factor(model.input$occupancy.info$season)))){
    stop("Number of basis functions cannot be larger than the smallest number of sites in any season")
  }

  if(!is.null(q)){
    basis = MakeBasis(q,model.input)
  }else{
    q = ceiling(.1*min(summary(as.factor(model.input$occupancy.info$season))))
    basis = MakeBasis(q,model.input)
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

  KtQK=basis[[1]]
  q=dim(KtQK)[1]
  K=basis[[2]]
  zseasonInd = model.input$occupancy.info$season
  zsite = model.input$occupancy.info$siteID

  T = max(zseasonInd) ### number of seasons

  if(T==1 & 'rho' %in% param2keep){
    param2keep <- param2keep[param2keep!='rho']
  }

  nT = nrow(X)

  S = dim(y)[2] #number of species
  y[onI==0,]=rep(0,S)
  ysum = matrix(0,nT,S)
  for(i in 1:nT){
    II = which(site==zsite[i] & seasonInd == zseasonInd[i])
    ysum[i,] = colSums(matrix(y[II,],length(II),S))
  }
  z = 1*(ysum>0) ### if y>0 for any site visit then z=1

  ap = dim(X)[2]
  bp = dim(W)[2]
  acov = solve(t(X)%*%X)

  ### set initial values for MCMC
  alpha = matrix(0,ap,S)
  beta = matrix(0,bp,S)
  Sig = diag(S)
  if(T>1){
    rho = 0.9*diag(S)
  }else{
    rho=0 ### no time dependence in the single season setting...just make 0
  }
  M.rho = kronecker(rho,diag(q)) ##temporal propagator matrix

  ### KtQKstar and Sigi get updated when Sig gets updated
  Sigi = solve(Sig)
  KtQKstar = array(0,c(S*q,S*q,T))
  for(t in 1:T){
    KtQKstar[,,t] = kronecker(Sigi,KtQK[,,t])
  }

  gg = matrix(0,q*T,S) ### gg = gamma in paper
  gvec = as.vector(gg)
  ztilde = z-0.5 ### need initial value negative if z<0 and positive if z>0.
  ytilde = y-0.5
  psi = pnorm(X%*%alpha + K%*%gg)
  p = pnorm(W%*%beta)

  loglik = matrix(0,nrow(W),S)
  for(s in 1:S){
    p.temp = p[,s]
    ## This requires expanding psi
    psi.temp = psi[,s]
    hold = data.frame(zsite, zseasonInd, psi.temp)
    colnames(hold) = c("siteID", "season", "psi.temp")
    expanded.psi = merge(model.input$detection.info, hold, by=c("siteID","season"), all.x=TRUE)
    expanded.psi = expanded.psi$psi.temp
    loglik[,s] = y[,s]*log(onI*expanded.psi*p.temp)+(rep(1,nrow(W))-y[,s])*log(1-onI*expanded.psi*p.temp)
    II = which(is.na(loglik[,s]))
    loglik[,s] = 0
  }

  ### hyperparameters
  Signu = S+2
  SigPsi = diag(S)

  ### create arrays to store output for all variables in param2keep

  output = list()
  for(v in 1:length(param2keep)){
    if(param2keep[v]=="sigma"){
      output[[v]] = matrix(0,S^2,(M.iter-M.burn)/M.thin)
    }else if(param2keep[v]=="gamma"){
      output[[v]] = matrix(0,q*T*S,(M.iter-M.burn)/M.thin)
    }else if(param2keep[v]=="rho"){
      output[[v]] = matrix(0,S,(M.iter-M.burn)/M.thin)
    }else
    {
      output[[v]] = matrix(0,S*dim(eval(as.name(param2keep[v])))[1],(M.iter-M.burn)/M.thin)
    }
  }
  names(output) <- param2keep

  if(WAIC==TRUE){
    loglik0 = matrix(0,S*nrow(W),(M.iter-M.burn)/M.thin)
  }

  #### Note we only need to update the z's where the corresponding y=0
  #### If y=1, then z has to be 1. This is already true from the initial
  #### condition, so those z's never need to be updated...

  nyz = dim(ysum)[1]-colSums(1*(ysum)>0) ### how many of the y's are not 1
  zInd = list()
  Indy0 = list()
  Indy1 = list()
  for(s in 1:S){
    zInd[[s]] = which(ysum[,s]==0)
    Indy0[[s]] = which(y[,s]==0)
    Indy1[[s]] = which(y[,s]==1)
  }


  progress_bar = txtProgressBar(min=0, max=M.iter, style = 3, char="=")
  st = Sys.time()

  for(m in 1:M.iter){

    ### note psi only changes when alpha, gamma updated
    ### p only changes when beta updated. Compute these once and then update as needed
    psi0 = X%*%alpha + K%*%gg
    psi = pnorm(psi0)
    p0 = W%*%beta
    p = pnorm(p0)

    ####### update the z's where the corresponding y are 0 (or missing) for each site visit
    ## for each species
    for (s in 1:S){
      ### need the product of the 1-p0 for each site/season combo
      temp  = p[,s]*onI
      temp[onI==0]=0
      pcit = aggregate(1-temp,by=list(site,seasonInd),FUN=prod)[,3]
      psibar0 = (psi[,s]*pcit)/(1-psi[,s]+psi[,s]*pcit)
      z[zInd[[s]],s] = (log(runif(nyz[s]))<log(psibar0[zInd[[s]]]))
    }

    ######## update each ztilde
    for (s in 1:S){
      Ind0 = which(z[,s]==0)
      Ind1 = which(z[,s]==1)
      ztilde[Ind0,s] = truncnorm::rtruncnorm(n=length(Ind0), a=-Inf, b=0, mean = psi0[Ind0,s], sd = 1)
      ztilde[Ind1,s] = truncnorm::rtruncnorm(n=length(Ind1), a=0, b=Inf, mean = psi0[Ind1,s], sd = 1)
    }

    ######## update each ytilde
    for (s in 1:S){
      ytilde[Indy0[[s]],s] = truncnorm::rtruncnorm(n=length(Indy0[[s]]), a=-Inf, b=0, mean = p0[Indy0[[s]],s], sd = 1)
      ytilde[Indy1[[s]],s] = truncnorm::rtruncnorm(n=length(Indy1[[s]]), a=0, b=Inf, mean = p0[Indy1[[s]],s], sd = 1)
    }

    ######## update gvec and gg
    if(T>1){
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
      gp = diag(S*q)+KtQKstar[,,t]+t(M.rho)%*%KtQKstar[,,t+1]%*%M.rho
      gcov=chol2inv(chol(gp))
      gmu = gcov%*%(t(kronecker(diag(S),K[Tind,((t-1)*q+1):(t*q)]))%*%as.vector(ztilde[Tind,]-X[Tind,]%*%alpha)+t(M.rho)%*%KtQKstar[,,t+1]%*%gvec[TIDg1])
      gstar = MASS::mvrnorm(1,gmu,gcov)
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
          gp = diag(S*q)+KtQKstar[,,t]+t(M.rho)%*%KtQKstar[,,t+1]%*%M.rho
          gcov=chol2inv(chol(gp))
          gmu = gcov%*%(t(kronecker(diag(S),K[Tind,((t-1)*q+1):(t*q)]))%*%as.vector(ztilde[Tind,]-X[Tind,]%*%alpha)+t(M.rho)%*%KtQKstar[,,t+1]%*%gvec[TIDg1]+KtQKstar[,,t]%*%M.rho%*%gvec[TIDg0])
          gstar = MASS::mvrnorm(1,gmu,gcov)
          for (s in 1:S){
            gvec[TIDg[((s-1)*q+1):(s*q)]] =  gstar[((s-1)*q+1):(s*q)]-mean(gstar[((s-1)*q+1):(s*q)])
          }
        }
      }
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
      gp = diag(S*q)+KtQKstar[,,t]
      gcov=chol2inv(chol(gp))
      gmu = gcov%*%(t(kronecker(diag(S),K[Tind,((t-1)*q+1):(t*q)]))%*%as.vector(ztilde[Tind,]-X[Tind,]%*%alpha)+KtQKstar[,,t]%*%M.rho%*%gvec[TIDg0])
      gstar = MASS::mvrnorm(1,gmu,gcov)
      for (s in 1:S){
        gvec[TIDg[((s-1)*q+1):(s*q)]] =  gstar[((s-1)*q+1):(s*q)]-mean(gstar[((s-1)*q+1):(s*q)])
      }
    }else{
      ### if T=1 case
      t=1
      TID = c() ### observations from same time period
      TIDg = c() ### rows of g from same time period
      TIDg1 = c() ### rows of g from next time period
      Tind = which(zseasonInd==t)
      for(s in 1:S){
        TID = c(TID, (s-1)*nT+Tind)
        TIDg = c(TIDg, ((s-1)*q*T+(t-1)*q+1):((s-1)*q*T+t*q))
      }
      ### compute precision matrix for full conditional of gamma_t
      gp = diag(S*q)+KtQKstar[,,t]
      gcov=chol2inv(chol(gp))
      gmu = gcov%*%(t(kronecker(diag(S),K[Tind,((t-1)*q+1):(t*q)]))%*%as.vector(ztilde[Tind,]-X[Tind,]%*%alpha))
      gstar = MASS::mvrnorm(1,gmu,gcov)
      for (s in 1:S){
        gvec[TIDg[((s-1)*q+1):(s*q)]] =  gstar[((s-1)*q+1):(s*q)]-mean(gstar[((s-1)*q+1):(s*q)])
      }
    }
    gg = matrix(gvec,q*T,S)
    psi0 = X%*%alpha + K%*%gg
    psi = pnorm(psi0)

    #### update beta
    ## This requires expanding z
    z.temp = z
    colnames(z.temp) = paste("z",1:dim(z.temp)[2],sep=".")
    hold = data.frame(zsite, zseasonInd, z.temp)
    colnames(hold) = c("siteID", "season",colnames(z.temp))
    expanded.z = merge(model.input$detection.info, hold, by=c("siteID","season"), all.x=TRUE)
    expanded.z = as.matrix(expanded.z[,colnames(z.temp)])

    for (s in 1:S){
      #### pull out the rows of y for which the corresponding z and onI is 1
      zstar = which(expanded.z[,s]*onI==1)
      bcov = solve(t(W[zstar,])%*%W[zstar,])
      bmu = bcov%*%t(W[zstar,])%*%ytilde[zstar,s]
      beta[,s] = MASS::mvrnorm(1,bmu,bcov)
    }
    p0 = W%*%beta
    p = pnorm(p0)

    ### update alpha
    for (s in 1:S){
      amu = acov%*%(t(X)%*%(ztilde[,s]-K%*%gg[,s]))
      alpha[,s]=MASS::mvrnorm(1,amu,acov)
    }
    psi0 = X%*%alpha + K%*%gg
    psi = pnorm(psi0)

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
        tempv = tempv + t(gs0)%*%KtQKstar[,,t]%*%gs0
        tempm = tempm + t(gs0)%*%KtQKstar[,,t]%*%gvec[TIDg]
      }
      vv = solve(tempv)
      mm = vv%*%tempm
      ### simulate from truncated normal
      rhovec = as.vector(tmvtnorm::rtmvnorm(n=1, mean = as.vector(mm), sigma=vv,lower=rep(0,S), upper=rep(1,S)))

      if (length(rhovec) == 1){rho = rhovec} else {
        rho = diag(rhovec)
      }
      M.rho = kronecker(rho,diag(q))
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
    if(T>1){
      for(t in 2:T){
        gt = as.vector(gg[((t-1)*q+1):(t*q),])
        gt0 = as.vector(gg[((t-2)*q+1):((t-1)*q),])
        C = (gt-M.rho%*%gt0)%*%t((gt-M.rho%*%gt0))
        for(i in 1:S){
          for(j in 1:S){
            At[i,j] = At[i,j]+sum(diag(C[((i-1)*q+1):(i*q),((j-1)*q+1):(j*q)]%*%KtQK[,,t]))
          }
        }
      }
    }
    nunew = Signu+q*T
    psinew = SigPsi+At
    Sig = MCMCpack::riwish(nunew,psinew)

    Sigi = solve(Sig)
    for(t in 1:T){
      KtQKstar[,,t] = kronecker(Sigi,KtQK[,,t])
    }

    if (m>M.burn & m %% M.thin == 0){
      if("loglik" %in% param2keep || WAIC==TRUE){
        for(s in 1:S){
          p.temp = pnorm(W%*%beta[,s])
          ## This requires expanding psi
          psi.temp = pnorm(X%*%alpha[,s]+K%*%gg[,s])
          hold = data.frame(zsite, zseasonInd, psi.temp)
          colnames(hold) = c("siteID", "season", "psi.temp")
          expanded.psi = merge(model.input$detection.info, hold, by=c("siteID","season"), all.x=TRUE)
          expanded.psi = expanded.psi$psi.temp
          loglik0[((s-1)*nrow(W)+1):(s*nrow(W)),(m-M.burn)/M.thin] = y[,s]*log(onI*expanded.psi*p.temp)+(rep(1,nrow(W))-y[,s])*log(1-onI*expanded.psi*p.temp)
          II = which(is.na(loglik0[((s-1)*nrow(W)+1):(s*nrow(W)),(m-M.burn)/M.thin]))
          loglik0[((s-1)*nrow(W)+II),(m-M.burn)/M.thin] = 0
        }
        loglik = loglik0[,(m-M.burn)/M.thin]
      }
      for(v in 1:length(param2keep)){
        if(param2keep[v]=="sigma"){
          output[[v]][,(m-M.burn)/M.thin] = matrix(Sig,nrow=S^2,byrow=T)
        }else if(param2keep[v]=="gamma"){
          output[[v]][,(m-M.burn)/M.thin] = gvec
        }else if(param2keep[v]=="rho"){
          output[[v]][,(m-M.burn)/M.thin] = rhovec
        }else{
          output[[v]][,(m-M.burn)/M.thin] = as.vector(get(param2keep[v]))
        }
      }
    }


    ## Save output every 1000 iterations.  Write over file.
    if (m %% every == 0 & m>M.burn & sv=="TRUE") {
      save(file="MCMCoutput.Rdata", output)
    }

    setTxtProgressBar(progress_bar, value = m)

  }

  close(progress_bar)
  run.time = Sys.time()-st

  ### make list of MCMC objects for the output
  alpha.names <- rep("NA",S*ap)
  beta.names <- rep("NA",S*bp)
  gamma.names <- rep("NA",S*q)
  rho.names <- rep("NA",S)
  sigma.names <- rep("NA",S^2)
  psi.names <- rep("NA",nrow(X)*S)
  z.names <- rep("NA",nrow(X)*S)
  p.names <- rep("NA",nrow(X)*S)
  loglik.names <- rep("NA",nrow(W)*S)
  for(s in 1:S){
    alpha.names[((s-1)*ap+1):(s*ap)] <- paste(model.input$DataNames$species[s],c("Int",model.input$DataNames$occupancy), sep=" ")
    beta.names[((s-1)*bp+1):(s*bp)] <- paste(model.input$DataNames$species[s],c("Int",model.input$DataNames$detection), sep=" ")
    gamma.names[((s-1)*q*T+1):(s*q*T)] <- paste(model.input$DataNames$species[s],paste("gamma",kronecker(rep(1,T),1:q),"season",kronecker(1:T,rep(1,q)),sep=" "),sep=" " )
    rho.names[s] <- paste(model.input$DataNames$species[s],"rho",sep=" ")
    sigma.names[((s-1)*S+1):(s*S)] <- paste(model.input$DataNames$species,"x", model.input$DataNames$species[s],"Covariance", sep=" ")
    psi.names[((s-1)*nrow(X)+1):(s*nrow(X))] <- paste(model.input$DataNames$species[s],paste("site",model.input$occupancy.info$site, "season", model.input$occupancy.info$season,sep=" "),"psi", sep=" ")
    z.names[((s-1)*nrow(X)+1):(s*nrow(X))] <- paste(model.input$DataNames$species[s],paste("site",model.input$occupancy.info$site, "season", model.input$occupancy.info$season,sep=" "),"z", sep=" ")
    p.names[((s-1)*nrow(W)+1):(s*nrow(W))] <- paste(model.input$DataNames$species[s],paste("site",model.input$detection.info$site, "season", model.input$detection.info$season,"survey",model.input$detection.info$survey,sep=" "),"p", sep=" ")
    loglik.names[((s-1)*nrow(W)+1):(s*nrow(W))] <- paste(model.input$DataNames$species[s],paste("site",model.input$detection.info$site, "season", model.input$detection.info$season,"survey",model.input$detection.info$survey,sep=" "),"loglik", sep=" ")
  }

  for(v in 1:length(param2keep)){
    output[[v]] = coda::as.mcmc(t(output[[v]]))
    coda::varnames(output[[v]]) <- get(paste(param2keep[v],".names",sep=""))
  }


  if(WAIC==TRUE){
    WAIC.out = rep(0,S)
    for(s in 1:S){
      lppd = sum(log(rowMeans(exp(loglik0[((s-1)*nrow(W)+1):(s*nrow(W)),]))))
      pWAIC2 = sum(apply(loglik0[((s-1)*nrow(W)+1):(s*nrow(W)),],2,FUN=var ))
      WAIC.out[s] = -2*(lppd-pWAIC2)
    }
    out = list("samples"=output, "run.time"=run.time,"WAIC"=WAIC.out,"occupancy.info"=model.input$occupancy.info,"detection.info"=model.input$detection.info[,c("siteID","site","season","survey")],"X"=X,"W"=W,"y"=y, "basis.K" = basis[[2]])
  }else{
    out = list("samples"=output, "run.time"=run.time,"occupancy.info"=model.input$occupancy.info,"detection.info"=model.input$detection.info[,c("siteID","site","season","survey")],"X"=X,"W"=W,"y"=y, "basis.K" = basis[[2]])
  }

  return(out)
}
