#' This function constructs basis functions.  It assumes coordinates form a metric.
#'
#' @param q The desired number of basis functions.  Must be an integer greater than or equal to 1.
#' @param model.input A list of output created by running the create.data.R function
#' @return A list with
#' \itemize{
#'   \item K spatial basis functions
#'   \item KtQK which is literally the matrix operation transpose(K) times Q times K, and is the variance of the multivariate #' random effect gamma.
#' }
#' @export

MakeBasis = function(q=q, model.input){

  coords = model.input$occupancy.info[,c("x","y")]
  X = model.input$X
  A = model.input$A
  seasonInd = model.input$occupancy.info$season

  T = length(unique(seasonInd))
  nx = 40
  ny = 40

  xo = seq(min(coords$x),max(coords$x),length.out=nx)
  yo = seq(min(coords$y),max(coords$y),length.out=ny)

  K = matrix(0,dim(A)[1],q*T)
  for(t in 1:T){
    Tind = which(seasonInd==t)
    P = diag(length(Tind))-X[Tind,]%*%MASS::ginv(t(X[Tind,])%*%X[Tind,])%*%t(X[Tind,])
    OM = P%*%A[Tind,Tind]%*%P
    S0 = Re(eigen(OM)$values)
    ## eig2keep = which(S0>0)[1:q]
    eig2keep = 1:q
    UU = Re(eigen(OM)$vectors)[,eig2keep]
    temp = UU
    #### note eigenvectors only unique up to multiplicative constant of +1 or -1
    #### so we need to make sure the vectors are chosen to vary smoothly over time and do not oscillate
    if(t>1){
      Tind0= which(seasonInd==(t-1))
      for(j in 1:q){
        ### check to see if it oscillates from the previous time period or not
        ### need to expand to all locations
        Ugridold = interp::interp(coords$x[Tind0],coords$y[Tind0],tempold[,j],xo=xo,yo=yo,linear=TRUE,extrap=TRUE)
        Ugrid = interp::interp(coords$x[Tind],coords$y[Tind],temp[,j],xo=xo,yo=yo,linear=TRUE,extrap=TRUE)
        Ugridold$z[which(is.na(Ugridold$z))]=0
        Ugrid$z[which(is.na(Ugrid$z))]=0

        if(norm(Ugrid$z-Ugridold$z)>norm(-1*Ugrid$z-Ugridold$z)){
          temp[,j]=-1*temp[,j]
        }
      }
    }
    tempold = temp
    K[Tind,((t-1)*q+1):(t*q)]=tempold
  }

  KtQK = array(0,dim=c(q,q,T))
  for(t in 1:T){
    Tind = which(seasonInd==t)
    At = A[Tind,Tind]
    Dt = diag(rowSums(At))
    KtQK[,,t] = t(K[Tind,((t-1)*q+1):(t*q)])%*%(Dt-At)%*%K[Tind,((t-1)*q+1):(t*q)]
  }

  out=list(KtQK, K)
  return(out)
}
