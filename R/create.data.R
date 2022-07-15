#####################################################################################################
##### write a function that will take two data frames - one for occupancy and one for detection #####
##### and output all matrices and vectors needed for the MCMC algorithm in the correct order ########
##### take a matrix of the coordinates for all sites in the study and use it to output adjacency ####
##### input a list with three elements - names of species, detection covariates, occupancy covs #####
#####################################################################################################

#' This function creates model.input for the run.mcmc() function
#'
#' @param detection a data frame that is one row for every site X season X survey combination
#' contains columns for the site, season, survey within season, all covariates to be used
#' in the detection process of the model, and binary indicators of detections for all species
#' to be modeled
#' @param occupancy A data frame that is one row for every site x season combination
#' contains columns for the site, season, and all covariates to be used in the latent occupancy
#' process of the model
#' @param coords is a data frame that is one row for every site included in the study
#' contains columns for the site, x, and y
#' output adjacency based on Euclidean distance threshold the user provides as a function input
#' @param names is a list with elements "species", "detection", and "occupancy"
#' names$species is a vector with the name of every species in the study
#' names$detection is a vector with the names of the detection covariates
#' names$occupancy is a vector with the names of the occupancy covariates
#' @param threshold is the distance which determines if two locations are neighbors
#' in the adjacency matrix or not.
#' @return model.input a list with
#' \itemize{
#' \item X design matrix for occupancy
#' \item W design matriw for detection
#' \item y observed data
#' \item A adjacency matrix
#' \item detection.info details for detection
#' \item occupancy.info details for occupancy
#' }
#' @export

create.data <- function(detection,occupancy,coords,names,threshold){

  detection$site = as.factor(detection$site)
  occupancy$site = as.factor(occupancy$site)
  coords$site = as.factor(coords$site)
  detection$siteID = as.numeric(detection$site)
  detloc = detection[,c("site","siteID")]
  detloc = detloc[!duplicated(detloc),]
  occupancy = merge(occupancy,detloc,by="site",all.x=TRUE,all.y=FALSE)
  coords = merge(coords,detloc,by="site",all.x=TRUE,all.y=FALSE)

  ### perform some checks
  ### should only be one row for each site x season combination in occupancy and one row for each site x season x survey combination in detection

  detection.obs = detection[,c("siteID","site","season","survey")]
  occ.obs = occupancy[,c("siteID","site","season")]
  if(max(as.numeric(duplicated(detection.obs)))==1){
    stop("Duplicated entries in detection")
  }
  if(max(as.numeric(duplicated(occ.obs)))==1){
    stop("Duplicated entries in occupancy")
  }

  #### note that the model cannot handle NA in W or X

  #### identify the W where y is not missing but W is.
  #### Remove those site x season combos from all 3 data frames
  dp = length(names$detection)
  na.site.season.survey = data.frame("site"=rep(NA,0),"season"=rep(NA,0),"survey"=rep(NA,0))
  na.site.season = data.frame("site"=rep(NA,0),"season"=rep(NA,0))
  #a0 = 1
  for(i in 1:dp){
    na.W = which(is.na(detection[,names$detection[i]]))
    II = which(is.na(rowSums(detection[na.W,names$species]))==FALSE) ## rows with observed responses
    na.obs=na.W[II]
    if(length(na.obs)>0){
      na.site.season.survey=rbind(na.site.season.survey,data.frame("site"=detection$site[na.obs],"season"=detection$season[na.obs],"survey"=detection$survey[na.obs]))
  #    a0 = a0+length(na.obs)
    }
  }

  ## Check for fully missing cases, meaning all surveys are missing for a site/season combination.
  ## Need to remove those site/season rows from X if this is the case.
  hold = unique(na.site.season.survey[,c("site","season")])
  if(dim(hold)[1]>0){
    for (i in 1:dim(hold)[1]){
      num.sur = dim(na.site.season.survey[na.site.season.survey$site == hold[i,1] & na.site.season.survey$season == hold[i,2],])[1]
      num.tot = dim(detection[detection$site == hold[i,1] & detection$season == hold[i,2],])[1]
      if (num.sur == num.tot){na.site.season = rbind(na.site.season,hold[i,])}
    }
  }

  ### remove these site x season combos from detection and occupancy
  if(dim(na.site.season.survey)[1]>0){
    toremovedet = list() ## this will be all site/season/survey combos with NAs as covariates
    for(i in 1:dim(na.site.season.survey)[1]){
      toremovedet[[i]] = which(detection$site==na.site.season.survey$site[i] & detection$season==na.site.season.survey$season[i] & detection$survey==na.site.season.survey$survey[i])
    }
    detection = detection[-unlist(toremovedet),]
    message("Warning: Rows in detection with missing covariates have been removed.")
  }

  if(dim(na.site.season)[1]>0){
    toremoveocc = rep(NA,dim(na.site.season)[1]) ## this will only be site/season combos for which *all* surveys have NAs as covariates
    for(i in 1:dim(na.site.season)[1]){
      toremoveocc[i] = which(occupancy$site==na.site.season$site[i] & occupancy$season==na.site.season$season[i])
    }
    occupancy = occupancy[-toremoveocc,]
    message("Warning: Rows in occupancy have been removed corresponding to missing rows in detection.")
  }

  ### repeat for X matrix
  op = length(names$occupancy)
  na.site.season=data.frame("site"=rep(NA,0),"season"=rep(NA,0))
  #a0 = 1
  for(i in 1:op){
    na.obs = which(is.na(occupancy[,names$occupancy[i]]))
    if(length(na.obs)>0){
      na.site.season=rbind(na.site.season,data.frame("site"=occupancy$site[na.obs],"season"=occupancy$season[na.obs]))
  #    a0 = a0+length(na.obs)
    }
  }
  ### remove these site x season combos from detection and occupancy
  if(dim(na.site.season)[1]>0){
    toremovedet = list()
    toremoveocc = rep(NA,dim(na.site.season)[1])
    for(i in 1:dim(na.site.season)[1]){
      toremovedet[[i]] = which(detection$site==na.site.season$site[i] & detection$season==na.site.season$season[i])
      toremoveocc[i] = which(occupancy$site==na.site.season$site[i] & occupancy$season==na.site.season$season[i])
    }
    detection = detection[-unlist(toremovedet),]
    occupancy = occupancy[-toremoveocc,]
    message("Warning: Rows in occupancy with missing covariates have been removed.  Corresponding rows in detection have also been removed.")
  }

  ### should have a row in occupancy for every site x season combination in detection - if not, then give an error
  detection.obs = detection[,c("siteID","site","season")]
  occ.obs = occupancy[,c("siteID","site","season")]
  detection.site.season = detection.obs[!duplicated(detection.obs[,c("site","season")]),c("site","season")]
  if(nrow(detection.site.season) != nrow(occ.obs)){
    stop("Different site/season combinations in detection and occupancy")
  }

  ### merge coords to occupancy data
  occupancy = merge(occupancy,coords[,c("site","x","y")],by="site")

  ####################################################
  ## Was a site/season/survey combo visited?
  detection$observations = 1-is.na(rowMeans(detection[,names$species]))

  ### order detection, occupancy, coords as desired
  detection = detection[order(detection$season,detection$site,detection$survey),]
  occupancy = occupancy[order(occupancy$season,occupancy$site),]
  y = detection[,names$species]

  if (max(y[is.na(y)==FALSE])>1) {stop("Forbidden values in the detection matrix.  Only 1, 0, and NA are permissible entries for detection.  Counts should be converted to binary detections.")}


  ### create design matrices.
  X = cbind(rep(1,nrow(occupancy)),as.matrix(occupancy[,names$occupancy]))
  W = cbind(rep(1,nrow(detection)),as.matrix(detection[,names$detection]))

  ### make site info matrix that has same number of rows as X
  dist.mat = as.matrix(dist(as.matrix(occupancy[,c("x","y")])))

  A = 1*(dist.mat<=threshold & dist.mat>0)
  ### need 0 if in different seasons
  time.mat = matrix(0,nrow(A),nrow(A))
  for(t in 1:max(occupancy$season)){
    II = which(occupancy$season==t)
    time.mat[II,II]=1
  }
  A = A*time.mat

  ### output a list that contains X,W,y,A,season,site,survey,coords
  model.input = list("X"=X,"W"=W,"y"=y,"A"=A,"detection.info"=detection[,c("siteID","site","season","survey","observations")],"occupancy.info"=occupancy[,c("siteID","site","season","x","y")])

  #.GlobalEnv$model.input = model.input
  return(model.input)
}
