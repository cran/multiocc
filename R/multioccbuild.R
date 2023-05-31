#####################################################################################################
##### write a function that will take two data frames - one for occupancy and one for detection #####
##### and output all matrices and vectors needed for the MCMC algorithm in the correct order ########
##### take a matrix of the coordinates for all sites in the study and use it to output adjacency ####
##### input a list with three elements - names of species, detection covariates, occupancy covs #####
#####################################################################################################

#' This function creates model.input for the GibbsSampler() function
#'
#' @param detection A data frame that has one row for every site X season X survey combination.
#' Must contain columns exactly named 'site', 'season', and 'survey' within season.
#' Must also contain all covariates in the detection process of the model, and binary indicators of detections for all species
#' to be modeled.  It is permissible for this data frame to have columns for species and/or variables
#' that will not be used in model.
#'
#' @param occupancy A data frame that is one row for every site x season combination.
#' Must contain columns for the 'site' and 'season', and these must be named 'site' and 'season' exactly.
#' Also must contain all covariates to be used in the latent occupancy process of the model.
#' It is permissible for this data frame to have columns for species and/or variables
#' that will not be used in model.
#'
#' @param coords A data frame that is one row for every site included in the study.
#' Contains columns for the 'site', and location coordinates x and y.
#' These are used to output the adjacency matrix A based on Euclidean distance threshold the user provides
#' as an input in the 'DataNames' argument.
#'
#' @param DataNames A list with elements "species", "detection", and "occupancy"
#' DataNames$species is a vector with the name of every species to be included in the model.  Must be a subset of names of columns of 'detection'.
#' DataNames$detection is a vector with the names of the detection covariates to be included in the model.  These names must be a subset of column names of 'detection'.
#' DataNames$occupancy is a vector with the names of the occupancy covariates to be included in the model.  These names must be a subset of column names of 'occupancy'.
#'
#' This list 'DataNames' is required because it:
#' (1) allows for modeling subsets of species and/or variables in varioys input data frames, which means the user does not need to modify either data frame for different runs of the model.
#' (2) this list also determines the order of covariates in X and W.
#'
#' @param threshold The distance which determines if two locations are neighbors
#' in the adjacency matrix or not.  This threshold is the Euclidean distance based on the x and y coordinates input in 'coords'.
#'
#' @return model.input a list with
#' \itemize{
#' \item DataNames, the list with elements "species", "detection", and "occupancy".
#' \item X, the design matrix for occupancy.  Contains a column of 1s for the intercept, and one column for each variable in names$occupancy.
#' \item W, design matrix for detection.  Contains a column of 1s for the intercept, and one column for each variable in names$detection.
#' \item y, observed occupancy data.  Contains one column for each species listed in names$species.
#' \item A, the adjacency matrix containing either a 0 or 1 indicating if two locations are neighbors based on their distances and the 'threshold' argument.
#' \item detection.info, details for detection.
#' \item occupancy.info, details for occupancy.
#' }
#' @export

multioccbuild <- function(detection,occupancy,coords,DataNames,threshold){

  ### ensure detection and occupancy has columns named
  if (!"site" %in% colnames(detection)){
    stop("detection does not have a column named 'site'")
  }
  if (!"season" %in% colnames(detection)){
    stop("detection does not have a column named 'season'")
  }
  if (!"survey" %in% colnames(detection)){
    stop("detection does not have a column named 'survey'")
  }
  if (!"site" %in% colnames(occupancy)){
    stop("occupancy does not have a column named 'site'")
  }
  if (!"season" %in% colnames(occupancy)){
    stop("occupancy does not have a column named 'season'")
  }

  detection$site = as.factor(detection$site)
  occupancy$site = as.factor(occupancy$site)
  coords$site = as.factor(coords$site)
  occupancy$siteID = as.numeric(occupancy$site)
  occloc = occupancy[,c("site","siteID")]
  occloc = occloc[!duplicated(occloc),]
  occloc2 = occupancy[,c("site","siteID","season")]
  occloc2 = occloc2[!duplicated(occloc2),]
  detection = merge(detection,occloc2,by=c("site","season"),all.y=TRUE,all.x=TRUE)

  if(length(which(is.na(detection$survey)))>0){
    detection$survey[which(is.na(detection$survey))]<-1
    message("Warning: Rows with missing values in detection have been added to correspond to additional site/season combinations present in occupancy.")
  }

  coords = merge(coords,occloc,by="site",all.y=TRUE,all.x=FALSE)

  ### perform some checks
  ### should only be one row for each site x season combination
  ### in occupancy and one row for each site x season x survey combination in detection

  detection.obs = detection[,c("siteID","site","season","survey")]
  occ.obs = occupancy[,c("siteID","site","season")]
  if(max(as.numeric(duplicated(detection.obs)))==1){
    stop("Duplicated entries in detection")
  }
  if(max(as.numeric(duplicated(occ.obs)))==1){
    stop("Duplicated entries in occupancy")
  }

  #### identify the W where y is not missing but W is.
  dp = length(DataNames$detection)
  for(i in 1:dp){
    na.W = which(is.na(detection[,DataNames$detection[i]]))
    if(length(na.W)>0){
      II = which(is.na(rowSums(as.matrix(detection[na.W,DataNames$species])))==FALSE) ## rows with observed responses
      if(length(II)>0){
        na.obs=na.W[II]
        #### for these rows, change the y to missing
        detection[na.obs,DataNames$species]=NA
        message("Warning: Rows in detection with missing covariates have been removed for purposes of fitting the model, but the site/season combination is retained in occupancy and therefore predictions will be outputted.")
      }
    }
  }

  ### if site/season combos in detection but not occupancy, remove and give warning
  detection.remov <- c()
  for(i in 1:dim(detection)[1]){
    II=which(occupancy$siteID==detection$siteID[i] & occupancy$season==detection$season[i])
    if(length(II)<1){
      detection.remov = c(detection.remov, i)
    }
  }
  if(length(detection.remov)>0){
    message("Warning: Rows in detection have been removed corresponding to missing rows in occupancy.")
    detection <- detection[-detection.remov,]
  }


  ### if site/season combos in occupancy but not detection, add rows with missing to detection
  #change<-0
  #for(i in 1:dim(occupancy)[1]){
  #  II=which(detection$siteID==occupancy$siteID[i] & detection$season==occupancy$season[i])
  #  if(length(II)<1){
  #    detection = rbind(detection,rep(NA,dim(detection)[2]))
  #    detection$site[dim(detection)[1]]=occupancy$site[i]
  #    detection$siteID[dim(detection)[1]]=occupancy$siteID[i]
  #    detection$season[dim(detection)[1]]=occupancy$season[i]
  #    detection$survey[dim(detection)[1]]=1
  #    change <- change+1
  #  }
  #}
  #if(change>0){
  #  message("Warning: Rows with missing values in detection have been added to correspond to additional site/season combinations present in occupancy.")
  #}

  #### remove site/season combinations that have missing occupancy covariate data
  op = length(DataNames$occupancy)
  na.site.season=data.frame("site"=rep(NA,0),"season"=rep(NA,0))
  for(i in 1:op){
    na.obs = which(is.na(occupancy[,DataNames$occupancy[i]]))
    if(length(na.obs)>0){
      na.site.season=rbind(na.site.season,data.frame("site"=occupancy$site[na.obs],"season"=occupancy$season[na.obs]))
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
    occupancy = occupancy[-toremoveocc,]
    if(length(unlist(toremovedet))>0){
      detection = detection[-unlist(toremovedet),]
    }
    message("Warning: Rows in occupancy with missing covariates have been removed.  Corresponding rows in detection have also been removed.")
  }

  ### merge coords to occupancy data
  occupancy = merge(occupancy,coords[,c("site","x","y")],by="site")

  ## Was a site/season/survey combo visited?
  detection$observations = 1-is.na(rowMeans(as.matrix(detection[,DataNames$species])))

  ### order detection, occupancy, coords as desired
  detection = detection[order(detection$season,detection$site,detection$survey),]
  occupancy = occupancy[order(occupancy$season,occupancy$site),]
  y = detection[,DataNames$species]

  if (max(y[is.na(y)==FALSE])>1) {stop("Forbidden values in the detection matrix.  Only 1, 0, and NA are permissible entries for detection.  Counts should be converted to binary detections.")}

  ### create design matrices.
  X = cbind(rep(1,nrow(occupancy)),as.matrix(occupancy[,DataNames$occupancy]))
  W = cbind(rep(1,nrow(detection)),as.matrix(detection[,DataNames$detection]))

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
  model.input = list("DataNames"=DataNames,"X"=X,"W"=W,"y"=y,"A"=A,"detection.info"=detection[,c("siteID","site","season","survey","observations")],"occupancy.info"=occupancy[,c("siteID","site","season","x","y")])

  #.GlobalEnv$model.input = model.input
  return(model.input)
}
