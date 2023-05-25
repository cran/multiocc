test_that("Missing column site in detection gives error", {
  data(detection)
  data(occupancy)
  data(coords)
  detection <- detection[,-1]
  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  expect_error(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000),
    "detection does not have a column named 'site'")
})

test_that("Missing column season in detection gives error", {
  data(detection)
  data(occupancy)
  data(coords)
  detection <- detection[,-2]
  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  expect_error(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000),
    "detection does not have a column named 'season'")
})

test_that("Missing column survey in detection gives error", {
  data(detection)
  data(occupancy)
  data(coords)
  detection <- detection[,-3]
  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  expect_error(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000),
    "detection does not have a column named 'survey'")
})

test_that("Missing column site in occupancy gives error", {
  data(detection)
  data(occupancy)
  data(coords)
  occupancy <- occupancy[,-1]
  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  expect_error(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000),
    "occupancy does not have a column named 'site'")
})

test_that("Missing column season in occupancy gives error", {
  data(detection)
  data(occupancy)
  data(coords)
  occupancy <- occupancy[,-2]
  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  expect_error(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000),
    "occupancy does not have a column named 'season'")
})

test_that("Duplicated entries in detection gives error", {
  data(detection)
  data(occupancy)
  data(coords)
  detection <- rbind(detection[1,],detection)
  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  expect_error(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000),
    "Duplicated entries in detection")
})

test_that("Duplicated entries in occupancy gives error", {
  data(detection)
  data(occupancy)
  data(coords)
  occupancy <- rbind(occupancy[7,],occupancy)
  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  expect_error(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000),
    "Duplicated entries in occupancy")
})

test_that("Impermissible values in y", {
  data(detection)
  data(occupancy)
  data(coords)
  detection[2,4]<-2
  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  expect_error(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000),
    "Forbidden values in the detection matrix.  Only 1, 0, and NA are permissible entries for detection.  Counts should be converted to binary detections."
  )
})

test_that("If all NAs occur in y", {
  data(detection)
  data(occupancy)
  data(coords)
  DataNames <- list("species"=colnames(detection)[4:9],
                  "detection"=c("duration"),"occupancy"=c("forest","elev"))

  wh.na.detection <- which(is.na(detection$duration)==TRUE)
  wh.na.birds.1 <- which(is.na(detection$Great.tit)==TRUE)
  wh.na.birds.2 <- which(is.na(detection$Blue.tit)==TRUE)
  wh.na.birds.3 <- which(is.na(detection$Coal.tit)==TRUE)
  wh.na.birds.4 <- which(is.na(detection$Crested.tit)==TRUE)
  wh.na.birds.5 <- which(is.na(detection$Marsh.tit)==TRUE)
  wh.na.birds.6 <- which(is.na(detection$Willow.tit)==TRUE)
  wh.na.birds <- unique(c(wh.na.birds.1,wh.na.birds.2,wh.na.birds.3,wh.na.birds.4,wh.na.birds.5,wh.na.birds.6))

  wh.na.detection.only <- wh.na.detection[which(wh.na.detection %in% wh.na.birds==FALSE)]
  detection<-detection[-wh.na.detection.only,]

  occupancy <- occupancy[-which(occupancy$site == "Q125" & occupancy$season == 2),]

  expect_silent(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000)
  )
})

test_that("If there are NAs are in y, with additional NAs in W for non-missing y", {
  data(detection)
  data(occupancy)
  data(coords)

  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  expect_message(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000),
    "Warning: Rows in detection with missing covariates have been removed for purposes of fitting the model, but the site/season combination is retained in occupancy and therefore predictions will be outputted.")
})

test_that("A site/season combination is in detection, but not occupancy", {
  data(detection)
  data(occupancy)
  data(coords)
  wh.na.detection <- which(is.na(detection$duration)==TRUE)
  wh.na.birds <- which(is.na(detection$Great.tit)==TRUE)
  wh.na.detection.only <- wh.na.detection[which(wh.na.detection %in% wh.na.birds==FALSE)]
  detection<-detection[-wh.na.detection.only,]

  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  occupancy <- occupancy[-c(1,2670),]

  expect_message(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000),
    "Warning: Rows in detection have been removed corresponding to missing rows in occupancy.")
})

test_that("site/season combination is in occupancy, but not detection", {
  data(detection)
  data(occupancy)
  data(coords)
  wh.na.detection <- which(is.na(detection$duration)==TRUE)
  wh.na.birds <- which(is.na(detection$Great.tit)==TRUE)
  wh.na.detection.only <- wh.na.detection[which(wh.na.detection %in% wh.na.birds==FALSE)]
  detection<-detection[-wh.na.detection.only,]

  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  detection <- detection[-c(1,268,535),]

  expect_message(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000),
    "Warning: Rows with missing values in detection have been added to correspond to additional site/season combinations present in occupancy.")
})

test_that("remove site/season combinations that have missing occupancy covariate data", {
  data(detection)
  data(occupancy)
  data(coords)
  wh.na.detection <- which(is.na(detection$duration)==TRUE)
  wh.na.birds <- which(is.na(detection$Great.tit)==TRUE)
  wh.na.detection.only <- wh.na.detection[which(wh.na.detection %in% wh.na.birds==FALSE)]
  detection<-detection[-wh.na.detection.only,]

  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  occupancy[2,3]<-NA

  expect_message(
    multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000),
    "Warning: Rows in occupancy with missing covariates have been removed.  Corresponding rows in detection have also been removed.")
})



