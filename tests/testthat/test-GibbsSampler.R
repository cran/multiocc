test_that("Incorrect param2keep specification", {
  data(detection)
  data(occupancy)
  data(coords)
  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  model.input <- multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000)

  expect_error(
    GibbsSampler(M.iter=10, M.burn=1, M.thin=1, model.input, q=10, sv=FALSE, param2keep=c("alpha","beta","gamma","rho","sigma","pi")),
    "param2keep must be a character vector, whose only permissible entries are `alpha`, `beta`, `gamma`, `rho`, `sigma`, `psi`, `z`, `p`, and `loglik`")
})

test_that("Trying to specify too many basis functions", {
  data(detection)
  data(occupancy)
  data(coords)
  DataNames <- list("species"=colnames(detection)[4:9],
                    "detection"=c("duration"),"occupancy"=c("forest","elev"))
  model.input <- multioccbuild(detection, occupancy, coords, DataNames, threshold = 15000)

  expect_error(
    GibbsSampler(M.iter=10, M.burn=1, M.thin=1, model.input, q=300, sv=FALSE, param2keep=c("alpha","beta","gamma","rho","sigma","psi")),
    "Number of basis functions cannot be larger than the smallest number of sites in any season")
})
