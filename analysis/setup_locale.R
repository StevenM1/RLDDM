hostname <- system('hostname', intern=TRUE)
if(hostname == 'tux11psy') {
  # tux
  dataDir <- '/home/stevenm/Documents/RLDDM/data'
  workDir <- '/home/stevenm/Documents/RLDDM'
  .libPaths(c('/home/stevenm/rpackages', .libPaths())) # for tux server
  lib.loc <- '/home/stevenm/rpackages'
} else if(Sys.getenv('USER')=='steven') {
  # running locale on MacBook
  dataDir <- '~/sync/PhDprojects/RLDDM/data'
  workDir <- '~/sync/PhDprojects/RLDDM'
  lib.loc <- .libPaths()
}