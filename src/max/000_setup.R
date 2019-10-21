# list.of.packages <- c("poibin","truncdist")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)


if (Sys.getenv()[["USER"]]=="domino"){
  
  dir = "/repos"
  
  # data_dir_1 = "/mnt/inputfiles_1/" # vcf files
  # data_dir_2 = "/mnt/inputfiles_2/" # sample to subject matching info

  # This is UCSC data?  
  data_dir_1 = "/stash/data/nonclin/DT-64/P-20160627-0001/" # vcf files
  # But this doesn't match UCSC data
  data_dir_2 = "/mnt/inputfiles_2/" # sample to subject matching info
  
  results_dir = "/mnt/results/"
  
} else {
  #dir = "/Users/golharr/workspace/P01465_CohortMatcher"
  #data_dir_1 = "/Users/golharr/workspace/P01465_CohortMatcher/data/P-20160627-0001/input/"
  #data_dir_2 = "/Users/golharr/workspace/P01465_CohortMatcher/data/P-20160627-0001/input2/"
  #results_dir = "/Users/golharr/workspace/P01465_CohortMatcher/data/P-20160627-0001/results/"
  #dir = "."
  code_dir = "/Users/golharr/workspace/NGS/cohort-matcher/src/max/"
  vcf_dir = "vcfs/"
  cache_dir = "cache/"
  #data_dir_1 = "."
  #data_dir_2 = "."
  results_dir = "results/"
}


require(tidyverse)
require(poisbinom)
require(parallel)
require(poibin)
require(reshape2)
require(truncdist)
require(grid)

