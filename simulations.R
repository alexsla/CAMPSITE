library(RPANDA)
library(tidyverse)
library(ggpubr)
library(phylolm)
library(ggtree)
library(tidytree)
library(foreach)
library(doParallel)

runCAMPSITE <- function(comp, selec){
  source("code/sim_CAMPSITE.R", local = T)
  cores = detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  set.seed(42)
  sim <- foreach(i=1:100, .packages = c("ape")) %dopar% {
    repeat {
      model <- sim_CAMPSITE(pars = c(0.25,
                                   0.01,
                                   0.6,
                                   0.2,
                                   0.2,
                                   0.01,
                                   0.4,
                                   0.4,
                                   0.02,
                                   comp,
                                   comp,
                                   0.5,
                                   20),
                          ou = list(c(0),
                                    c(selec),
                                    c(selec)),
                          bounds = c(-Inf,
                                     Inf),
                          root.value = 0,
                          age.max = 50,
                          step.size = 0.01,
                          full.sim = T,
                          plot = F
      )
      if (!isTRUE(model$gsp_extant$tree == c("process died")) && model$gsp_extant$tree$Nnode > 4) {break}
    }
    model
  }
  stopCluster(cl)
  saveRDS(sim, paste("output/sims/comp", comp, "_selec", selec, ".rds", sep = ""))
  rm(sim)
}

selec.vec <- c(0, 0.01, 0.025, 0.05, 0.075, 0.1)
comp.vec <- c(0, 0.01, 0.025, 0.05, 0.075, 0.1)

params <- crossing(selec.vec, comp.vec)

mapply(runCAMPSITE, params$comp.vec, params$selec.vec)