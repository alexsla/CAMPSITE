library(tidyverse)
library(phytools)
library(motmot)
library(RPANDA)
source("code/support functions.R")

clades <- str_remove(dir("empirical/data/"), ".csv")

sum.tab <- tibble(kurt = numeric(),
                  K = numeric(),
                  lambda = numeric(),
                  trait.model = character(),
                  rich = numeric(),
                  div.model = character(),
                  gamma = numeric())
test <- list()
for(i in 1:length(clades)){
  test[[1]] <- read.tree(paste("empirical/trees/", clades[[i]], ".tre", sep = ""))
  test[[2]] <- read_csv(paste("empirical/data/", clades[[i]], ".csv", sep = ""))$trait
  names(test[[2]]) <- read_csv(paste("empirical/data/", clades[[i]], ".csv", sep = ""))$species
  test[[2]] <- test[[2]][test[[1]]$tip.label]
  
  names(test) <- c("tree_extant", "tip_traits")
  
  sum.tab[i,1] <- e1071::kurtosis(test$tip_traits, na.rm = T)
  sum.tab[i,2] <- extractSignal(test, "K")
  sum.tab[i,3] <- extractSignal(test, "lambda")
  sum.tab[i,4] <- fitTraitModels(test)
  sum.tab[i,5] <- length(test$tree_extant$tip.label)
  sum.tab[i,6] <- fitDivModels(test)
  sum.tab[i,7] <- gammaStat(test$tree_extant)
}

sum.tab$clade <- clades
sum.tab %>% write_csv("output/empirical.csv")
