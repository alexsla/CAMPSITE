library(tidyverse)
library(phytools)
source("code/support functions.R")

clades <- str_remove(dir("empirical/data/"), ".csv")

sum.tab <- tibble(rich = numeric(),
                  sd = numeric(),
                  kurt = numeric(),
                  mnnd = numeric(),
                  vnnd = numeric(),
                  K = numeric(),
                  lambda = numeric(),
                  gamma = numeric())
test <- list()
for(i in 1:length(clades)){
  test[[1]] <- read.tree(paste("empirical/trees/", clades[[i]], ".tre", sep = ""))
  test[[2]] <- read_csv(paste("empirical/data/", clades[[i]], ".csv", sep = ""))$trait
  names(test[[2]]) <- read_csv(paste("empirical/data/", clades[[i]], ".csv", sep = ""))$species
  test[[2]] <- test[[2]][test[[1]]$tip.label]
  
  names(test) <- c("tree_extant", "tip_traits")
  
  sum.tab[i,1] <- length(test$tree_extant$tip.label)
  sum.tab[i,2] <- sd(test$tip_traits, na.rm = T)
  sum.tab[i,3] <- e1071::kurtosis(test$tip_traits, na.rm = T)
  sum.tab[i,4] <- MNND(test$tip_traits)
  sum.tab[i,5] <- VNND(test$tip_traits)
  sum.tab[i,6] <- extractSignal(test, "K")
  sum.tab[i,7] <- extractSignal(test, "lambda")
  sum.tab[i,8] <- gammaStat(test$tree_extant)
}

sum.tab$clade <- clades
sum.tab %>% write_csv("output/empirical.csv")

# run ABC
