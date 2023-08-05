library(tidyverse)
library(phytools)
library(motmot)
library(RPANDA)
library(abc)
library(foreach)
library(doParallel)
library(harrypotter)
library(PBD)
source("code/support functions.R")

### alcids ###
alcid_data <- read_csv("empirical/data/Alcidae1.csv") %>% as.data.frame()
alcid_tree <- read.tree("empirical/trees/Alcidae1.tre")

rownames(alcid_data) <- alcid_data$species

sortedData <- sortTraitData(phy = alcid_tree, y = alcid_data, data.name = "trait", log.trait = F)

alcid_age <- round(max(nodeHeights(alcid_tree)), 0)

alcid_root <- transformPhylo.ML(y = sortedData$trait, phy = sortedData$phy, model = "bm")$root.state

# estimate BM
alcid_bm <- transformPhylo.ML(y = sortedData$trait, phy = sortedData$phy, model = "bm")$brownianVariance[[1]]

# estimate OU
alcid_ou <- transformPhylo.ML(y = sortedData$trait, phy = sortedData$phy, model = "ou")$Alpha[[1]]

# estimate MC
alcid_mc <- abs(fit_t_comp(sortedData$phy, setNames(as.numeric(t(sortedData$trait)), rownames(sortedData$trait)), model = "MC")$s)

brts <- branching.times(sortedData$phy)  # branching times
init_b <- 0.2  # speciation-initiation rate
init_mu_1 <- 0.05  # extinction rate of good species
init_la_1 <- 0.3 # speciation-completion rate
init_mu_2 <- 0.05  # extinction rate of incipient species  # not used

# The initial values of the parameters that must be optimized
initparsopt <- c(init_b, init_mu_1, init_la_1, init_mu_2)

# The first element of the branching times is the crown age (and not the stem age)
soc <- 2

# Conditioning on non-extinction of the phylogeny
# as I actively selected for a nice phylogeny
cond <- 1

# Give the likelihood of the phylogeny (instead of the likelihood of the branching times)
btorph <- 1

r <- pbd_ML(
  brts = brts,
  initparsopt = initparsopt, 
  exteq = F,
  soc = soc, 
  cond = cond,
  btorph = btorph,
  verbose = FALSE
)


# # run simulations for ABC:
set.seed(42)
sim.vals <- tibble(lambda1 = sapply(round(abs(rnorm(100, r$b, sd = 0.01)), 3), max, 0.001),
                   tau0 = sapply(round(abs(rnorm(100, r$lambda_1, sd = 0.01)), 3), max, 0.001),
                   mubg = sapply(round(abs(rnorm(100, r$mu_1, sd = 0.01)), 3), max, 0.001),
                   muibg = sapply(round(abs(rnorm(100, r$mu_2, sd = 0.01)), 3), max, 0.001),
                   bm.vals = sapply(round(abs(rnorm(100, alcid_bm*100, sd = 0.05)), 3), max, 0.001),
                   mc.vals = sapply(round(abs(rnorm(100, alcid_mc, sd = 0.05)), 3), max, 0.001),
                   ou.vals = sapply(round(abs(rnorm(100, alcid_ou, sd = 0.05)), 3), max, 0.001)) %>%
  unique()

runCAMPSITE.ABC <- function(nsim, lambda1, tau0, mubg, muibg, bm, mc, ou, path){
  source("code/sim_CAMPSITE.R", local = T)
  cores = detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  set.seed(42)
  sim <- foreach(i=1:nsim, .packages = c("ape"), .export = c("alcid_root", "alcid_age")) %dopar% {
    repeat {
      model <- sim_CAMPSITE(pars = c(lambda1,
                                     tau0,
                                     0.6,
                                     0.2,
                                     0.2,
                                     mubg,
                                     0.4,
                                     0.4,
                                     muibg,
                                     mc,
                                     mc,
                                     bm,
                                     20),
                            ou = list(c(alcid_root),
                                      c(ou),
                                      c(ou)),
                            bounds = c(-Inf,
                                       Inf),
                            root.value = alcid_root,
                            age.max = alcid_age,
                            step.size = 0.01,
                            full.sim = T,
                            plot = F
      )
      if (!isTRUE(model$gsp_extant$tree == c("process died")) && model$gsp_extant$tree$Nnode > 4) {break}
    }
    model
  }
  stopCluster(cl)
  saveRDS(sim, paste(path, lambda1, "_", tau0, "_", mubg, "_", muibg, "_", bm, "_", mc, "_", ou, ".rds", sep = ""))
  rm(sim)
}

# # simulate "BM":
# mapply(runCAMPSITE.ABC, rep(100, 100), sim.vals$lambda1, sim.vals$tau0, sim.vals$mubg, sim.vals$muibg, sim.vals$bm.vals, rep(0, 100), rep(0, 100), rep("output/emp_sims/alcids/", 100))
# 
# # simulate competition:
# mapply(runCAMPSITE.ABC, rep(100, 100), sim.vals$lambda1, sim.vals$tau0, sim.vals$mubg, sim.vals$muibg, sim.vals$bm.vals, sim.vals$mc.vals, rep(0, 100), rep("output/emp_sims/alcids/", 100))
# 
# # simulate selection:
# mapply(runCAMPSITE.ABC, rep(100, 100), sim.vals$lambda1, sim.vals$tau0, sim.vals$mubg, sim.vals$muibg, sim.vals$bm.vals, rep(0, 100), sim.vals$ou.vals, rep("output/emp_sims/alcids/", 100))
# 
# # simulate competition+selection:
# mapply(runCAMPSITE.ABC, rep(100, 100), sim.vals$lambda1, sim.vals$tau0, sim.vals$mubg, sim.vals$muibg, sim.vals$bm.vals, sim.vals$mc.vals, sim.vals$ou.vals, rep("output/emp_sims/alcids/", 100))
# 
# summarise
models <- list.files("output/emp_sims/alcids", full.names = T)
model.names <- list.files("output/emp_sims/alcids", full.names = F)
# 
# cores = detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
# 
# results <- foreach(i=1:length(models)) %dopar% {
#   model <- readRDS(models[[i]])
#   result <- lapply(model, sumModels)
#   result
# }
# stopCluster(cl)
# 
# names(results) <- stringr::str_replace(model.names, "\\.rds", "")
# 
# bm_ind <- which(sapply(str_split(str_remove(model.names, ".rds"), "_"), function(x) ifelse(length(which(x == 0)) == 2, TRUE, FALSE)))
# 
# comp_ind <- which(sapply(str_split(str_remove(model.names, ".rds"), "_"), function(x) ifelse(x[[6]] > 0 & x[[7]] == 0, TRUE, FALSE)))
# 
# selec_ind <- which(sapply(str_split(str_remove(model.names, ".rds"), "_"), function(x) ifelse(x[[6]] == 0 & x[[7]] > 0, TRUE, FALSE)))
# 
# comp_selec_ind <- which(sapply(str_split(str_remove(model.names, ".rds"), "_"), function(x) ifelse(length(which(x == 0)) == 0, TRUE, FALSE)))
# 
# 
# tip_traits_multi <- list()
# for (i in 1:length(results)){
#   tip_traits_multi[[i]] <- tipTraitPlot(lapply(results[[i]], "[[", "tip_traits"),
#                                         as.numeric(str_split(names(results)[i], "_")[[1]][6]),
#                                         as.numeric(str_split(names(results)[i], "_")[[1]][7])) %>%
#     add_column(bm = as.numeric(str_split(names(results)[i], "_")[[1]][5]),
#                lambda1 = as.numeric(str_split(names(results)[i], "_")[[1]][1]),
#                tau0 = as.numeric(str_split(names(results)[i], "_")[[1]][2]),
#                mubg = as.numeric(str_split(names(results)[i], "_")[[1]][3]),
#                muibg = as.numeric(str_split(names(results)[i], "_")[[1]][4]))
# }
# tip_traits_multi <- bind_rows(tip_traits_multi)
# 
# MNND_multi <- list()
# for (i in 1:length(results)){
#   MNND_multi[[i]] <- as.tibble(matrix(sapply(results[[i]], "[[", "MNND"), ncol = 100)) %>%
#     add_column(time = c(1:5000),
#                competition = as.numeric(str_split(names(results)[i], "_")[[1]][6]),
#                selection = as.numeric(str_split(names(results)[i], "_")[[1]][7]),
#                bm = as.numeric(str_split(names(results)[i], "_")[[1]][5]),
#                lambda1 = as.numeric(str_split(names(results)[i], "_")[[1]][1]),
#                tau0 = as.numeric(str_split(names(results)[i], "_")[[1]][2]),
#                mubg = as.numeric(str_split(names(results)[i], "_")[[1]][3]),
#                muibg = as.numeric(str_split(names(results)[i], "_")[[1]][4]))
# }
# MNND_multi <- bind_rows(MNND_multi)
# 
# 
# VNND_multi <- list()
# for (i in 1:length(results)){
#   VNND_multi[[i]] <- as.tibble(matrix(sapply(results[[i]], "[[", "VNND"), ncol = 100)) %>%
#     add_column(time = c(1:5000),
#                competition = as.numeric(str_split(names(results)[i], "_")[[1]][6]),
#                selection = as.numeric(str_split(names(results)[i], "_")[[1]][7]),
#                bm = as.numeric(str_split(names(results)[i], "_")[[1]][5]),
#                lambda1 = as.numeric(str_split(names(results)[i], "_")[[1]][1]),
#                tau0 = as.numeric(str_split(names(results)[i], "_")[[1]][2]),
#                mubg = as.numeric(str_split(names(results)[i], "_")[[1]][3]),
#                muibg = as.numeric(str_split(names(results)[i], "_")[[1]][4]))
# }
# VNND_multi <- bind_rows(VNND_multi)
# 
# k_multi <- list()
# lambda_multi <- list()
# gamma_multi <- list()
# for (i in 1:length(results)){
#   model <- results[[i]]
#   k_multi[[i]] <- numeric()
#   lambda_multi[[i]] <- numeric()
#   gamma_multi[[i]] <- numeric()
#   for (k in 1:length(model)){
#     k_multi[[i]][k] <- extractSignal(model[[k]], "K")
#     lambda_multi[[i]][k] <- extractSignal(model[[k]], "lambda")
#     gamma_multi[[i]][k] <- gammaStat(model[[k]]$tree_extant)
#   }
#   k_multi[[i]][101] <- as.numeric(str_split(names(results)[i], "_")[[1]][6])
#   k_multi[[i]][102] <- as.numeric(str_split(names(results)[i], "_")[[1]][7])
#   k_multi[[i]][103] <- as.numeric(str_split(names(results)[i], "_")[[1]][5])
#   k_multi[[i]][104] <- as.numeric(str_split(names(results)[i], "_")[[1]][1])
#   k_multi[[i]][105] <- as.numeric(str_split(names(results)[i], "_")[[1]][2])
#   k_multi[[i]][106] <- as.numeric(str_split(names(results)[i], "_")[[1]][3])
#   k_multi[[i]][107] <- as.numeric(str_split(names(results)[i], "_")[[1]][4])
#   lambda_multi[[i]][101] <- as.numeric(str_split(names(results)[i], "_")[[1]][6])
#   lambda_multi[[i]][102] <- as.numeric(str_split(names(results)[i], "_")[[1]][7])
#   lambda_multi[[i]][103] <- as.numeric(str_split(names(results)[i], "_")[[1]][5])
#   lambda_multi[[i]][104] <- as.numeric(str_split(names(results)[i], "_")[[1]][1])
#   lambda_multi[[i]][105] <- as.numeric(str_split(names(results)[i], "_")[[1]][2])
#   lambda_multi[[i]][106] <- as.numeric(str_split(names(results)[i], "_")[[1]][3])
#   lambda_multi[[i]][107] <- as.numeric(str_split(names(results)[i], "_")[[1]][4])
#   gamma_multi[[i]][101] <- as.numeric(str_split(names(results)[i], "_")[[1]][6])
#   gamma_multi[[i]][102] <- as.numeric(str_split(names(results)[i], "_")[[1]][7])
#   gamma_multi[[i]][103] <- as.numeric(str_split(names(results)[i], "_")[[1]][5])
#   gamma_multi[[i]][104] <- as.numeric(str_split(names(results)[i], "_")[[1]][1])
#   gamma_multi[[i]][105] <- as.numeric(str_split(names(results)[i], "_")[[1]][2])
#   gamma_multi[[i]][106] <- as.numeric(str_split(names(results)[i], "_")[[1]][3])
#   gamma_multi[[i]][107] <- as.numeric(str_split(names(results)[i], "_")[[1]][4])
# }
# 
# k_multi <- as_tibble(do.call(rbind, k_multi))
# names(k_multi) <- c(1:100, "competition", "selection", "bm", "lambda1", "tau0", "mubg", "muibg")
# k_multi$signal <- "K"
# 
# lambda_multi <- as_tibble(do.call(rbind, lambda_multi))
# names(lambda_multi) <- c(1:100, "competition", "selection", "bm", "lambda1", "tau0", "mubg", "muibg")
# lambda_multi$signal <- "lambda"
# 
# gamma_multi <- as_tibble(do.call(rbind, gamma_multi))
# names(gamma_multi) <- c(1:100, "competition", "selection", "bm", "lambda1", "tau0", "mubg", "muibg")
# gamma_multi$signal <- "gamma"
# 
# signal_multi <- bind_rows(k_multi, lambda_multi, gamma_multi)
# 
# richness_multi <- list()
# for (i in 1:length(results)){
#   richness_multi[[i]] <- as.tibble(sapply(results[[i]], "[[", "Nnode_extant")) %>%
#     add_column(competition = as.numeric(str_split(names(results)[i], "_")[[1]][6]),
#                selection = as.numeric(str_split(names(results)[i], "_")[[1]][7]),
#                bm = as.numeric(str_split(names(results)[i], "_")[[1]][5]),
#                lambda1 = as.numeric(str_split(names(results)[i], "_")[[1]][1]),
#                tau0 = as.numeric(str_split(names(results)[i], "_")[[1]][2]),
#                mubg = as.numeric(str_split(names(results)[i], "_")[[1]][3]),
#                muibg = as.numeric(str_split(names(results)[i], "_")[[1]][4]))
# }
# richness_multi <- bind_rows(richness_multi)
# 
# 
# traits_models <- tip_traits_multi %>%
#   gather("iteration", "var", -c(competition, selection, bm, lambda1, tau0, mubg, muibg)) %>%
#   group_by(bm, competition, selection, lambda1, tau0, mubg, muibg, iteration) %>%
#   summarise(mean = mean(var, na.rm = T),
#             sd = sd(var, na.rm = T),
#             kurt = e1071::kurtosis(var, na.rm = T)) %>%
#   left_join(MNND_multi %>%
#               filter(time == alcid_age*100) %>%
#               gather("iteration", "mnnd", -c(competition, selection, bm, lambda1, tau0, mubg, muibg)) %>%
#               mutate(iteration = str_remove(iteration, "V"))) %>%
#   left_join(VNND_multi %>%
#               filter(time == alcid_age*100) %>%
#               gather("iteration", "vnnd", -c(competition, selection, bm, lambda1, tau0, mubg, muibg)) %>%
#               mutate(iteration = str_remove(iteration, "V"))) %>%
#   left_join(signal_multi %>%
#               gather("iteration", "var", -c(competition, selection, bm, lambda1, tau0, mubg, muibg, signal)) %>%
#               spread(signal, var)) %>%
#   left_join(richness_multi %>%
#               add_column(iteration = as.character(rep(1:100, times = 400))) %>%
#               rename(richness = value)) %>%
#   mutate(model = case_when(competition == 0 & selection == 0 ~ "bm",
#                            competition > 0 & selection == 0 ~ "comp",
#                            competition == 0 & selection > 0 ~ "selec",
#                            T ~ "comp_selec"))
# 
# traits_models %>% write_csv("output/summary_alcids_empirical_models.csv")

traits_models <- read_csv("output/summary_alcids_empirical_models.csv") %>% as.data.frame()
models <- traits_models$model
traits <- traits_models[,10:17]
# 
# boxplot(traits[,"sd"]~models, main="SD of traits")
# boxplot(traits[,"kurt"]~models, main="Kurtosis of traits")
# boxplot(traits[,"mnnd"]~models, main="MNND of traits")
# boxplot(traits[,"vnnd"]~models, main="VNND of traits")
# boxplot(traits[,"gamma"]~models, main="Gamma")
# boxplot(traits[,"K"]~models, main="Blomberg's K")
# boxplot(traits[,"lambda"]~models, main="Pagel's lambda")
# boxplot(traits[,"richness"]~models, main="Species richness")

traits_models %>%
  pivot_longer(cols = 10:17) %>%
  mutate(name = case_when(name == "sd" ~ "SD of traits",
                          name == "kurt" ~ "Kurtosis of traits",
                          name == "mnnd" ~ "MNND of traits",
                          name == "vnnd" ~ "VNND of traits",
                          name == "gamma" ~ "Gamma",
                          name == "K" ~ "Blomberg's K",
                          name == "lambda" ~ "Pagel's lambda",
                          name == "richness" ~ "Species richness"),
         model = case_when(model == "bm" ~ "Brownian Motion",
                           model == "comp" ~ "Competition",
                           model == "selec" ~ "Selection",
                           model == "comp_selec" ~ "Competition + Selection")) %>%
  ggplot(aes(x = model, y = value, fill = model)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~name, scale = "free") +
  scale_fill_manual(values = hp(n = 4, option = "Ravenclaw")) +
  labs(x = "Model", y = "") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

set.seed(42)
cv.modsel <- cv4postpr(models, traits, nval=5, tol=.01, method="mnlogistic")
s <- summary(cv.modsel)

# plot(cv.modsel, names.arg=c("Brownian Motion", "Competition", "Competition + Selection", "Selection"))
cv.modsel$model.probs %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(true_model = rep(c("Brownian Motion", "Competition", "Competition + Selection", "Selection"), each = 5)) %>%
  rename(`Brownian Motion` = `tol0.01.bm`,
         `Competition` = `tol0.01.comp`,
         `Competition + Selection` = `tol0.01.comp_selec`,
         `Selection` = `tol0.01.selec`) %>%
  pivot_longer(cols = 1:4) %>%
  ggplot(aes(x = true_model, y = value, fill = name)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  scale_fill_manual(values = hp(n = 4, option = "Ravenclaw")) +
  labs(x = "Model", y = "") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

# model selection
stat.alcid <- data.frame(sd = sd(sortedData$trait),
                         kurt = e1071::kurtosis(sortedData$trait),
                         mnnd = MNND(sortedData$trait),
                         vnnd = VNND(sortedData$trait),
                         gamma = gammaStat(sortedData$phy),
                         K = as.numeric(phylosig(sortedData$phy, sortedData$trait, "K")),
                         lambda = phylosig(sortedData$phy, sortedData$trait, "lambda")$lambda,
                         richness = sortedData$phy$Nnode+1)

set.seed(42)
modsel.alcid <- postpr(stat.alcid[1,],
                       models,
                       traits,
                       tol=.025,
                       method="mnlogistic")

summary(modsel.alcid)

# goodness-of-fit

set.seed(42)
alcid.gfit.selec = gfit(target=stat.alcid[1,],
                        sumstat=traits[models=="selec",],
                        statistic=mean,
                        nb.replicate=100)

summary(alcid.gfit.selec)
plot(alcid.gfit.selec, main="Histogram under H0; Selection")

set.seed(42)
alcid.gfit.bm = gfit(target=stat.alcid[1,],
                     sumstat=traits[models=="bm",],
                     statistic=mean,
                     nb.replicate=100)

summary(alcid.gfit.bm)
plot(alcid.gfit.bm, main="Histogram under H0; Brownian Motion")

set.seed(42)
alcid.gfit.comp = gfit(target=stat.alcid[1,],
                       sumstat=traits[models=="comp",],
                       statistic=mean,
                       nb.replicate=100)

summary(alcid.gfit.comp)
plot(alcid.gfit.comp, main="Histogram under H0; Competition")

set.seed(42)
alcid.gfit.comp_selec = gfit(target=stat.alcid[1,],
                             sumstat=traits[models=="comp_selec",],
                             statistic=mean,
                             nb.replicate=100)

summary(alcid.gfit.comp_selec)
plot(alcid.gfit.comp_selec, main="Histogram under H0; Competition + Selection")

# PCA
res.prcomp = prcomp(traits, scale = T, center = T)

dc <- list()
for (i in unique(models)){
  alcid.kde <- kde2d(res.prcomp$x[which(models == i),1], res.prcomp$x[which(models == i),2])
  dx <- diff(alcid.kde$x[1:2])  # lifted from emdbook::HPDregionplot()
  dy <- diff(alcid.kde$y[1:2])
  sz <- sort(alcid.kde$z)
  c1 <- cumsum(sz) * dx * dy
  dimnames(alcid.kde$z) <- list(alcid.kde$x, alcid.kde$y)
  dc[[i]] <- reshape2::melt(alcid.kde$z)
  dc[[i]]$prob <- approx(sz,1-c1,dc[[i]]$value)$y
  dc[[i]]$model <- i
}

prob <- 0.9

ggplot(bind_rows(dc),aes(x=Var1,y=Var2))+
  geom_point(aes(x = PC1, y = PC2, color = model), alpha = .05, data = res.prcomp$x %>%
               as_tibble() %>%
               mutate(model = models)) +
  geom_contour(aes(z=prob, color=model),breaks=prob) +
  geom_point(aes(x = predict(res.prcomp, stat.alcid[1,])[1],
                 y = predict(res.prcomp, stat.alcid[1,])[2]), colour = "darkred", size = 1.5) +
  theme_bw() +
  scale_color_manual(values = hp(n = 4, option = "Ravenclaw"),
                     name = "Model",
                     labels = c("Brownian Motion", "Competition", "Competition + Selection", "Selection")) +
  ggtitle("90% envelope of the 2 Principal Components obtained with each model")

# cross-validation
stat.selec.sim <- subset(traits,subset=models=="selec")
par.selec.sim <- subset(traits_models[,1:7],subset=models=="selec")

set.seed(42)
cv.res.rej <- cv4abc(data.frame(alpha=par.selec.sim[, "selection"],
                                sigma2=par.selec.sim[, "bm"],
                                lambda1=par.selec.sim[, "lambda1"],
                                tau0=par.selec.sim[, "tau0"],
                                mubg=par.selec.sim[, "mubg"],
                                muibg=par.selec.sim[, "muibg"]),
                     stat.selec.sim,
                     nval=10,
                     tols=c(.005,.01,0.05),
                     method="rejection")
set.seed(42)
cv.res.reg <- cv4abc(data.frame(alpha=par.selec.sim[, "selection"],
                                sigma2=par.selec.sim[, "bm"],
                                lambda1=par.selec.sim[, "lambda1"],
                                tau0=par.selec.sim[, "tau0"],
                                mubg=par.selec.sim[, "mubg"],
                                muibg=par.selec.sim[, "muibg"]),
                     stat.selec.sim,
                     nval=10,
                     tols=c(.005,.01,0.05),
                     method="loclinear")
set.seed(42)
cv.res.nn <- cv4abc(data.frame(alpha=par.selec.sim[, "selection"],
                               sigma2=par.selec.sim[, "bm"],
                               lambda1=par.selec.sim[, "lambda1"],
                               tau0=par.selec.sim[, "tau0"],
                               mubg=par.selec.sim[, "mubg"],
                               muibg=par.selec.sim[, "muibg"]),
                    stat.selec.sim,
                    nval=10,
                    tols=c(.005,.01,0.05),
                    method="neuralnet")

summary(cv.res.rej)
summary(cv.res.reg)
summary(cv.res.nn)

# plot(cv.res.rej,caption="Rejection")
# plot(cv.res.reg,caption="Local linear regression")
# plot(cv.res.nn,caption="Neural network")

# parameter estimation

set.seed(42)
res_params <- abc(target=stat.alcid[1,],
                  param=data.frame(alpha=par.selec.sim[, "selection"],
                                   sigma2=par.selec.sim[, "bm"],
                                   lambda1=par.selec.sim[, "lambda1"],
                                   tau0=par.selec.sim[, "tau0"],
                                   mubg=par.selec.sim[, "mubg"],
                                   muibg=par.selec.sim[, "muibg"]),
                  sumstat=stat.selec.sim,
                  tol=0.05,
                  transf=c("log"),
                  method="rejection")
summary(res_params)
hist(res_params)

# posterior predictive tests:

# simulate:
post_params <- res_params$unadj.values %>% as_tibble() %>% group_by(alpha, sigma2, lambda1, tau0, mubg, muibg) %>% count() %>% rename(rep = n)

mapply(runCAMPSITE.ABC, post_params$rep, post_params$lambda1, post_params$tau0, post_params$mubg, post_params$muibg, post_params$sigma2, rep(0, dim(post_params)[1]), post_params$alpha, "output/emp_ppc/alcids/")

# summarise
models <- list.files("output/emp_ppc/alcids", full.names = T)
model.names <- list.files("output/emp_ppc/alcids", full.names = F)

cores = detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

results <- foreach(i=1:length(models)) %dopar% {
  model <- readRDS(models[[i]])
  result <- lapply(model, sumModels)
  result
}
stopCluster(cl)

names(results) <- stringr::str_replace(model.names, "\\.rds", "")

tip_traits_multi <- list()
for (i in 1:length(results)){
  tip_traits_multi[[i]] <- tipTraitPlot(lapply(results[[i]], "[[", "tip_traits"),
                                        as.numeric(str_split(names(results)[i], "_")[[1]][6]),
                                        as.numeric(str_split(names(results)[i], "_")[[1]][7])) %>%
    add_column(bm = as.numeric(str_split(names(results)[i], "_")[[1]][5]),
               lambda1 = as.numeric(str_split(names(results)[i], "_")[[1]][1]),
               tau0 = as.numeric(str_split(names(results)[i], "_")[[1]][2]),
               mubg = as.numeric(str_split(names(results)[i], "_")[[1]][3]),
               muibg = as.numeric(str_split(names(results)[i], "_")[[1]][4]))
}
tip_traits_multi <- bind_rows(tip_traits_multi)

MNND_multi <- list()
for (i in 1:length(results)){
  MNND_multi[[i]] <- as.tibble(matrix(sapply(results[[i]], "[[", "MNND"), ncol = length(results[[i]]))) %>%
    add_column(time = c(1:5000),
               competition = as.numeric(str_split(names(results)[i], "_")[[1]][6]),
               selection = as.numeric(str_split(names(results)[i], "_")[[1]][7]),
               bm = as.numeric(str_split(names(results)[i], "_")[[1]][5]),
               lambda1 = as.numeric(str_split(names(results)[i], "_")[[1]][1]),
               tau0 = as.numeric(str_split(names(results)[i], "_")[[1]][2]),
               mubg = as.numeric(str_split(names(results)[i], "_")[[1]][3]),
               muibg = as.numeric(str_split(names(results)[i], "_")[[1]][4]))
}
MNND_multi <- bind_rows(MNND_multi)


VNND_multi <- list()
for (i in 1:length(results)){
  VNND_multi[[i]] <- as.tibble(matrix(sapply(results[[i]], "[[", "VNND"), ncol = length(results[[i]]))) %>%
    add_column(time = c(1:5000),
               competition = as.numeric(str_split(names(results)[i], "_")[[1]][6]),
               selection = as.numeric(str_split(names(results)[i], "_")[[1]][7]),
               bm = as.numeric(str_split(names(results)[i], "_")[[1]][5]),
               lambda1 = as.numeric(str_split(names(results)[i], "_")[[1]][1]),
               tau0 = as.numeric(str_split(names(results)[i], "_")[[1]][2]),
               mubg = as.numeric(str_split(names(results)[i], "_")[[1]][3]),
               muibg = as.numeric(str_split(names(results)[i], "_")[[1]][4]))
}
VNND_multi <- bind_rows(VNND_multi)

k_multi <- list()
lambda_multi <- list()
gamma_multi <- list()
for (i in 1:length(results)){
  model <- results[[i]]
  k_multi[[i]] <- numeric()
  lambda_multi[[i]] <- numeric()
  gamma_multi[[i]] <- numeric()
  for (k in 1:length(model)){
    k_multi[[i]][k] <- extractSignal(model[[k]], "K")
    lambda_multi[[i]][k] <- extractSignal(model[[k]], "lambda")
    gamma_multi[[i]][k] <- gammaStat(model[[k]]$tree_extant)
  }
  k_multi[[i]][101] <- as.numeric(str_split(names(results)[i], "_")[[1]][6])
  k_multi[[i]][102] <- as.numeric(str_split(names(results)[i], "_")[[1]][7])
  k_multi[[i]][103] <- as.numeric(str_split(names(results)[i], "_")[[1]][5])
  k_multi[[i]][104] <- as.numeric(str_split(names(results)[i], "_")[[1]][1])
  k_multi[[i]][105] <- as.numeric(str_split(names(results)[i], "_")[[1]][2])
  k_multi[[i]][106] <- as.numeric(str_split(names(results)[i], "_")[[1]][3])
  k_multi[[i]][107] <- as.numeric(str_split(names(results)[i], "_")[[1]][4])
  lambda_multi[[i]][101] <- as.numeric(str_split(names(results)[i], "_")[[1]][6])
  lambda_multi[[i]][102] <- as.numeric(str_split(names(results)[i], "_")[[1]][7])
  lambda_multi[[i]][103] <- as.numeric(str_split(names(results)[i], "_")[[1]][5])
  lambda_multi[[i]][104] <- as.numeric(str_split(names(results)[i], "_")[[1]][1])
  lambda_multi[[i]][105] <- as.numeric(str_split(names(results)[i], "_")[[1]][2])
  lambda_multi[[i]][106] <- as.numeric(str_split(names(results)[i], "_")[[1]][3])
  lambda_multi[[i]][107] <- as.numeric(str_split(names(results)[i], "_")[[1]][4])
  gamma_multi[[i]][101] <- as.numeric(str_split(names(results)[i], "_")[[1]][6])
  gamma_multi[[i]][102] <- as.numeric(str_split(names(results)[i], "_")[[1]][7])
  gamma_multi[[i]][103] <- as.numeric(str_split(names(results)[i], "_")[[1]][5])
  gamma_multi[[i]][104] <- as.numeric(str_split(names(results)[i], "_")[[1]][1])
  gamma_multi[[i]][105] <- as.numeric(str_split(names(results)[i], "_")[[1]][2])
  gamma_multi[[i]][106] <- as.numeric(str_split(names(results)[i], "_")[[1]][3])
  gamma_multi[[i]][107] <- as.numeric(str_split(names(results)[i], "_")[[1]][4])
}

k_multi <- as_tibble(do.call(rbind, k_multi))
names(k_multi) <- c(1:100, "competition", "selection", "bm", "lambda1", "tau0", "mubg", "muibg")
k_multi$signal <- "K"

lambda_multi <- as_tibble(do.call(rbind, lambda_multi))
names(lambda_multi) <- c(1:100, "competition", "selection", "bm", "lambda1", "tau0", "mubg", "muibg")
lambda_multi$signal <- "lambda"

gamma_multi <- as_tibble(do.call(rbind, gamma_multi))
names(gamma_multi) <- c(1:100, "competition", "selection", "bm", "lambda1", "tau0", "mubg", "muibg")
gamma_multi$signal <- "gamma"

signal_multi <- bind_rows(k_multi, lambda_multi, gamma_multi)

richness_multi <- list()
for (i in 1:length(results)){
  richness_multi[[i]] <- as.tibble(sapply(results[[i]], "[[", "Nnode_extant")) %>%
    add_column(competition = as.numeric(str_split(names(results)[i], "_")[[1]][6]),
               selection = as.numeric(str_split(names(results)[i], "_")[[1]][7]),
               bm = as.numeric(str_split(names(results)[i], "_")[[1]][5]),
               lambda1 = as.numeric(str_split(names(results)[i], "_")[[1]][1]),
               tau0 = as.numeric(str_split(names(results)[i], "_")[[1]][2]),
               mubg = as.numeric(str_split(names(results)[i], "_")[[1]][3]),
               muibg = as.numeric(str_split(names(results)[i], "_")[[1]][4]))
}
richness_multi <- bind_rows(richness_multi)


ppc_models <- tip_traits_multi %>%
  gather("iteration", "var", -c(competition, selection, bm, lambda1, tau0, mubg, muibg)) %>%
  group_by(bm, competition, selection, lambda1, tau0, mubg, muibg, iteration) %>%
  summarise(mean = mean(var, na.rm = T),
            sd = sd(var, na.rm = T),
            kurt = e1071::kurtosis(var, na.rm = T)) %>%
  drop_na() %>%
  left_join(MNND_multi %>%
              filter(time == alcid_age) %>%
              gather("iteration", "mnnd", -c(competition, selection, bm, lambda1, tau0, mubg, muibg)) %>%
              mutate(iteration = str_remove(iteration, "V")) %>%
              drop_na()) %>%
  left_join(VNND_multi %>%
              filter(time == alcid_age) %>%
              gather("iteration", "vnnd", -c(competition, selection, bm, lambda1, tau0, mubg, muibg)) %>%
              mutate(iteration = str_remove(iteration, "V")) %>%
              drop_na()) %>%
  left_join(signal_multi %>%
              gather("iteration", "var", -c(competition, selection, bm, lambda1, tau0, mubg, muibg, signal)) %>%
              spread(signal, var) %>%
              drop_na()) %>%
  left_join(richness_multi %>%
              add_column(iteration = as.character(unlist(sapply(lengths(sapply(results, lengths)),
                                                                function(x) c(1:x))))) %>%
              rename(richness = value))

ppc_models %>% write_csv("output/summary_alcids_empirical_ppc.csv")

ppc_models <- read_csv("output/summary_alcids_empirical_ppc.csv")

ppc_models %>%
  pivot_longer(cols = 9:17,
               names_to = "parameter") %>%
  filter(!parameter == "mean") %>%
  ggplot(aes(x = value)) + 
  geom_histogram() +
  geom_vline(aes(xintercept = value),
             colour = "darkred",
             data = stat.alcid %>%
               as_tibble() %>%
               pivot_longer(cols = 1:8,
                            names_to = "parameter")) +
  facet_wrap(~parameter, scales = "free") +
  theme_bw()
