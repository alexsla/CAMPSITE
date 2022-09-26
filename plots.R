library(tidyverse)
library(ggpubr)
library(cowplot)
library(phytools)
library(ape)
library(geiger)
library(motmot)
library(RPANDA)
library(foreach)
library(doParallel)
library(harrypotter)
library(Hmisc)
library(paleotree)
library(adephylo)
library(apTreeshape)
library(waffle)
library(arbutus)

source("code/support functions.R")

models <- list.files("output/sims", full.names = T)
model.names <- list.files("output/sims", full.names = F)
comp <- as.numeric(stringr::str_match(model.names, "comp\\s*(.*?)\\s*_")[,2])
selec <- as.numeric(stringr::str_match(model.names, "selec\\s*(.*?)\\s*.rds")[,2])

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

###TRAIT DISTRIBUTIONS AT TIPS
tip_traits_multi <- list()
for (i in 1:length(results)){
  tip_traits_multi[[i]] <- tipTraitPlot(lapply(results[[i]], "[[", "tip_traits"), comp[[i]], selec[[i]])
}
tip_traits_multi <- bind_rows(tip_traits_multi)

pdf("plots/Figure S1.pdf", useDingbats = F, height = 10, width = 10)
p_hist <- tip_traits_multi %>%
  #select(-model) %>%
  gather("iteration", "var", -c(competition, selection)) %>%
  ggplot(aes(x = var)) +
  geom_density(aes(group = as.factor(iteration)), color = "light grey", alpha = 0.1) +
  geom_density(aes(color = as.factor(competition)), alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  # geom_vline(xintercept = -10, color = "red", linetype = "dashed") +
  # geom_vline(xintercept = 20, color = "blue") +
  # geom_vline(xintercept = -20, color = "blue") +
  scale_color_manual(values = hp(n = 6, option = "Ravenclaw")) +
  labs(x = "Trait value at tip") +
  #scale_y_continuous(trans = "log", breaks = c(0, 5, 10, 20, 50)) +
  #ylab("richness") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()) +
  facet_grid(rows = vars(competition),
             cols = vars(selection),
             labeller = label_both)
p_hist
dev.off()

pdf("plots/Figure 1.pdf", useDingbats = F, height = 5, width = 5)
p_hist <- tip_traits_multi %>%
  #select(-model) %>%
  gather("iteration", "var", -c(competition, selection)) %>%
  filter(competition %in% c(0, 0.05, 0.1),
         selection %in% c(0, 0.05, 0.1)) %>%
  ggplot(aes(x = var)) +
  geom_density(aes(group = as.factor(iteration)), color = "light grey", alpha = 0.1) +
  geom_density(aes(color = as.factor(competition)), alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  # geom_vline(xintercept = -10, color = "red", linetype = "dashed") +
  # geom_vline(xintercept = 20, color = "blue") +
  # geom_vline(xintercept = -20, color = "blue") +
  scale_color_manual(values = hp(n = 3, option = "Ravenclaw"), name = "Competition") +
  labs(x = "Trait value at tip") +
  #scale_y_continuous(trans = "log", breaks = c(0, 5, 10, 20, 50)) +
  #ylab("richness") +
  theme_pubr() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()) +
  facet_grid(rows = vars(competition),
             cols = vars(selection),
             labeller = label_both)
p_hist
dev.off()



pdf("plots/traits_tips_mean.pdf", useDingbats = F, width = 14)
p_mean <- tip_traits_multi %>%
  gather("iteration", "var", -c(competition, selection)) %>%
  group_by(competition, selection, iteration) %>%
  summarise(mean = mean(var, na.rm = T),
            sd = sd(var, na.rm = T),
            kurt = e1071::kurtosis(var, na.rm = T)) %>%
  ggplot(aes(y = mean, x = competition, fill = as.factor(competition))) +
  geom_boxplot() +
  geom_smooth(method = "loess", aes(group=1), colour = "dark red", show.legend = F) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
  labs(x = "Competition", y = "Mean Trait at Tips") +
  #ylab("richness") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(cols = vars(selection),
             labeller = label_both)
p_mean
dev.off()

pdf("plots/traits_tips_sd.pdf", useDingbats = F, width = 14)
p_sd <- tip_traits_multi %>%
  gather("iteration", "var", -c(competition, selection)) %>%
  group_by(competition, selection, iteration) %>%
  summarise(mean = mean(var, na.rm = T),
            sd = sd(var, na.rm = T),
            kurt = e1071::kurtosis(var, na.rm = T)) %>%
  ggplot(aes(y = sd, x = competition, fill = as.factor(competition))) +
  geom_boxplot() +
  geom_smooth(method = "loess", span = 1.5, aes(group=1), colour = "dark red", show.legend = F) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
  scale_x_continuous(breaks = c(0, 0.01, 0.025, 0.05, 0.075, 0.1)) +
  labs(x = "Competition", y = "SD Trait at Tips") +
  #ylab("richness") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(cols = vars(selection),
             labeller = label_both)
p_sd
dev.off()

pdf("plots/traits_tips_kurt.pdf", useDingbats = F, width = 14)
p_kurt <- tip_traits_multi %>%
  gather("iteration", "var", -c(competition, selection)) %>%
  group_by(competition, selection, iteration) %>%
  summarise(mean = mean(var, na.rm = T),
            sd = sd(var, na.rm = T),
            kurt = e1071::kurtosis(var, na.rm = T)) %>%
  ggplot(aes(y = kurt, x = competition, fill = as.factor(competition))) +
  geom_boxplot() +
  geom_smooth(method = "loess", span = 1.5, aes(group=1), colour = "dark red", show.legend = F) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
  scale_x_continuous(breaks = c(0, 0.01, 0.025, 0.05, 0.075, 0.1)) +
  labs(x = "Competition", y = "Kurtosis Trait at Tips") +
  #ylab("richness") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(cols = vars(selection),
             labeller = label_both)
p_kurt
dev.off()

trait_legend <- get_legend(p_mean +
                             theme(legend.position = "bottom") +
                             guides(fill = guide_legend(nrow = 1,
                                                        title = "Competition")))

p_trait <- plot_grid(p_sd, p_kurt, ncol = 1,
                     labels = "AUTO")
p_trait <- plot_grid(p_trait, trait_legend, ncol = 1, rel_heights = c(19,1))

ggsave("plots/Figure S2.pdf",
       p_trait,
       width = 4,
       height = 5,
       units = "in",
       scale = 2)

###TRAIT DISPARITY THROUGH TIME
traits_multi_var <- list()
for (i in 1:length(results)){
  traits_multi_var[[i]] <- as.tibble(sapply(results[[i]], "[[", "traits_var")) %>%
    add_column(gen = c(1:5000), 
               competition = comp[[i]], 
               selection = selec[[i]])
}
traits_multi_var <- bind_rows(traits_multi_var)

pdf("plots/traits_var.pdf", useDingbats = F)
p_disp_v <- traits_multi_var %>%
  gather("sim", "var", -c(gen, selection, competition)) %>%
  mutate(gen = cut(gen, 50, ordered_result = T, labels = c(1:50))) %>%
  group_by(selection, gen, competition) %>%
  summarise(var.mean = mean(var, na.rm = T), var.max = quantile(var, 0.975, na.rm = T), var.min = quantile(var, 0.025, na.rm = T)) %>%
  ggplot(aes(x = as.numeric(gen), y = var.mean, fill = as.factor(competition), colour = as.factor(competition))) +
  geom_ribbon(aes(ymax = var.max, ymin = var.min), alpha = 0.1, colour = NA) +
  stat_summary(aes(y = var.mean), fun = mean, size = 1, geom = "line") +
  scale_color_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  labs(x = "Time", y = "Trait disparity") +
  theme_pubr() + 
  facet_grid(rows = vars(selection),
             labeller = label_both)
p_disp_v
dev.off()

traits_multi_mean <- list()
for (i in 1:length(results)){
  traits_multi_mean[[i]] <- as.tibble(sapply(results[[i]], "[[", "traits_mean")) %>%
    add_column(gen = c(1:5000), 
               competition = comp[[i]], 
               selection = selec[[i]])
}
traits_multi_mean <- bind_rows(traits_multi_mean)

pdf("plots/traits_mean.pdf", useDingbats = F)
p_disp_m <- traits_multi_mean %>%
  gather("sim", "mean", -c(gen, selection, competition)) %>%
  mutate(gen = cut(gen, 50, ordered_result = T, labels = c(1:50))) %>%
  group_by(selection, gen, competition) %>%
  summarise(mean.mean = mean(mean, na.rm = T), mean.max = quantile(mean, 0.975, na.rm = T), mean.min = quantile(mean, 0.025, na.rm = T)) %>%
  ggplot(aes(x = as.numeric(gen), y = mean.mean, fill = as.factor(competition), colour = as.factor(competition))) +
  geom_ribbon(aes(ymax = mean.max, ymin = mean.min), alpha = 0.1, colour = NA) +
  stat_summary(aes(y = mean.mean), fun = mean, size = 1, geom = "line") +
  scale_color_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  labs(x = "Time", y = "Trait mean") +
  theme_pubr() + 
  facet_grid(rows = vars(selection),
             labeller = label_both)
p_disp_m
dev.off()

disp_legend <- get_legend(p_disp_v +
                             theme(legend.position = "bottom") +
                             guides(fill = guide_legend(nrow = 1,
                                                        title = "Competition")))

p_disp <- plot_grid(p_disp_m + theme(legend.position = "none"), p_disp_v + theme(legend.position = "none"), nrow = 1,
                     labels = "AUTO")
p_disp <- plot_grid(p_disp, disp_legend, ncol = 1, rel_heights = c(19,1))

ggsave("plots/Figure S3.pdf",
       p_disp,
       width = 8,
       height = 5,
       units = "in",
       scale = 2)

###MNND
MNND_multi <- list()
for (i in 1:length(results)){
  MNND_multi[[i]] <- as.tibble(matrix(sapply(results[[i]], "[[", "MNND"), ncol = 100)) %>%
    add_column(time = c(1:5000),
               competition = comp[[i]],
               selection = selec[[i]])
}
MNND_multi <- bind_rows(MNND_multi)

pdf("plots/mnnd_full.pdf", useDingbats = F)
p_mnnd <- MNND_multi %>%
  gather("sim", "MNND", -c(time, competition, selection)) %>%
  group_by(time, competition, selection) %>%
  summarise(MNND = mean(MNND)) %>%
  ggplot(aes(y = MNND, x = time/100, colour = as.factor(competition))) +
  geom_line() +
  scale_color_manual(values = hp(n = 6, option = "Ravenclaw")) +
  labs(x = "Time", y = "MNND") +
  #ylab("richness") +
  theme_pubr() + 
  theme(legend.position = "none") +
  facet_grid(rows = vars(selection),
             scales = "free",
             labeller = label_both)
p_mnnd
dev.off()

pdf("plots/mnnd.pdf", useDingbats = F)
p_mnnd <- MNND_multi %>%
  gather("sim", "MNND", -c(time, competition, selection)) %>%
  filter(competition %in% c(0, 0.05, 0.1),
         selection %in% c(0, 0.05, 0.1)) %>%
  group_by(time, competition, selection) %>%
  summarise(MNND = mean(MNND)) %>%
  ggplot(aes(y = MNND, x = time/100, colour = as.factor(competition))) +
  geom_line() +
  scale_color_manual(values = hp(n = 3, option = "Ravenclaw"), name = "Competition") +
  labs(x = "Time", y = "MNND") +
  #ylab("richness") +
  theme_pubr() + 
  theme(legend.position = "none") +
  facet_grid(rows = vars(selection),
             scales = "free",
             labeller = label_both)
p_mnnd
dev.off()

pdf("plots/mnnd_tips.pdf", useDingbats = F, width = 14)
p_mnnd_t <- MNND_multi %>%
  gather("sim", "MNND", -c(time, competition, selection)) %>%
  filter(time == 5000) %>%
  ggplot(aes(y = MNND, x = competition, fill = as.factor(competition))) +
  geom_boxplot() +
  geom_smooth(method = "loess", se=TRUE, aes(group=1), colour = "dark red") +
  scale_x_continuous(breaks = c(0, 0.01, 0.025, 0.05, 0.075, 0.1)) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
  labs(x = "Competition", y = "MNND") +
  #ylab("richness") +
  theme_pubr() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(cols = vars(selection),
             labeller = label_both)
p_mnnd_t
dev.off()

###VNND
VNND_multi <- list()
for (i in 1:length(results)){
  VNND_multi[[i]] <- as.tibble(matrix(sapply(results[[i]], "[[", "VNND"), ncol = 100)) %>%
    add_column(time = c(1:5000),
               competition = comp[[i]],
               selection = selec[[i]])
}
VNND_multi <- bind_rows(VNND_multi)

pdf("plots/VNND_full.pdf", useDingbats = F)
p_vnnd <- VNND_multi %>%
  gather("sim", "VNND", -c(time, competition, selection)) %>%
  group_by(time, competition, selection) %>%
  summarise(VNND = mean(VNND)) %>%
  ggplot(aes(y = VNND, x = time/100, colour = as.factor(competition))) +
  geom_line() +
  scale_color_manual(values = hp(n = 6, option = "Ravenclaw")) +
  labs(x = "Time", y = "VNND") +
  #ylab("richness") +
  theme_pubr() + 
  theme(legend.position = "none") +
  facet_grid(rows = vars(selection),
             scales = "free",
             labeller = label_both)
p_vnnd
dev.off()

pdf("plots/VNND.pdf", useDingbats = F)
p_vnnd <- VNND_multi %>%
  gather("sim", "VNND", -c(time, competition, selection)) %>%
  group_by(time, competition, selection) %>%
  filter(competition %in% c(0, 0.05, 0.1),
         selection %in% c(0, 0.05, 0.1)) %>%
  summarise(VNND = mean(VNND)) %>%
  ggplot(aes(y = VNND, x = time/100, colour = as.factor(competition))) +
  geom_line() +
  scale_color_manual(values = hp(n = 3, option = "Ravenclaw"), name = "Competition") +
  labs(x = "Time", y = "VNND") +
  #ylab("richness") +
  theme_pubr() + 
  theme(legend.position = "none") +
  facet_grid(rows = vars(selection),
             scales = "free",
             labeller = label_both)
p_vnnd
dev.off()

pdf("plots/vnnd_tips.pdf", useDingbats = F, width = 14)
p_vnnd_t <- VNND_multi %>%
  gather("sim", "VNND", -c(time, competition, selection)) %>%
  filter(time == 5000) %>%
  ggplot(aes(y = VNND, x = competition, fill = as.factor(competition))) +
  geom_boxplot() +
  geom_smooth(method = "loess", se=TRUE, aes(group=1), colour = "dark red") +
  scale_x_continuous(breaks = c(0, 0.01, 0.025, 0.05, 0.075, 0.1)) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
  labs(x = "Competition", y = "VNND") +
  #ylab("richness") +
  theme_pubr() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(cols = vars(selection),
             labeller = label_both)
p_vnnd_t
dev.off()

trait_legend <- get_legend(p_mnnd +
                             theme(legend.position = "bottom") +
                             guides(fill = guide_legend(nrow = 1,
                                                        title = "Competition")))

p_disparity <- plot_grid(p_mnnd, p_vnnd,
                         labels = "AUTO")
p_disparity <- plot_grid(p_disparity, trait_legend, ncol = 1, rel_heights = c(19,1))

ggsave("plots/Figure 2.pdf",
       p_disparity,
       width = 10,
       height = 5,
       units = "cm",
       scale = 2)

disparity_legend <- get_legend(p_mnnd_t +
                             theme(legend.position = "bottom") +
                             guides(fill = guide_legend(nrow = 1,
                                                        title = "Competition")))

p_disparity <- plot_grid(p_mnnd_t, p_vnnd_t, ncol = 1,
                     labels = "AUTO")
p_disparity <- plot_grid(p_disparity, disparity_legend, ncol = 1, rel_heights = c(19,1))

ggsave("plots/Figure S4.pdf",
       p_disparity,
       width = 4,
       height = 5,
       units = "in",
       scale = 2)

###PHYLOGENETIC SIGNAL
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
  k_multi[[i]][101] <- comp[[i]]
  k_multi[[i]][102] <- selec[[i]]
  lambda_multi[[i]][101] <- comp[[i]]
  lambda_multi[[i]][102] <- selec[[i]]
  gamma_multi[[i]][101] <- comp[[i]]
  gamma_multi[[i]][102] <- selec[[i]]
}

k_multi <- as_tibble(do.call(rbind, k_multi))
names(k_multi) <- c(1:100, "competition", "selection")
k_multi$signal <- "K"

lambda_multi <- as_tibble(do.call(rbind, lambda_multi))
names(lambda_multi) <- c(1:100, "competition", "selection")
lambda_multi$signal <- "lambda"

gamma_multi <- as_tibble(do.call(rbind, gamma_multi))
names(gamma_multi) <- c(1:100, "competition", "selection")
gamma_multi$signal <- "gamma"

signal_multi <- bind_rows(k_multi, lambda_multi, gamma_multi)

pdf("plots/Figure 3.pdf", useDingbats = F, height = 5, width = 5)
signal_multi %>% 
  gather("sim", "value", -c("competition", "selection", "signal")) %>%
  filter(competition %in% c(0, 0.05, 0.1),
         selection %in% c(0, 0.05, 0.1)) %>%
  mutate(sim = as.numeric(sim)) %>%
  ggboxplot(x = "competition",
            y = "value",
            fill = "competition",
            bxp.errorbar = T) +
  facet_grid(cols = vars(selection),
             rows = vars(signal),
             scales = "free_y",
             labeller = labeller(.cols = label_both)) +
  scale_fill_manual(values = hp(n = 3, option = "Ravenclaw"))
dev.off()

pdf("plots/gamma.pdf", useDingbats = F)
signal_multi %>% 
  gather("sim", "value", -c("competition", "selection", "signal")) %>%
  mutate(sim = as.numeric(sim)) %>%
  filter(signal == "gamma",
         !competition == 0.01,
         !selection == 0.01) %>%
  ggplot(aes(x = competition, y = selection, fill = value)) +
  geom_raster(interpolate = T) +
  xlab("Competition") +
  ylab ("Selection") +
  scale_fill_gradient2(low = "dark red", mid = "grey99", high = "dark blue", name = "gamma") +
  scale_x_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1)) +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

pdf("plots/K.pdf", useDingbats = F)
signal_multi %>% 
  gather("sim", "value", -c("competition", "selection", "signal")) %>%
  mutate(sim = as.numeric(sim)) %>%
  filter(signal == "K",
         !competition == 0.01,
         !selection == 0.01) %>%
  ggplot(aes(x = competition, y = selection, fill = value)) +
  geom_raster(interpolate = T) +
  xlab("Competition") +
  ylab ("Selection") +
  scale_fill_gradient2(low = "dark red", mid = "grey99", high = "dark blue", name = "K") +
  scale_x_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1)) +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

pdf("plots/lambda.pdf", useDingbats = F)
signal_multi %>% 
  gather("sim", "value", -c("competition", "selection", "signal")) %>%
  mutate(sim = as.numeric(sim)) %>%
  filter(signal == "lambda",
         !competition == 0.01,
         !selection == 0.01) %>%
  ggplot(aes(x = competition, y = selection, fill = value)) +
  geom_raster(interpolate = T) +
  xlab("Competition") +
  ylab ("Selection") +
  scale_fill_gradient2(low = "dark red", mid = "grey99", high = "dark blue", name = "lambda") +
  scale_x_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1)) +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075, 0.1)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

###TRAIT EVOLUTIONARY MODELS
cores = detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

trait_models_multi <- foreach(i=1:length(results), .packages = c("ape", "motmot", "RPANDA")) %dopar% {
  model <- results[[i]]
  model_multi <- character()
  for (k in 1:length(model)){
    if (model[[k]]$tree_extant$Nnode < 300) {
      write.table(c("started"), file=paste("output/trait models/comp_",comp[[i]],"selec_",selec[[i]],"model_",k,".txt", sep=""),
                  sep="\t", col.names=F)
      model_multi[k] <- fitTraitModels(model[[k]])
      write.table(model_multi[k], file=paste("output/trait models/comp_",comp[[i]],"selec_",selec[[i]],"model_",k,".txt", sep=""),
                  sep="\t", col.names=F)
    } else next
  }
  model_multi[101] <- comp[[i]]
  model_multi[102] <- selec[[i]]
  model_multi
}

stopCluster(cl)

trait_models_multi <- as_tibble(do.call(rbind, trait_models_multi))
names(trait_models_multi) <- c(1:100, "competition", "selection")
trait_models_multi <- trait_models_multi %>% mutate(competition = as.numeric(competition),
                                                    selection = as.numeric(selection))


pdf("plots/Figure S5.pdf", useDingbats = F, width = 14, height = 10)
p_trait_m <- trait_models_multi %>% 
  gather("sim", "model", -c("competition", "selection")) %>%
  group_by(model, competition, selection) %>%
  count() %>%
  drop_na() %>%
  ggplot(aes(values = n, fill = model)) +
  geom_waffle(n_rows = 10, colour = "white") +
  facet_grid(cols = vars(selection),
             rows = vars(competition),
             labeller = label_both) +
  scale_fill_manual(values = hp(n = 6, option = "ronweasley2"), name = "Model") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_pubr() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        axis.title = element_blank())
p_trait_m
dev.off()

###RICHNESS AT TIPS
richness_multi <- list()
for (i in 1:length(results)){
  richness_multi[[i]] <- as.tibble(sapply(results[[i]], "[[", "Nnode_extant")) %>%
    add_column(competition = comp[[i]],
               selection = selec[[i]])
}
richness_multi <- bind_rows(richness_multi)

pdf("plots/Figure S6.pdf", useDingbats = F, width = 10)
p_rich <- richness_multi %>%
  ggplot(aes(y = value, x = competition, fill = as.factor(competition))) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.data = mean_sdl, 
               geom = "pointrange", colour = "black") +
  scale_y_continuous(trans = "log10") +
  xlab("Competition") +
  ylab("Richness") +
  geom_smooth(method = "loess", se=TRUE, aes(group=1), colour = "dark red", show.legend = F) +
  scale_x_continuous(breaks = c(0, 0.01, 0.025, 0.05, 0.075, 0.1)) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  facet_grid(cols = vars(selection),
             labeller = label_both) +
  theme_pubr() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  guides(fill = guide_legend(nrow = 1))
p_rich
dev.off()

###LINEAGES THROUGH TIME
cores = detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

lin_df_multi <- foreach(i=1:length(results), .packages = c("tidyverse")) %dopar% {
  lin_df_multi_i <- as_tibble(do.call(rbind, unname(Map(cbind, sim = 1:length(results[[i]]), lapply(results[[i]], "[[", "lineages")))),
                              .name_repair = "unique") %>%
    add_column(competition = comp[[i]],
               selection = selec[[i]])
  names(lin_df_multi_i)[1:8] <- c("sim",
                                  "P_node",
                                  "D_node",
                                  "start_t",
                                  "end_t",
                                  "status",
                                  "comp_t",
                                  "sp_comp_t")
  lin_df_multi_i
}

bl_df_multi <- foreach(i=1:length(results), .packages = c("adephylo", "paleotree", "tidyverse")) %dopar% {
  bl_df_multi_i <- as_tibble(do.call(rbind, unname(Map(cbind, sim = 1:length(results[[i]]), lapply(lapply(results[[i]], "[[", "tree_fossil"), sumBL)))),
                             .name_repair = "unique") %>%
    add_column(competition = comp[[i]],
               selection = selec[[i]])
  names(bl_df_multi_i)[1:3] <- c("sim",
                                 "time_bin",
                                 "bl")
  bl_df_multi_i
}

stopCluster(cl)

lin_df_multi <- bind_rows(lin_df_multi)
bl_df_multi <- bind_rows(bl_df_multi)

lin_df_sp <- lin_df_multi %>%
  drop_na(sp_comp_t) %>%
  mutate(time_bin = cut(sp_comp_t, 50, ordered_result = T, labels = c(1:50))) %>% 
  group_by(competition, selection, sim, time_bin) %>%
  summarise(rate = n()) %>%
  add_column(process = "speciation")

lin_df_ext <- lin_df_multi %>%
  drop_na(sp_comp_t) %>%
  filter(status == -2) %>%
  mutate(time_bin = cut(comp_t, 50, ordered_result = T, labels = c(1:50))) %>% 
  group_by(competition, selection, sim, time_bin) %>%
  summarise(rate = n()) %>%
  add_column(process = "extinction")

lin_df_div <- rbind(lin_df_sp, lin_df_ext) %>%
  spread(process, rate) %>%
  mutate(speciation = replace_na(speciation, 0),
         extinction = replace_na(extinction, 0),
         time_bin = as.numeric(time_bin)) %>%
  left_join(bl_df_multi, by = c("sim", "time_bin", "competition", "selection")) %>%
  mutate(speciation = speciation/bl,
         extinction = extinction/bl,
         diversification = speciation-extinction,
         turnover = replace_na(na_if(extinction/speciation, Inf),0))

p_div2 <- lin_df_div %>%
  mutate(speciation = speciation*bl,
         extinction = extinction*bl,
         diversification = speciation-extinction) %>%
  gather("process", "rate", c(extinction, speciation, diversification)) %>%
  filter(competition %in% c(0, 0.05, 0.1),
         selection %in% c(0, 0.05, 0.1)) %>%
  group_by(time_bin, competition, selection, process) %>%
  summarise(rate = mean(rate)) %>%
  ggplot(aes(x = as.numeric(time_bin), y = rate, color = as.factor(process))) +
  #geom_smooth(aes(group = as.factor(sim)), color = "light grey", alpha = 0.1, se = F) +
  #geom_point(color = "light grey", size = 0.1, alpha = 0.1) +
  geom_line() +
  # scale_y_continuous(trans = "log", breaks = c(0, 1, 5, 10, 50)) +
  scale_color_manual(values = hp(n = 3, option = "DracoMalfoy"), name = "Process") +
  labs(x = "Time", y = "Rate") +
  theme_pubr() +
  theme(axis.text.x = element_blank()) +
  facet_grid(rows = vars(competition),
             cols = vars(selection),
             labeller = label_both)

p_div3 <- lin_df_div %>%
  gather("process", "rate", c(extinction, speciation, diversification)) %>%
  filter(competition %in% c(0, 0.05, 0.1),
         selection %in% c(0, 0.05, 0.1)) %>%
  group_by(time_bin, competition, selection, process) %>%
  summarise(rate = mean(rate)) %>%
  ggplot(aes(x = as.numeric(time_bin), y = rate, color = as.factor(process))) +
  #geom_smooth(aes(group = as.factor(sim)), color = "light grey", alpha = 0.1, se = F) +
  #geom_point(color = "light grey", size = 0.1, alpha = 0.1) +
  geom_line() +
  # scale_y_continuous(trans = "log", breaks = c(0, 1, 5, 10, 50)) +
  scale_color_manual(values = hp(n = 3, option = "DracoMalfoy"), name = "Process") +
  labs(x = "Time", y = "Rate") +
  theme_pubr() +
  theme(axis.text.x = element_blank()) +
  facet_grid(rows = vars(competition),
             cols = vars(selection),
             labeller = label_both)

div_legend <- get_legend(p_div2 +
                             theme(legend.position = "bottom") +
                             guides(fill = guide_legend(nrow = 1,
                                                        title = "Competition")))

p_div <- plot_grid(p_div2 + theme(legend.position = "none"), p_div3 + theme(legend.position = "none"),
                         labels = "AUTO")
p_div <- plot_grid(p_div, div_legend, ncol = 1, rel_heights = c(19,1))

ggsave("plots/Figure 4.pdf",
       p_div,
       width = 10,
       height = 5,
       units = "cm",
       scale = 2)

p_ltt <- lin_df_div %>%
  mutate(speciation = speciation*bl,
         extinction = extinction*bl,
         diversification = speciation-extinction) %>%
  group_by(sim, time_bin, competition, selection) %>%
  summarise(richness = cumsum(diversification)) %>%
  # group_by(time_bin, competition, selection) %>%
  # summarise(richness = mean(richness)) %>%
  ggplot(aes(x = as.numeric(time_bin), y = richness, colour = as.factor(competition))) +
  geom_smooth(aes(group = as.factor(sim)), colour = "light grey", se = F, alpha = 0.1, size = .2) +
  geom_smooth(se = F) +
  scale_y_continuous(trans = "log", breaks = c(0, 1, 10, 100)) +
  scale_color_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  labs(x = "Time", y = "Richness") +
  theme_pubr() +
  theme(axis.text.x = element_blank()) +
  facet_grid(cols = vars(selection),
             rows = vars(competition),
             labeller = label_both)

ggsave("plots/Figure S7.pdf",
       p_ltt,
       width = 8.25,
       height = 8.25,
       units = "in",
       scale = 1.5)

p_turnover <- lin_df_div %>%
  group_by(time_bin, competition, selection) %>%
  summarise(rate = mean(turnover)) %>%
  ggplot(aes(x = as.numeric(time_bin), y = rate)) +
  #geom_smooth(aes(group = as.factor(sim)), color = "light grey", alpha = 0.1, se = F) +
  #geom_point(color = "light grey", size = 0.1, alpha = 0.1) +
  geom_line(colour = hp(n = 1, option = "HarryPotter")) +
  # scale_y_continuous(trans = "log", breaks = c(0, 1, 5, 10, 50)) +
  #scale_color_manual(values = hp(n = 3, option = "DracoMalfoy"), name = "Process") +
  labs(x = "Time", y = "Turnover Rate") +
  theme_pubr() +
  theme(axis.text.x = element_blank()) +
  facet_grid(rows = vars(competition),
             cols = vars(selection),
             labeller = label_both)

ggsave("plots/Figure S8.pdf",
       p_turnover,
       width = 8.25,
       height = 8.25,
       units = "in",
       scale = 1.5)

p_ext <- lin_df_multi %>%
  drop_na(sp_comp_t) %>%
  mutate(time_to_ext = end_t - sp_comp_t,
         time_to_sp = sp_comp_t - start_t) %>%
  filter(status == -2,
         start_t <= 40) %>%
  ggplot(aes(y = time_to_ext, x = start_t, colour = as.factor(competition), fill = as.factor(competition))) +
  # geom_point(size = .1, alpha = .2, show.legend = F) +
  geom_smooth(se = T) +
  scale_y_continuous(trans = "log10") +
  xlab("Start Time") +
  ylab("Time to Extinction") +
  scale_colour_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  facet_grid(cols = vars(selection),
             labeller = label_both) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  guides(fill = guide_legend(nrow = 1))

p_sp <- lin_df_multi %>%
  drop_na(sp_comp_t) %>%
  mutate(time_to_ext = end_t - sp_comp_t,
         time_to_sp = sp_comp_t - start_t) %>%
  filter(start_t <= 40,
         time_to_sp > 0) %>%
  ggplot(aes(y = time_to_sp, x = start_t, colour = as.factor(competition), fill = as.factor(competition))) +
  # geom_point(size = .1, alpha = .2, show.legend = F) +
  geom_smooth(se = T) +
  scale_y_continuous(trans = "log10") +
  xlab("Start Time") +
  ylab("Speciation Duration") +
  scale_colour_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  facet_grid(cols = vars(selection),
             labeller = label_both) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  guides(fill = guide_legend(nrow = 1))

time_legend <- get_legend(p_ext +
                            theme(legend.position = "bottom") +
                            guides(fill = guide_legend(nrow = 1,
                                                       title = "Competition")))

p_time <- plot_grid(p_sp, p_ext,
                    ncol = 1,
                    labels = "AUTO")

p_time <- plot_grid(p_time, time_legend,
                    ncol = 1, rel_heights = c(19,1))

ggsave("plots/Figure S9.pdf",
       p_time,
       width = 8.25,
       height = 8.25,
       units = "in",
       scale = 1.5)

pdf("plots/Figure S10.pdf", useDingbats = F, width = 10)
lin_df_div %>%
  mutate(speciation = speciation*bl,
         extinction = extinction*bl,
         diversification = speciation-extinction) %>%
  gather("process", "rate", c(extinction, speciation, diversification)) %>%
  group_by(time_bin, competition, selection, process) %>%
  summarise(rate = mean(rate)) %>%
  ggplot(aes(x = as.numeric(time_bin), y = rate, color = as.factor(competition))) +
  #geom_smooth(aes(group = as.factor(sim)), color = "light grey", alpha = 0.1, se = F) +
  #geom_point(color = "light grey", size = 0.1, alpha = 0.1) +
  geom_line() +
  # scale_y_continuous(trans = "log", breaks = c(0, 1, 5, 10, 50)) +
  scale_color_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  labs(x = "Time", y = "Rate") +
  theme_pubr() +
  theme(axis.text.x = element_blank()) +
  facet_grid(rows = vars(process),
             cols = vars(selection),
             labeller = label_both)
dev.off()


###DIVERSIFICATION EVOLUTIONARY MODELS
cores = detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer

div_models_multi <- list()
registerDoParallel(cl)
for (i in 1:length(results)){
  model <- results[[i]]
  model_multi <- foreach(k = 1:length(model), .packages = c("ape", "motmot", "RPANDA")) %dopar% {
    if (model[[k]]$tree_extant$Nnode < 300) {
      best_mod <- fitDivModels(model[[k]])
      write.table(best_mod, file=paste("output/div models/comp_",comp[[i]],"selec_",selec[[i]],"model_",k,".txt", sep=""),
                  sep="\t", col.names=F)
      best_mod
    } else next
  }
  model_multi <- as.vector(unlist(model_multi, use.names=FALSE))
  model_multi[101] <- comp[[i]]
  model_multi[102] <- selec[[i]]
  div_models_multi[[i]] <- model_multi
  write.table(div_models_multi[[i]], file=paste("output/div models/comp_",comp[[i]],"selec_",selec[[i]],".txt", sep=""),
              sep="\t", col.names=F)
}
stopCluster(cl)

div_models_multi <- as_tibble(do.call(rbind, div_models_multi))
names(div_models_multi) <- c(1:100, "competition", "selection")
div_models_multi <- div_models_multi %>% mutate(competition = as.numeric(competition),
                                                selection = as.numeric(selection))


pdf("plots/Figure S11.pdf", useDingbats = F, width = 14, height = 10)
p_div_m <- div_models_multi %>% 
  gather("sim", "model", -c("competition", "selection")) %>%
  group_by(model, competition, selection) %>%
  count() %>%
  drop_na() %>%
  ggplot(aes(values = n, fill = model)) +
  geom_waffle(n_rows = 10, colour = "white") +
  facet_grid(cols = vars(selection),
             rows = vars(competition),
             labeller = label_both) +
  scale_fill_manual(values = hp(n = 9, option = "LunaLovegood"), name = "Model") +
  theme_pubr() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        axis.title = element_blank())
p_div_m
dev.off()
