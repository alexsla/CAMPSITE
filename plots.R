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

source("code/support functions.R")

models <- list.files("output/sims", full.names = T)
model.names <- list.files("output/sims", full.names = F)
comp <- as.numeric(stringr::str_match(model.names, "comp\\s*(.*?)\\s*_")[,2])
selec <- as.numeric(stringr::str_match(model.names, "selec\\s*(.*?)\\s*.rds")[,2])

results <- list()
for (i in 1:length(models)){
  model <- readRDS(models[[i]])
  results[[i]] <- lapply(model, sumModels)
}

names(results) <- stringr::str_replace(model.names, "\\.rds", "")

###TRAIT DISTRIBUTIONS AT TIPS
tip_traits_multi <- list()
for (i in 1:length(results)){
  tip_traits_multi[[i]] <- tipTraitPlot(lapply(results[[i]], "[[", "tip_traits"), comp[[i]], selec[[i]])
}
tip_traits_multi <- bind_rows(tip_traits_multi)

pdf("plots/trait_hist.pdf", useDingbats = F, height = 10, width = 10)
p_hist <- tip_traits_multi %>%
  #select(-model) %>%
  gather("iteration", "var", -c(competition, selection)) %>%
  ggplot(aes(x = var)) +
  #geom_density(aes(group = as.factor(iteration)), color = "light grey", alpha = 0.1) +
  geom_density(aes(color = as.factor(competition)), alpha = 0.7) +
  geom_vline(xintercept = 10, color = "red", linetype = "dashed") +
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
  geom_smooth(method = "loess", aes(group=1), colour = "dark red", show.legend = F) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
  scale_x_continuous(breaks = c(0, 0.01, 0.02, 0.04, 0.07, 0.1)) +
  labs(x = "Competition", y = "SD Trait at Tips") +
  #ylab("richness") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(cols = vars(selection),
             labeller = label_both)
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
  geom_smooth(method = "loess", aes(group=1), colour = "dark red", show.legend = F) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
  scale_x_continuous(breaks = c(0, 0.01, 0.02, 0.04, 0.07, 0.1)) +
  labs(x = "Competition", y = "Kurtosis Trait at Tips") +
  #ylab("richness") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(cols = vars(selection),
             labeller = label_both)
dev.off()

trait_legend <- get_legend(p_mean +
                             theme(legend.position = "bottom") +
                             guides(fill = guide_legend(nrow = 1,
                                                        title = "Competition")))

p_trait <- plot_grid(p_hist, p_sd, p_kurt, ncol = 1,
                     rel_heights = c(3,1,1),
                     labels = "AUTO")
p_trait <- plot_grid(p_trait, trait_legend, ncol = 1, rel_heights = c(19,1))

ggsave("plots/traits_composite.pdf",
       p_trait,
       width = 8.25,
       height = 11.75,
       units = "in",
       scale = 1.5)

###TRAIT DISPARITY THROUGH TIME
traits_multi <- list()
for (i in 1:length(results)){
  traits_multi[[i]] <- as.tibble(sapply(results[[i]], "[[", "traits")) %>%
    add_column(gen = c(1:5000), 
               competition = comp[[i]], 
               selection = selec[[i]])
}
traits_multi <- bind_rows(traits_multi)

pdf("plots/traits.pdf", useDingbats = F)
p_disp <- traits_multi %>%
  #select(-c(a, a3)) %>%
  gather("sim", "var", -c(gen, selection, competition)) %>%
  mutate(gen = cut(gen, 50, ordered_result = T, labels = c(1:50))) %>%
  group_by(sim, selection, gen, competition) %>%
  summarise(var = mean(var)) %>%
  ggplot(aes(x = as.numeric(gen), y = var, colour = as.factor(competition))) +
  geom_point(alpha = 0.1, size = 0.2, position = "jitter") +
  #geom_smooth(aes(group = as.factor(sim)), color = "light grey", alpha = 0.1, se = F) +
  stat_summary(fun = mean, size = 1, geom = "line") +
  scale_color_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  labs(x = "Time", y = "Trait disparity") +
  #scale_y_continuous(trans = "log", breaks = c(0, 5, 10, 20, 50)) +
  #ylab("richness") +
  theme_pubr() + 
  facet_grid(rows = vars(selection),
             scales = "free",
             labeller = label_both)
p_disp
dev.off()

###MNND
MNND_multi <- list()
for (i in 1:length(results)){
  MNND_multi[[i]] <- as.tibble(matrix(sapply(results[[i]], "[[", "MNND"), ncol = 100)) %>%
    add_column(time = c(1:5000),
               competition = comp[[i]],
               selection = selec[[i]])
}
MNND_multi <- bind_rows(MNND_multi)

pdf("plots/mnnd.pdf", useDingbats = F)
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
dev.off()

pdf("plots/mnnd_tips.pdf", useDingbats = F, width = 14)
p_mnnd_t <- MNND_multi %>%
  gather("sim", "MNND", -c(time, competition, selection)) %>%
  filter(time == 5000) %>%
  ggplot(aes(y = MNND, x = competition, fill = as.factor(competition))) +
  geom_boxplot() +
  geom_smooth(method = "loess", se=TRUE, aes(group=1), colour = "dark red") +
  scale_x_continuous(breaks = c(0, 0.01, 0.02, 0.04, 0.07, 0.1)) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
  labs(x = "Competition", y = "MNND") +
  #ylab("richness") +
  theme_pubr() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(cols = vars(selection),
             labeller = label_both)
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

pdf("plots/VNND.pdf", useDingbats = F)
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
dev.off()

pdf("plots/vnnd_tips.pdf", useDingbats = F, width = 14)
p_vnnd_t <- VNND_multi %>%
  gather("sim", "VNND", -c(time, competition, selection)) %>%
  filter(time == 5000) %>%
  ggplot(aes(y = VNND, x = competition, fill = as.factor(competition))) +
  geom_boxplot() +
  geom_smooth(method = "loess", se=TRUE, aes(group=1), colour = "dark red") +
  scale_x_continuous(breaks = c(0, 0.01, 0.02, 0.04, 0.07, 0.1)) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
  labs(x = "Competition", y = "VNND") +
  #ylab("richness") +
  theme_pubr() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(cols = vars(selection),
             labeller = label_both)
dev.off()

p_disparity <- plot_grid(p_mnnd, p_mnnd_t, p_vnnd, p_vnnd_t,
                         ncol = 1, rel_heights = c(3,1,3,1),
                         labels = "AUTO")
p_disparity <- plot_grid(p_disparity, trait_legend, ncol = 1, rel_heights = c(19,1))

ggsave("plots/disparity_composite.pdf",
       p_disparity,
       width = 8.25,
       height = 11.75,
       units = "in",
       scale = 1.5)

###RICHNESS AT TIPS
richness_multi <- list()
for (i in 1:length(results)){
  richness_multi[[i]] <- as.tibble(sapply(results[[i]], "[[", "Nnode_extant")) %>%
    add_column(competition = comp[[i]],
               selection = selec[[i]])
}
richness_multi <- bind_rows(richness_multi)

p_rich <- richness_multi %>%
  ggplot(aes(y = value, x = competition, fill = as.factor(competition))) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.data = mean_sdl, 
               geom = "pointrange", colour = "black") +
  scale_y_continuous(trans = "log10") +
  xlab("Competition") +
  ylab("Richness") +
  geom_smooth(method = "loess", se=TRUE, aes(group=1), colour = "dark red", show.legend = F) +
  scale_x_continuous(breaks = c(0, 0.01, 0.02, 0.04, 0.07, 0.1)) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
  facet_grid(cols = vars(selection),
             labeller = label_both) +
  theme_pubr() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  guides(fill = guide_legend(nrow = 1))

###LINEAGES THROUGH TIME
lin_df_multi <- list()
for (i in 1:length(results)){
  lin_df_multi[[i]] <- as_tibble(do.call(rbind, unname(Map(cbind, sim = 1:length(results[[i]]), lapply(results[[i]], "[[", "lineages")))),
                                 .name_repair = "unique") %>%
    add_column(competition = comp[[i]],
               selection = selec[[i]])
  names(lin_df_multi[[i]])[1:9] <- c("sim",
                                     "P_node",
                                     "D_node",
                                     "start_t",
                                     "end_t",
                                     "status",
                                     "comp_t",
                                     "sp_comp_t",
                                     "loc")
  
}
lin_df_multi <- bind_rows(lin_df_multi)

bl_df_multi <- list()
for (i in 1:length(results)){
  bl_df_multi[[i]] <- as_tibble(do.call(rbind, unname(Map(cbind, sim = 1:length(results[[i]]), lapply(lapply(results[[i]], "[[", "tree_fossil"), sumBL)))),
                                .name_repair = "unique") %>%
    add_column(competition = comp[[i]],
               selection = selec[[i]])
  names(bl_df_multi[[i]])[1:3] <- c("sim",
                                    "time_bin",
                                    "bl")
}
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

p_div <- lin_df_div %>%
  gather("process", "rate", c(extinction, speciation, diversification)) %>%
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

p_divrich <- plot_grid(p_div, p_rich, ncol = 1,
                       rel_heights = c(2,1),
                       labels = "AUTO")

ggsave("plots/div_composite.pdf",
       p_divrich,
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

ggsave("plots/turnover.pdf",
       p_turnover,
       width = 8.25,
       height = 8.25,
       units = "in",
       scale = 1.5)

pdf("plots/div2.pdf", useDingbats = F, height = 10, width = 10)
lin_df_div %>%
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
  facet_grid(cols = vars(selection),
             rows = vars(process),
             scales = "free_y")
dev.off()

p_div2 <- lin_df_div %>%
  mutate(speciation = speciation*bl,
         extinction = extinction*bl,
         diversification = speciation-extinction) %>%
  gather("process", "rate", c(extinction, speciation, diversification)) %>%
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
  facet_grid(cols = vars(selection),
             rows = vars(process),
             scales = "free_y")

p_div4 <- plot_grid(p_div2, p_div3, ncol = 1,
                    rel_heights = c(2,1),
                    labels = "AUTO")

ggsave("plots/div2_composite.pdf",
       p_div4,
       width = 8.25,
       height = 8.25,
       units = "in",
       scale = 1.5)

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

ggsave("plots/ltt.pdf",
       p_ltt,
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
  ggplot(aes(y = time_to_ext, x = start_t, colour = as.factor(competition))) +
  geom_point(size = .1, alpha = .2, show.legend = F) +
  geom_smooth(se = F) +
  scale_y_continuous(trans = "log10") +
  xlab("Start Time") +
  ylab("Time to Extinction") +
  scale_colour_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
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
  ggplot(aes(y = time_to_sp, x = start_t, colour = as.factor(competition))) +
  geom_point(size = .1, alpha = .2, show.legend = F) +
  geom_smooth(se = F) +
  scale_y_continuous(trans = "log10") +
  xlab("Start Time") +
  ylab("Speciation Duration") +
  scale_colour_manual(values = hp(n = 6, option = "Ravenclaw"), name = "Competition") +
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

ggsave("plots/times.pdf",
       p_time,
       width = 8.25,
       height = 8.25,
       units = "in",
       scale = 1.5)

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

pdf("plots/signal.pdf", useDingbats = F, height = 10, width = 14)
signal_multi %>% 
  gather("sim", "value", -c("competition", "selection", "signal")) %>%
  mutate(sim = as.numeric(sim)) %>%
  ggboxplot(x = "competition",
            y = "value",
            fill = "competition",
            bxp.errorbar = T) +
  facet_grid(cols = vars(selection),
             rows = vars(signal),
             scales = "free_y") +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw"))
dev.off()

###TREE BALANCE
sackin_multi <- list()
for (i in 1:length(results)){
  sackin_multi[[i]] <- tibble(sackin = sapply(lapply(lapply(results[[i]], "[[", "tree_extant"), FUN = as.treeshape), FUN = sackin, norm = "yule")) %>%
    add_column(competition = comp[[i]],
               selection = selec[[i]])
}
sackin_multi <- bind_rows(sackin_multi)

pdf("plots/balance.pdf", useDingbats = F, height = 10, width = 14)
sackin_multi %>% 
  ggboxplot(x = "competition",
            y = "sackin",
            fill = "competition",
            bxp.errorbar = T,
            xlab = "Competition",
            ylab = "Sackin Index") +
  scale_y_continuous(trans = "log10") +
  facet_grid(cols = vars(selection),
             scales = "free_y",
             labeller = label_both) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  guides(fill = guide_legend(nrow = 1,
                             title = "Competition"))
dev.off()

###TRAIT EVOLUTIONARY MODELS
cores = detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

trait_models_multi <- foreach(i=1:length(results), .packages = c("ape", "motmot", "RPANDA")) %dopar% {
  model <- results[[i]]
  model_multi <- character()
  for (k in 1:length(model)){
    write.table(c("started"), file=paste("output/trait models/comp_",comp[[i]],"selec_",selec[[i]],"model_",k,".txt", sep=""),
                sep="\t", col.names=F)
    model_multi[k] <- fitTraitModels(model[[k]])
    write.table(model_multi[k], file=paste("output/trait models/comp_",comp[[i]],"selec_",selec[[i]],"model_",k,".txt", sep=""),
                sep="\t", col.names=F)
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


pdf("plots/trait_models.pdf", useDingbats = F, width = 14)
p_trait_m <- trait_models_multi %>% 
  gather("sim", "model", -c("competition", "selection")) %>%
  group_by(model, competition, selection) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = factor(1), y = count, fill = model)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
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
dev.off()

p_trait_s <- signal_multi %>% 
  gather("sim", "value", -c("competition", "selection", "signal")) %>%
  mutate(sim = as.numeric(sim)) %>%
  filter(signal %in% c("K", "lambda")) %>%
  ggboxplot(x = "competition",
            y = "value",
            fill = "competition",
            bxp.errorbar = T,
            xlab = "Competition") +
  geom_hline(yintercept = 1,
             colour = "red",
             linetype = "dashed") +
  facet_grid(cols = vars(selection),
             rows = vars(signal),
             scales = "free_y",
             labeller = label_both) +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  guides(fill = guide_legend(nrow = 1,
                             title = "Competition"))

p_trait_signal <- plot_grid(p_trait_m, p_trait_s,
                            ncol = 1, rel_heights = c(2,1),
                            labels = "AUTO")

ggsave("plots/trait_signal_composite.pdf",
       p_trait_signal,
       width = 8,
       height = 11,
       units = "in",
       scale = 1.5)

###DIVERSIFICATION EVOLUTIONARY MODELS
cores = detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer

div_models_multi <- list()
registerDoParallel(cl)
for (i in 1:length(results)){
  model <- results[[i]]
  model_multi <- foreach(k = 1:length(model), .packages = c("ape", "motmot", "RPANDA")) %dopar% {
    best_mod <- fitDivModels(model[[k]])
    write.table(best_mod, file=paste("output/div models/comp_",comp[[i]],"selec_",selec[[i]],"model_",k,".txt", sep=""),
                sep="\t", col.names=F)
    best_mod
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


pdf("plots/div_models.pdf", useDingbats = F, width = 14)
p_div_m <- div_models_multi %>% 
  gather("sim", "model", -c("competition", "selection")) %>%
  group_by(model, competition, selection) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = factor(1), y = n, fill = model)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) +
  facet_grid(cols = vars(selection),
             rows = vars(competition),
             labeller = label_both) +
  scale_fill_manual(values = hp(n = 9, option = "LunaLovegood"), name = "Model") +
  theme_pubr() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        axis.title = element_blank())
dev.off()

p_div_s <- signal_multi %>% 
  gather("sim", "value", -c("competition", "selection", "signal")) %>%
  mutate(sim = as.numeric(sim)) %>%
  filter(signal == "gamma") %>%
  ggboxplot(x = "competition",
            y = "value",
            fill = "competition",
            bxp.errorbar = T,
            xlab = "Competition") +
  facet_grid(cols = vars(selection),
             scales = "free_y",
             labeller = label_both) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = hp(n = 6, option = "Ravenclaw")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  guides(fill = guide_legend(nrow = 1,
                             title = "Competition"))

p_div_signal <- plot_grid(p_div_m, p_div_s,
                          ncol = 1, rel_heights = c(2.5,1),
                          labels = "AUTO")

ggsave("plots/div_signal_composite.pdf",
       p_div_signal,
       width = 8,
       height = 11,
       units = "in",
       scale = 1.5)
