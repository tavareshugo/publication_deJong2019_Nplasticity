##
# Figure 8
##

#### setup ####

library(tidyverse)
library(patchwork)

# Change ggplot2 default aesthetics
theme_set(theme_bw() + 
            theme(panel.grid = element_blank(), 
                  text = element_text(size = 8),
                  plot.title = element_text(hjust = -0.07, size = 10, face = "bold", 
                                            margin = margin(t = -5, b = 1)),
                  plot.subtitle = element_text(size = 10, hjust = 0.5)))


#### Read data ####

# QTL scan results
qtl_results <- read_csv("../../data/qtl_magic/qtl2_scans_lm.csv")
qtl_thresh <- read_csv("../../data/qtl_magic/qtl2_scans_lm_perm.csv") %>% 
  filter(level == 0.05)

# Linear mixed model results
qtl_lmm <- read_csv("../../data/qtl_magic/qtl_scans_lmm.csv") %>% 
  select(-LOD_full) %>% 
  gather("test", "LOD", matches("LOD")) %>% 
  mutate(test = str_remove(test, "LOD_"))
qtl_lmm_thresh <- read_csv("../../data/qtl_magic/qtl2_scans_lmm_perm.csv") %>% 
  select(-LOD_full, -seed) %>% 
  summarise_all(funs(threshold = quantile(., probs = 0.95))) %>% 
  gather("test", "threshold") %>% 
  mutate(test = str_remove(test, "LOD_")) %>% 
  separate(test, c("test", "level"), sep = "_") 


#### Get peaks above thresholds ####

# For simple linear model
sig_lm <- qtl_results %>% 
  # Join with thresholds table 
  full_join(qtl_thresh, by = c("covariate", "trait", "nitrate")) %>% 
  # Keep those above threshold
  filter(LOD >= threshold) %>%
  # Take top LOD score for each sign level, test and chromosome 
  group_by(chromosome, trait, nitrate, covariate) %>% 
  arrange(desc(LOD)) %>% 
  slice(1)

# For mixed model
sig_lmm <- qtl_lmm %>% 
  # Join with thresholds table 
  full_join(qtl_lmm_thresh, by = "test") %>% 
  # Keep those above threshold
  filter(LOD >= threshold) %>%
  # Take top LOD score for each sign level, test and chromosome 
  group_by(level, test, chrom) %>% 
  arrange(desc(LOD)) %>% 
  slice(1) %>% 
  ungroup()


#### Panel A

p1 <- qtl_results %>% 
  filter(grepl("bolt", trait)) %>% 
  ggplot(aes(bp/1e6, LOD)) +
  geom_hline(data = qtl_thresh %>% filter(grepl("bolt", trait)),
             aes(yintercept = threshold), colour = "grey") +
  geom_line() +
  geom_point(data = sig_lm %>% filter(grepl("bolt", trait)), aes(colour = covariate)) +
  facet_grid(toupper(nitrate) ~ chromosome, scales = "free", space = "free") +
  labs(x = "Mb", tag = "A") +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 20, 4)) +
  scale_x_continuous(breaks = seq(0, 30, 10))


#### Panel B

p2 <- qtl_results %>% 
  filter(grepl("totalbr", trait)) %>% 
  ggplot(aes(bp/1e6, LOD)) +
  geom_hline(data = qtl_thresh %>% filter(grepl("totalbr", trait)),
             aes(yintercept = threshold, linetype = covariate), colour = "grey") +
  geom_line(aes(linetype = covariate)) +
  geom_point(data = sig_lm %>% filter(grepl("totalbr", trait)), aes(colour = covariate)) +
  facet_grid(factor(toupper(nitrate), levels = c("HN", "LN", "D")) ~ chromosome, scales = "free", space = "free") +
  labs(x = "Mb", tag = "B") +
  theme(panel.grid = element_blank()) +
  scale_linetype_manual(values = c(1, 3)) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 20, 4), expand = c(0, 1.5)) +
  scale_x_continuous(breaks = seq(0, 30, 10))


#### Panel C

# Build the plot
p3 <- qtl_lmm %>% 
  mutate(test = str_remove(test, "LOD_")) %>% 
  ggplot(aes(bp/1e6, LOD)) +
  geom_hline(data = qtl_lmm_thresh,
             aes(yintercept = threshold), colour = "grey") +
  geom_line() +
  geom_point(data = sig_lmm, colour = "brown", size = 2) +
  facet_grid(test ~ chrom, scales = "free", space = "free") +
  labs(x = "Mb", tag = "C") +
  scale_y_continuous(breaks = seq(0, 20, 4), expand = c(0, 1.5)) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  scale_shape_manual(values = c(3, 19)) +
  theme(panel.grid = element_blank(), 
        legend.position = "none")


#### Assemble panels

# This was further annotated in Inkscape
pdf("./figures/figure8.pdf", width = 7, height = 8)
p1 + p2 + p3 + plot_layout(ncol = 1, heights = c(2/9, 3/9, 2/9))
dev.off()


