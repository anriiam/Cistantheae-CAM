##########################CAM plots##########################
rm(list = ls()) # clear environment
setwd("~/Dropbox/Yale/Research/CAM")
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(ggpattern)
library(RColorBrewer)
set.seed(22)

#load data
cam <- read.csv("Fig2_Titrations_CAM_Diversity.csv")
View(cam)

cam2 <- cam %>%
  dplyr::select(PlantID,
                Genus,
                Species,
                Voucher_ID,
                Time,
                Treatment,
                H,
                Experiment) %>%
  dplyr::filter(!is.na(H))%>%
  pivot_wider(id_cols = c(PlantID, Genus, Species, Voucher_ID, Treatment, Experiment), names_from = Time , values_from = H)%>%
  mutate(deltaH = AM - PM)
View(cam2)

#summary tables
cam2 %>%
  group_by(Genus) %>% 
  summarise(number = length(PlantID))

#re-order how I want it
cam2$Species <- factor(cam2$Species, levels = c("Cistanthe_cachinelensis",
                                              "Cistanthe_crassifolia",
                                              "Cistanthe_cymosa",
                                              "Cistanthe_litoralis",
                                              "Montiopsis_gilliesii",
                                              "Montiopsis_umbellatum",
                                              "Philippiamra_calycina",
                                              "Calyptridium_parryi",
                                              "Thingia_ambigua"
                                              ))
cam2$Treatment <- factor(cam2$Treatment, levels = c("Watered","Drought"))

#filter out the NAs and extremely negative values
cam3 <- cam2 %>%
  dplyr::filter(deltaH > - 10)

#new column names
cam_new <- cam3 %>%
  mutate(New_Label = case_when(
    Voucher_ID == "AC121" ~ "T. ambigua",
    Voucher_ID == "LPH260" ~ "P. calycina",
    Voucher_ID == "LPH261" ~ "C. aff. cymosa 1",
    Voucher_ID == "LPH267" ~ "C. crassifolia",
    Voucher_ID == "LPH268" ~ "C. cachinalensis",
    Voucher_ID == "LPH280" ~ "C. aff. litoralis",
    Voucher_ID == "LPH300" ~ "C. aff. cymosa 2",
    Voucher_ID == "gilli" ~ "M. gilliesii",
    Voucher_ID == "parryi" ~ "C. parryi",
    Voucher_ID == "umbe" ~ "M. umbellata",
    TRUE ~ as.character(Voucher_ID)
  ))

cam_new$New_Label <- factor(cam_new$New_Label, levels = c("C. cachinalensis",
                                                "C. crassifolia",
                                                "C. aff. cymosa 1",
                                                "C. aff. cymosa 2",
                                                "C. aff. litoralis",
                                                "M. gilliesii",
                                                "M. umbellata",
                                                "P. calycina",
                                                "C. parryi",
                                                "T. ambigua"))

cam_new$Genus <- factor(cam_new$Genus, levels = c("Cistanthe",
                                                "Montiopsis",
                                                "Philippiamra",
                                                "Calyptridium",
                                                "Thingia"))


#plot
patterns <- c('Watered' = 'circle', 'Drought' = 'none')

p <- ggplot(cam_new, aes(x = New_Label, y = deltaH, fill = Genus, pattern = Treatment))
p + geom_boxplot_pattern(alpha=0.7, aes(colour = Genus),
                         outlier.shape = 2,
                         pattern_fill = 'white',
                         pattern_color = "white",
                         pattern_angle = 45, pattern_size = 0.5,
                         pattern_spacing = 0.01) +
  #geom_jitter(width = 0.1, size = 0.5, aes(colour = Treatment)) + 
  scale_fill_manual(values=c("#CC00CC", "#388e3c", "#d32f2f", "#fbc02d", "#6C3BAA")) +
  scale_pattern_manual(values = patterns) +
  scale_color_manual(values=c("#CC00CC", "#388e3c", "#d32f2f", "#fbc02d", "#6C3BAA")) +
  ylab("Titratable acidity (∆H µmol/g FW)" ) +
  stat_compare_means(method= "wilcox.test",
                     method.args = list(alternative = "greater"),
                     label.y = 75,
                     size = 4,
                     aes(group = Treatment, label = paste0("p = ", after_stat(p.format)))) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        #panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 13),
        legend.position = "none")


p + geom_boxplot(alpha=0.7, outlier.shape = 1, aes(colour = Genus)) +
  #geom_jitter(width = 0.1, size = 1.5, aes(colour = Genus)) + 
  scale_fill_manual(values=c("#CC00CC", "#388e3c", "#d32f2f", "#fbc02d", "#6C3BAA")) +
  scale_color_manual(values=c("#CC00CC", "#388e3c", "#d32f2f", "#fbc02d", "#6C3BAA")) +
  ylab("Titratable acidity (∆H µmol/g FW)" ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        legend.position = "right")

p1 <- ggplot(cam_new, aes(x = New_Label, y = deltaH, pattern = Treatment))
  #geom_jitter(width = 0.1, size = 1.5, aes(colour = Genus)) + 
p1 + geom_boxplot_pattern(alpha=0.7,
                         outlier.shape = 8,
                         pattern_fill = 'white',
                         pattern_color = "white",
                         pattern_angle = 45, pattern_size = 0.7,
                         pattern_spacing = 0.01,
                         fill = "grey" ) +
  #geom_jitter(width = 0.1, size = 0.5, aes(colour = Treatment)) + 
  scale_pattern_manual(values = patterns) +
  ylab("Titratable acidity (∆H µmol/g FW)" ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        legend.position = "right")



#### Stats: One-sided t-test
cam_cac <- cam_new %>%
  dplyr::filter(New_Label=="C. cachinalensis")
cam_cra <- cam_new %>%
  dplyr::filter(New_Label=="C. crassifolia")
cam_cy1 <- cam_new %>%
  dplyr::filter(New_Label=="C. aff. cymosa 1")
cam_cy2 <- cam_new %>%
  dplyr::filter(New_Label=="C. aff. cymosa 2")
cam_lit <- cam_new %>%
  dplyr::filter(New_Label=="C. aff. litoralis")
cam_gil <- cam_new %>%
  dplyr::filter(New_Label=="M. gilliesii")
cam_umb <- cam_new %>%
  dplyr::filter(New_Label=="M. umbellata")
cam_cal <- cam_new %>%
  dplyr::filter(New_Label=="P. calycina")
cam_par <- cam_new %>%
  dplyr::filter(New_Label=="C. parryi")
cam_amb <- cam_new %>%
  dplyr::filter(New_Label=="T. ambigua")
View(cam_cac)

res <- wilcox.test(deltaH ~ Treatment, data = cam_cra,
                   alternative = "less")
res

cam_d <- cam_cra %>%
  dplyr::filter(Treatment=="Drought")
resd <- wilcox.test(cam_d$deltaH, mu = 0, alternative = "greater")
resd

cam_w <- cam_cra %>%
  dplyr::filter(Treatment=="Watered")
resw <- wilcox.test(cam_w$deltaH, mu = 0, alternative = "greater")
resw

#more summary tables
cam_new %>%
  group_by(New_Label) %>% 
  summarise(number = length(PlantID))

average_summary <- cam_new %>%
  group_by(New_Label, Treatment) %>%
  summarize(
    mean_deltaH = mean(deltaH, na.rm = TRUE), # Calculate mean, removing NAs
    .groups = "drop" # Drops the grouping structure after summarizing
  )
print(average_summary)

