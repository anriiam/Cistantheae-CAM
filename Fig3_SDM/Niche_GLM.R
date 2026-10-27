##########################  Species distribution modeling using GLM ########################## 
# Anri Chomentowska, 2025

#############  Set up ############# 
# I am working on R 4.4.0
rm(list = ls()) # clear environment

# CHANGE THIS to your working directory
setwd("~/Dropbox/Yale/Research/Map")

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(geodata)
library(terra)
library(sf)
library(dismo)
library(brglm2)
#library(maxnet) #I don't think this is used but just in case leaving this

# Load points & rasters
# Make sure these files are under the working directory
df <- read.csv("2025_data_elev_bioclim.csv")
pred <- rast(list.files("bioclim_output/", pattern="\\.asc$", full.names = TRUE)) #download bioclim from https://www.worldclim.org/data/bioclim.html

# Filter for species of interest
# 1. North American un
df <- df %>%
  dplyr::filter(species=="Calyptridium parryi" | species=="Thingia ambigua")

# 2. South American run
df <- df %>%
  dplyr::filter(species=="Montiopsis gilliesii" | species=="Cistanthe cachinalensis" | species=="Philippiamra calycina" | species=="Cistanthe cymosa" | species=="Cistanthe longiscapa")

df <- df %>%
  mutate(species = ifelse(species == "Cistanthe cymosa",
                          "C. longiscapa complex",
                          species))

df <- df %>%
  mutate(species = ifelse(species == "Cistanthe longiscapa",
                          "C. longiscapa complex",
                          species))

View(df)
table(df$species)

##### Determine geographic extent of our data
max_lat <- ceiling(max(df$latitude))
min_lat <- floor(min(df$latitude))
max_lon <- ceiling(max(df$longitude))
min_lon <- floor(min(df$longitude))
# Store boundaries in a single extent object
geographic_extent <- ext(x = c(min_lon, max_lon, min_lat, max_lat))

# We will just map things as a sanity check
# Download data with geodata's world function to use for base map
world_map <- world(resolution = 3)

# Crop the map to our area of interest
my_map <- crop(x = world_map, y = geographic_extent)

# Plot the base map
plot(my_map,
     axes = TRUE, 
     col = "grey95")

# Add the points for individual observations
points(x = df$longitude, 
       y = df$latitude, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)

# Make an extent that is 25% larger
sample_extent <- geographic_extent * 1.25

# Crop bioclim data to desired extent
bioclim_data <- crop(x = pred, y = sample_extent)

# Plot the first of the bioclim variables to check on cropping
plot(bioclim_data[[1]])


############# SDM using presence/absence data under GLM ############# 

##### Set up
# Build M (accessible area) as convex hull + buffer (in km)
make_M <- function(xy, crs_ll = "EPSG:4326", buf_km = 200){
  s <- st_as_sf(xy, coords = c("longitude","latitude"), crs = crs_ll)
  s_eq <- st_transform(s, "EPSG:6933")
  hull <- st_convex_hull(st_union(s_eq))
  M_eq <- st_buffer(hull, dist = buf_km * 1000)
  st_transform(M_eq, crs_ll)
}

M <- make_M(df[,c("longitude","latitude")], buf_km = 200)
M_v <- vect(M)
predM <- mask(crop(pred, M_v), M_v)

# Deduplicate presences per species by raster cell (avoid pseudo-rep)
df$cell <- terra::cellFromXY(predM[[1]], cbind(df$longitude, df$latitude))
df <- df[!is.na(df$cell), ]
df <- dplyr::group_by(df, species) |> dplyr::filter(!duplicated(cell)) |> dplyr::ungroup()
df$cell <- NULL

# Optional: drop highly correlated predictors using background sample in M
set.seed(42)
bg_for_corr <- spatSample(predM, size = 20000, method = "random", na.rm = TRUE, as.points = TRUE)
vals <- extract(predM, bg_for_corr)[, -1, drop = FALSE] #drop ID column
C <- cor(vals, use = "pairwise.complete.obs", method = "pearson")

d  <- as.dist(1 - abs(C))
hc <- hclust(d, method = "average")
groups <- cutree(hc, h = 1 - 0.8) #same as a 0.8 cutoff

# Pick one variable per cluster (by name or by some previous knowledge; here, by name)
keep_names <- unlist(tapply(names(groups), groups, function(v) v[1]))
keep_names #record which bioclim variables were kept!!
## North American kept bioclims
#    1       2       3       4       5       6       7       8 
# "bio1" "bio10" "bio12" "bio15" "bio18"  "bio2"  "bio3"  "bio8"
## South American kept bioclims
#    1       2       3       4       5       6       7       8       9      10 
#"bio1" "bio11" "bio12" "bio14" "bio15" "bio18" "bio19"  "bio2"  "bio3"  "bio4" 

predM <- predM[[keep_names]]

# Sample background points inside M
set.seed(123)
bg_pts <- spatSample(predM, size = 10000, method = "random", na.rm = TRUE, as.points = TRUE)
bg_xy  <- as.data.frame(crds(bg_pts)); names(bg_xy) <- c("longitude","latitude")

head(bg_xy)

# Plot the base map
plot(my_map,
     axes = TRUE, 
     col = "grey95")

# Add the background points
points(bg_xy,
       col = "grey30",
       pch = 1,
       cex = 0.75)

# Add the points for individual observations
points(x = df$longitude, 
       y = df$latitude, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)

##### Run through the model one species at a time (maybe there's a better way to do this)
sp <- "Calyptridium parryi"   #change this each time
#sp <- "Thingia ambigua"
#sp <- "Montiopsis gilliesii"
#sp <- "Philippiamra calycina"
#sp <- "Cistanthe cachinalensis"
#sp <- "Cistanthe longiscapa complex"

pres_xy <- df[df$species == sp, c("longitude","latitude")]

# De-duplicate by raster cell:
cells <- cellFromXY(predM[[1]], as.matrix(pres_xy))
pres_xy <- pres_xy[!duplicated(cells) & !is.na(cells), , drop=FALSE]
nrow(pres_xy)

# Presences
pres_vals <- terra::extract(predM, pres_xy)[, -1, drop = FALSE]
pres_vals <- pres_vals[complete.cases(pres_vals), , drop = FALSE]
pres_dat  <- data.frame(presence = 1, pres_vals)

# Pseudo-absences (up to 10x presences sampled from bg)
set.seed(42)
n_abs   <- min(nrow(bg_xy), 10 * nrow(pres_dat))
abs_xy  <- bg_xy[sample(seq_len(nrow(bg_xy)), n_abs), , drop = FALSE]
abs_vals <- terra::extract(predM, abs_xy)[, -1, drop = FALSE]
abs_vals <- abs_vals[complete.cases(abs_vals), , drop = FALSE]
abs_dat  <- data.frame(presence = 0, abs_vals)

# Assemble dataset and do a 70/30 hold-out
dat <- rbind(pres_dat, abs_dat)
dat <- dat[complete.cases(dat), , drop = FALSE]

set.seed(7)
idx_test <- sample(seq_len(nrow(dat)), size = floor(0.3 * nrow(dat)))
test  <- dat[idx_test, , drop = FALSE]
train <- dat[-idx_test, , drop = FALSE]

# Fit a plain GLM (binomial) with all kept predictors
form <- as.formula(paste("presence ~", paste(colnames(train)[-1], collapse = " + ")))

# Use either line (second one implements bias reduction, good for species with few points)
#fit <- glm(form, data = train, family = binomial())
fit <- glm(form, data = train, family = binomial(), method = brglmFit)  #penalized GLM

# Quick sanity check
summary(fit)$aic    

# Valuate on the test set
p_hat <- predict(fit, newdata = test[test$presence == 1, , drop = FALSE], type = "response")
a_hat <- predict(fit, newdata = test[test$presence == 0, , drop = FALSE], type = "response")
ev    <- dismo::evaluate(p = p_hat, a = a_hat) #you *could* save this for each species under different variable name
ev@auc
thr <- dismo::threshold(ev, stat = "spec_sens")
thr

# Predict map
# I would save the first map with and without points as PDFs
r_pred <- terra::predict(predM, fit, type = "response", na.rm = TRUE)
# Plot the base map
cols <- hcl.colors(50, "YlGnBu")
zlim_fixed <- c(0, 1)
plot(r_pred, main = paste0(sp, " | GLM probability"), zlim = zlim_fixed, plg = FALSE, col = cols)
points(pres_xy, pch = 20, cex = 0.5, col = "red")
fields::image.plot(legend.only=TRUE, zlim= zlim_fixed, col=cols,
                   smallplot=c(0.78,0.80,0.20,0.80))

r_bin  <- r_pred >= thr
plot(r_bin, main = paste0(sp, " | binary (spec_sens)"))

# Save GeoTIFFs too
writeRaster(r_pred, paste0(gsub(" ", "_", sp), "_glm_suitability_NEW.tif"), overwrite = TRUE)
writeRaster(r_bin,  paste0(gsub(" ", "_", sp), "_glm_binary_specsens_NEW.tif"), overwrite = TRUE)

# Save results
out <- list(
  species     = sp,
  keep_names  = keep_names,
  fit         = fit,
  aic         = summary(fit)$aic,
  auc         = ev@auc,
  threshold   = thr,
  coefficients = coef(summary(fit))
)
saveRDS(out, file = paste0(gsub(" ", "_", sp), "_glm_model_NEW.rds"))

tab <- data.frame(
  term = rownames(coef(summary(fit))),
  estimate = coef(summary(fit))[, "Estimate"],
  se = coef(summary(fit))[, "Std. Error"],
  z = coef(summary(fit))[, "z value"],
  p = coef(summary(fit))[, "Pr(>|z|)"]
)
tab$odds_ratio <- exp(tab$estimate)
write.csv(tab, paste0(gsub(" ", "_", sp), "_glm_coefficients_NEW.csv"), row.names = FALSE)

# What are the relative importance of bioclim variables?
train_sc <- train
train_sc[,-1] <- scale(train_sc[,-1])
fit_sc <- glm(form, data = train_sc, family = binomial(), method = brglmFit)
coef(summary(fit_sc))

coef_sc <- data.frame(
  term = rownames(coef(summary(fit_sc))),
  estimate = coef(summary(fit_sc))[, "Estimate"],
  std_error = coef(summary(fit_sc))[, "Std. Error"],
  z = coef(summary(fit_sc))[, "z value"],
  p = coef(summary(fit_sc))[, "Pr(>|z|)"]
)

write.csv(
  coef_sc,
  file = paste0(gsub(" ", "_", sp), "_glm_standardized_coefficients_NEW.csv"),
  row.names = FALSE
)

saveRDS(
  fit_sc,
  file = paste0(gsub(" ", "_", sp), "_glm_standardized_model_NEW.rds")
)

or_sc <- data.frame(
  term = rownames(coef(summary(fit_sc))),
  odds_ratio = exp(coef(fit_sc))
)

write.csv(
  or_sc,
  file = paste0(gsub(" ", "_", sp), "_glm_standardized_odds_ratios_NEW.csv"),
  row.names = FALSE
)

# re-load
GLM_model <- readRDS("Cistanthe_longiscapa_complex_glm_model_NEW.rds")
GLM_model

# Create continuous rasters for each species individually to use for next steps
r_parryi   <- r_pred
#r_gilliesii <- r_pred
#r_ambigua <- r_pred
#r_calycina <- r_pred
#r_cachinalensis <- r_pred
#r_long <- r_pred


############# Stats ############# 
# Example: C. parryi vs T. ambigua
# Use continuous rasters you already created right above
# Run this through every meaningful comparisons you can make

##### Geographic overlap in predicted range
# Choose threshold (can change; right now it's set at 50% predicted probability)
thr <- 0.5      
#or: thr <- dismo::threshold(ev_monandrum, "spec_sens")

# Convert to binary (1 = suitable, NA = unsuitable)
rb_parryi  <- terra::ifel(r_parryi >= thr, 1, NA)
rb_ambigua <- terra::ifel(r_ambigua >= thr, 1, NA)

# sanity check
plot(rb_parryi)
plot(rb_ambigua)

# Cell areas for grid
area_km <- terra::cellSize(predM[[1]], unit = "km")

# Areas of each species (A and B)
A <- terra::global( terra::mask(area_km, rb_parryi), "sum", na.rm = TRUE )[1,1]
B <- terra::global( terra::mask(area_km, rb_ambigua), "sum", na.rm = TRUE )[1,1]

# Intersection mask (1 only where both are 1)
inter_mask <- rb_parryi * rb_ambigua
Iarea <- terra::global( terra::mask(area_km, inter_mask), "sum", na.rm = TRUE )[1,1]

# Union mask (1 where either is 1)
rbA0 <- terra::ifel( is.na(rb_parryi),   0, 1 )
rbB0 <- terra::ifel( is.na(rb_ambigua), 0, 1 )
union_mask <- terra::ifel(rbA0 + rbB0 > 0, 1, NA)
Uarea <- terra::global( terra::mask(area_km, union_mask), "sum", na.rm = TRUE )[1,1]

A; B; Iarea; Uarea

### different ways of calculating overlap
# Jaccard (intersection over union; ∈ [0,1])
jaccard  <- Iarea / Uarea

# Sørensen–Dice (∈ [0,1])
sorensen <- 2 * Iarea / (A + B)

# "Percent overlap" as written by Sophie (∈ [0, 0.5] unless you multiply by 2)
# Maybe not a great measure
friend_pct <- Iarea / (A + B)

# Overlap coefficient (Szymkiewicz–Simpson): intersection over the smaller range (∈ [0,1])
overlap_coeff <- Iarea / min(A, B)

# Save the following output as your result
c(jaccard = jaccard, sorensen = sorensen, friend_pct = friend_pct, overlap_coeff = overlap_coeff)

##### D and I stats using dismo (climatic overlap coefficients)
# Convert terra to raster
rA <- raster::raster(r_parryi)
rB <- raster::raster(r_ambigua)

# Save the following output as your result
dismo::nicheOverlap(rA, rB, stat = "D")
dismo::nicheOverlap(rA, rB, stat = "I")


############ Between species ############ 
# Go through each meaningful comparison
pres_only <- df %>%
  filter(species %in% c("Calyptridium parryi","Thingia ambigua"))

pres_only$species <- factor(
  pres_only$species,
  levels = c("Calyptridium parryi",
             "Thingia ambigua"))

pres_only <- df %>%
  filter(species %in% c("Montiopsis gilliesii","C. longiscapa complex"))

#pres_only <- df %>%
#  filter(species %in% c("Montiopsis gilliesii","Philippiamra calycina"))

#pres_only <- df %>%
#  filter(species %in% c("Montiopsis gilliesii","Cistanthe cachinalensis","Philippiamra calycina"))

pres_only$species <- factor(
  pres_only$species,
  levels = c("Montiopsis gilliesii",
             "C. longiscapa complex"))

#pres_only$species <- factor(
#  pres_only$species,
#  levels = c("Montiopsis gilliesii",
#             "Cistanthe cachinalensis",
#             "Philippiamra calycina"))

X <- terra::extract(predM, pres_only[,c("longitude","latitude")])[, -1, drop=FALSE]
dat_sp <- cbind(species = pres_only$species, X)
dat_sp <- dat_sp[complete.cases(dat_sp), ]

dat_sp$species <- factor(dat_sp$species)
fit_diff <- glm(species ~ ., data = dat_sp, family = binomial())
summary(fit_diff)

pairs(dat_sp[, -1], col = dat_sp$species)

cor(dat_sp[, -1], use = "pairwise.complete.obs")

fit_diff_simple <- glm(
  species ~ bio10,
  data = dat_sp,
  family = binomial()
)
summary(fit_diff_simple)
exp(coef(fit_diff_simple))

# save results
saveRDS(fit_diff,
        file = "species_diff_full_model_THIG_CALYP.rds")

saveRDS(fit_diff_simple,
        file = "species_diff_simple_model_THIG_CALYP_1-10.rds")

sum_full <- coef(summary(fit_diff))

sum_full_df <- data.frame(
  term = rownames(sum_full),
  estimate = sum_full[, "Estimate"],
  std_error = sum_full[, "Std. Error"],
  z = sum_full[, "z value"],
  p = sum_full[, "Pr(>|z|)"]
)

write.csv(
  sum_full_df,
  "species_diff_full_model_coefficients_THIG_CALYP.csv",
  row.names = FALSE
)

sum_simple <- coef(summary(fit_diff_simple))

sum_simple_df <- data.frame(
  term = rownames(sum_simple),
  estimate = sum_simple[, "Estimate"],
  std_error = sum_simple[, "Std. Error"],
  z = sum_simple[, "z value"],
  p = sum_simple[, "Pr(>|z|)"]
)

write.csv(
  sum_simple_df,
  "species_diff_simple_model_coefficients_THIG_CALYP_1-10.csv",
  row.names = FALSE
)

or_all <- data.frame(
  term = names(coef(fit_diff)),
  odds_ratio = exp(coef(fit_diff))
)

write.csv(
  or_all,
  "species_diff_model_odds_ratios_THIG_CALYP.csv",
  row.names = FALSE
)

or_simple <- data.frame(
  term = names(coef(fit_diff_simple)),
  odds_ratio = exp(coef(fit_diff_simple))
)

write.csv(
  or_simple,
  "species_diff_simple_model_odds_ratios_THIG_CALYP_1-12.csv",
  row.names = FALSE
)

# to re-load already-run models
fit_diff_simple <- readRDS("species_diff_simple_model_CALYP_THINGIA.rds")
summary(fit_diff_simple)
fit_diff_simple
contrasts(model.frame(fit_diff_simple)$species)


### Boxplots
# Keep only species + selected bioclim variables
plot_df <- dat_sp[, c("species", keep_names)]

# Convert to long format
plot_long <- tidyr::pivot_longer(
  plot_df,
  cols = all_of(as.character(keep_names)),
  names_to = "bioclim",
  values_to = "value"
)

species_cols <- c(
  "Thingia ambigua"      = "#FF7900",  # muted orange
  "Calyptridium parryi"  = "#2F6B4F"   # dark green
)

species_cols <- c(
  "Montiopsis gilliesii"      = "#D04E4E",  
  "Cistanthe cachinalensis"  = "#3F6FA5",   
  "Philippiamra calycina" = "#DDCC77"
)


species_cols <- c(
  "Montiopsis gilliesii"      = "#D04E4E", 
  "C. longiscapa complex"  = "#805B87"
)

plot_long$bioclim <- factor(plot_long$bioclim, levels = c("bio1","bio2","bio3","bio8","bio10","bio12","bio15","bio18"))

plot_long$bioclim <- factor(plot_long$bioclim, levels = c("bio1","bio2","bio3","bio4","bio11","bio12","bio14","bio15","bio18","bio19"))

library(ggpubr)
ggplot(plot_long, aes(x = species, y = value, fill = species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  facet_wrap(~ bioclim, scales = "free_y") +
  scale_fill_manual(values = species_cols) +
  stat_compare_means(method = "anova",
                     label = "p.signif",
                     label.x.npc = 0.45,
                     label.y.npc = 0.95) +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
  ) +
  labs(
    x = NULL,
    y = "Bioclimatic value",
  ) +
  facet_wrap(~ bioclim, scales = "free_y", ncol = 4)


### more rigorous ANOVA
# one ANOVA per bioclim
pdat <- plot_long %>%
  group_by(bioclim) %>%
  summarise(
    p = summary(aov(value ~ species))[[1]][["Pr(>F)"]][1],
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    stars = case_when(
      p_adj <= 1e-4 ~ "****",
      p_adj <= 1e-3 ~ "***",
      p_adj <= 1e-2 ~ "**",
      p_adj <= 5e-2 ~ "*",
      TRUE ~ "ns"
    ),
    stars = ifelse(stars == "ns", "", stars)
  )

# y-position per facet (works with scales="free_y")
pos <- plot_long %>%
  group_by(bioclim) %>%
  summarise(
    ymin = min(value, na.rm = TRUE),
    ymax = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y = ymax + 0.06 * (ymax - ymin))

# place label at the rightmost species level in each panel
right_species <- levels(plot_long$species)[length(levels(plot_long$species))]

anno <- pdat %>%
  left_join(pos, by = "bioclim") %>%
  mutate(species = right_species)

ggplot(plot_long, aes(x = species, y = value, fill = species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
  facet_wrap(~ bioclim, scales = "free_y") +
  scale_fill_manual(values = species_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  geom_text(
    data = anno,
    aes(x = species, y = y, label = stars),
    inherit.aes = FALSE,
    hjust = 1, vjust = 0,
    size = 3
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
  ) +
  labs(x = NULL, y = "Bioclimatic value")

write.csv(pdat,
          "bioclim_species_ANOVA_BH_results_THIG_CALYP.csv",
          row.names = FALSE)

