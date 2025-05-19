# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

#  Black Eye Syndrome Spatiotemporal Analyses ----------------------------------
#  code by Reyn Yoshioka

# Data Analysis Script -----
# Exploring potential environmental drivers of BES based on ADF&G crab observer
# data and Bering10K model hindcasts.

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# 1: Getting Started -----------------------------------------------------------
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# > 1.1: Clear everything ----- 
# clear workspace
rm(list=ls())
# clear console
cat("\014")

# > 1.2: Load packages ---------------------------------------------------------
library(here)
library(sdmTMB)
library(fmesher)
library(tidyverse)
library(sf) # simple features for spatial handling
library(ggOceanMaps) # Make ocean/bathymetric maps. 
# ^ Note: Use version 2.1.1! Version 2.2.0 has a bug (open issue)
library(ggspatial) # spatial ggplot
library(pROC) # ROC-AUC
library(patchwork) # arrange plots

# > > 1.2.1: here::i_am() -----
# setting working directory
here::i_am("SpatiotemporalBES/BES_Analyses.R")

# > 1.3: Read in data ------
dfm_QO1821 = read.csv(here("BES_dfm_QO1821_20250325_1619.csv"))
df_envdat = read.csv(here("BES_df_envdat_20250325_1620.csv"))

df_statarea = read.csv(here("df_statarea.csv"))

# > 1.4: Minor data prep -----
dfm_QO1821$X = NULL
df_statarea$X = NULL

# remove shells with very few samples (<10)
dfm_QO1821 |>
  filter(shell2 != "new",
         shell2 != "old",
         shell2 != "v.old") |>
  group_by(shell2) |>
  summarize(count = length(shell2))
# 16 crab total in any of those classes

# Percent omitted crab based on shell condition
100 - ((length(dfm_QO1821$shell2) - 16) / length(dfm_QO1821$shell2)) * 100
# So remaining crab
((length(dfm_QO1821$shell2) - 16) / length(dfm_QO1821$shell2)) * 100

dfm_QO1821 = dfm_QO1821[dfm_QO1821$shell2 == "new" |
                          dfm_QO1821$shell2 == "old" |
                          dfm_QO1821$shell2 == "v.old" ,]

dfm_QO1821$shell2 = factor(dfm_QO1821$shell2,
                           levels = c("new",
                                      "old",
                                      "v.old"),
                           ordered = FALSE) # turn to TRUE to make ordinal

dfm_QO1821$shellind = as.numeric(dfm_QO1821$shell2)

# Final dataset size 
length(dfm_QO1821$shell2)

# > 1.5: Standardization and Derived Variables -----
# Means (only for temp)
vec_std_M = dfm_QO1821 |>
  distinct(ind, gridind, .keep_all = TRUE) |> # mean by location-time, not crab
  select(temp_mean24:temp_mean4) |> # only env var columns
  apply(2, mean)

# Standard Deviations (only for temp)
vec_std_SD = dfm_QO1821 |>
  distinct(ind, gridind, .keep_all = TRUE) |> # mean by location-time, not crab
  select(temp_mean24:temp_mean4) |> # only env var columns
  apply(2, sd)

# Calculate
dfm_QO1821 = 
  dfm_QO1821 |>
  mutate(temp_mean4_std = 
           (temp_mean4 - vec_std_M["temp_mean4"]) /
           vec_std_SD["temp_mean4"],
         temp_mean8_std = 
           (temp_mean8 - vec_std_M["temp_mean8"]) /
           vec_std_SD["temp_mean8"],
         temp_mean12_std = 
           (temp_mean12 - vec_std_M["temp_mean12"]) /
           vec_std_SD["temp_mean12"],
         temp_mean16_std = 
           (temp_mean16 - vec_std_M["temp_mean16"]) /
           vec_std_SD["temp_mean16"],
         temp_mean20_std = 
           (temp_mean20 - vec_std_M["temp_mean20"]) /
           vec_std_SD["temp_mean20"],
         temp_mean24_std = 
           (temp_mean24 - vec_std_M["temp_mean24"]) /
           vec_std_SD["temp_mean24"],
         # pH threshold is 7.75 for expected negative effects
         pH_thresh4 = ifelse(pH_mean4 > 7.75, 0, 1),
         pH_thresh8 = ifelse(pH_mean8 > 7.75, 0, 1),
         pH_thresh12 = ifelse(pH_mean12 > 7.75, 0, 1),
         pH_thresh16 = ifelse(pH_mean16 > 7.75, 0, 1),
         pH_thresh20 = ifelse(pH_mean20 > 7.75, 0, 1),
         pH_thresh24 = ifelse(pH_mean24 > 7.75, 0, 1))

df_envdat =
  df_envdat |>
  mutate(pH_thresh4 = ifelse(pH_mean4 > 7.75, 0, 1),
         pH_thresh8 = ifelse(pH_mean8 > 7.75, 0, 1),
         pH_thresh12 = ifelse(pH_mean12 > 7.75, 0, 1),
         pH_thresh16 = ifelse(pH_mean16 > 7.75, 0, 1),
         pH_thresh20 = ifelse(pH_mean20 > 7.75, 0, 1),
         pH_thresh24 = ifelse(pH_mean24 > 7.75, 0, 1))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# 2. Mesh construction ---------------------------------------------------------
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# mesh boundaries
obj_mesh_bounds = 
  fm_extensions(unique(cbind(dfm_QO1821$gridXCkm,
                             dfm_QO1821$gridYCkm)),
                convex = -0.05)

# plug in estimated range to construct mesh
# val_est_range = 100 # For testing
val_est_range = 34.52 # estimated range for final run

# construct mesh
obj_mesh =
  make_mesh(data = dfm_QO1821,
            xy_cols = c("gridXCkm",
                        "gridYCkm"),
            fmesher_func = fmesher::fm_mesh_2d_inla,
            boundary = obj_mesh_bounds,
            max.edge = c(val_est_range / 3, # adjust to 1/10-1/3 spatial range
                         2 * val_est_range / 3), # twice the inner max edge
            offset = c(val_est_range / 3, # about inner max edge
                       val_est_range), # adjust to approximate spatial range
            cutoff = val_est_range / 15)

plot(obj_mesh)

# Extract mesh segments
idx = rbind(obj_mesh$mesh$graph$tv[, 1:2, drop = FALSE], 
            obj_mesh$mesh$graph$tv[, 2:3, drop = FALSE], 
            obj_mesh$mesh$graph$tv[, c(3, 1), drop = FALSE])

df_mesh_segments = data.frame(obj_mesh$mesh$loc[idx[, 1], 1:2],
                              obj_mesh$mesh$loc[idx[, 2], 1:2])
colnames(df_mesh_segments) = c("X", "Y", "Xend", "Yend")

# Convert gridXCkm and gridXYkm back to 3571
# > gridX3571med
# [1] 441358.7
# > gridY3571med
# [1] -3478220

df_mesh_segments = 
  df_mesh_segments |>
  mutate(X = X * 1000 + 441358.7, # Add back median X location of grid
         Y = Y * 1000 - 3478220, # Add back median Y location of grid
         Xend = Xend * 1000 + 441358.7,  # Add back median X location of grid
         Yend = Yend * 1000 - 3478220) # Add back median Y location of grid

# Plot mesh
# Omit if you don't need to plot the map, it takes some time...
# basemap(limits = c(180, -164, 54, 63),
#         lon.interval = 5,
#         lat.interval = 2.5,
#         # crs = 3571,
#         rotate = TRUE,
#         bathy.style = "raster_continuous_blues",
#         # bathy.alpha = 0.5,
#         legends = TRUE,
#         land.col = "black") +
#   geom_spatial_segment(aes(x = X,
#                            y = Y,
#                            xend = Xend,
#                            yend = Yend),
#                        color = "black",
#                        crs = 3571,
#                        alpha = 0.5,
#                        linewidth = 0.1,
#                        data = df_mesh_segments) +
#   annotation_north_arrow(location = "tr",
#                          which_north = "true",
#                          style = north_arrow_orienteering(line_col = "grey75",
#                                                           fill = c("white",
#                                                                    "grey75"),
#                                                           text_col = "grey75",
#                                                           text_size = 6),
#                          height = unit(0.5, "cm"),
#                          width = unit(0.5, "cm")) +
#   annotation_scale(location = "bl",
#                    unit_category = "metric",
#                    style = "ticks",
#                    line_col = "grey25",
#                    text_col = "grey25",
#                    text_cex = 0.5) +
#   theme(text = element_text(family = "serif"))
# 
# ggsave("BES_plot_mesh.png",
#        device = "png",
#        dpi = 300,
#        height = 5,
#        width = 6.5,
#        units = "in")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# NOTE: You can pause here to skip all model (re)fitting and comparison ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# 3. Environmental Window Selection --------------------------------------------
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# Choosing the most predictive window size leading up to crab observation for 
# each assessed environmental variable

# We assume that all demographic variables are informative and maintain them
# in all models.

# Preliminary runs showed AR1 rho = 0.2, so low spatiotemporal auto-
# correlation and we treat the spatiotemporal process as iid. This makes sense
# to some extent, as we only have four fishery seasons to work with.

# > 3.1: Temperature -----
# > > 3.1.1: 4 weeks prior -----
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_temp4.rds"))) {
  mod_temp4 = read_rds(here("BES_mod_temp4.rds"))
} else {
  mod_temp4 =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2 +
             temp_mean4_std,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # save model
  write_rds(mod_temp4, here("BES_mod_temp4.rds"))
  
}

# Record Sanity and start building out cAIC table
mod_temp4_sanity = sanity(mod_temp4)
df_mod_cAIC_temp = data.frame(mod = "temp_mean4_std",
                              condAIC = cAIC(mod_temp4))

# remove model to save space
rm(mod_temp4)

# > > 3.1.2: 8 weeks prior -----
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_temp8.rds"))) {
  mod_temp8 = read_rds(here("BES_mod_temp8.rds"))
} else {
  mod_temp8 =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2 +
             temp_mean8_std,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # save model
  write_rds(mod_temp8, here("BES_mod_temp8.rds"))
  
}

# Record Sanity and start building out cAIC table
mod_temp8_sanity = sanity(mod_temp8)
df_mod_cAIC_temp = 
  df_mod_cAIC_temp |>
  add_row(mod = "temp_mean8_std",
          condAIC = cAIC(mod_temp8))

# remove model to save space
rm(mod_temp8)

# > > 3.1.3: 12 weeks prior -----
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_temp12.rds"))) {
  mod_temp12 = read_rds(here("BES_mod_temp12.rds"))
} else {
  mod_temp12 =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2 +
             temp_mean12_std,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # save model
  write_rds(mod_temp12, here("BES_mod_temp12.rds"))
  
}

# Record Sanity and start building out cAIC table
mod_temp12_sanity = sanity(mod_temp12)
df_mod_cAIC_temp = 
  df_mod_cAIC_temp |>
  add_row(mod = "temp_mean12_std",
          condAIC = cAIC(mod_temp12))

# remove model to save space
rm(mod_temp12)

# > > 3.1.4: 16 weeks prior -----
# Avoid refitting model if it already exists
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_temp16.rds"))) {
  mod_temp16 = read_rds(here("BES_mod_temp16.rds"))
} else {
  mod_temp16 =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2 +
             temp_mean16_std,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # save model
  write_rds(mod_temp16, here("BES_mod_temp16.rds"))
  
}

# Record Sanity and start building out cAIC table
mod_temp16_sanity = sanity(mod_temp16)
df_mod_cAIC_temp = 
  df_mod_cAIC_temp |>
  add_row(mod = "temp_mean16_std",
          condAIC = cAIC(mod_temp16))

# remove model to save space
rm(mod_temp16)

# > > 3.1.5: 20 weeks prior -----
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_temp20.rds"))) {
  mod_temp20 = read_rds(here("BES_mod_temp20.rds"))
} else {
  mod_temp20 =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2 +
             temp_mean20_std,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # save model
  write_rds(mod_temp20, here("BES_mod_temp20.rds"))
  
}

# Record Sanity and start building out cAIC table
mod_temp20_sanity = sanity(mod_temp20)
df_mod_cAIC_temp = 
  df_mod_cAIC_temp |>
  add_row(mod = "temp_mean20_std",
          condAIC = cAIC(mod_temp20))

# remove model to save space
rm(mod_temp20)

# > > 3.1.6: 24 weeks prior -----
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_temp24.rds"))) {
  mod_temp24 = read_rds(here("BES_mod_temp24.rds"))
} else {
  mod_temp24 =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2 +
             temp_mean24_std,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # save model
  write_rds(mod_temp24, here("BES_mod_temp24.rds"))
  
}

# Record Sanity and start building out cAIC table
mod_temp24_sanity = sanity(mod_temp24)
df_mod_cAIC_temp = 
  df_mod_cAIC_temp |>
  add_row(mod = "temp_mean24_std",
          condAIC = cAIC(mod_temp24))

# remove model to save space
rm(mod_temp24)

# > > 3.1.7: Compare temp models -----
df_mod_cAIC_temp = 
  df_mod_cAIC_temp |>
  mutate(dcondAIC = condAIC - min(condAIC))

write_csv(df_mod_cAIC_temp, here("df_mod_cAIC_temp.csv"))

var_temp_selected_window = 
  df_mod_cAIC_temp[which.min(df_mod_cAIC_temp$dcondAIC),]$mod

var_temp_selected_window

# temp 4 weeks preceding

# > 3.2: pH -----
# > > 3.2.1: 4 weeks prior -----
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_pH4.rds"))) {
  mod_pH4 = read_rds(here("BES_mod_pH4.rds"))
} else {
  mod_pH4 =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2 +
             pH_thresh4,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # save model
  write_rds(mod_pH4, here("BES_mod_pH4.rds"))
  
}

# Record Sanity and start building out cAIC table
mod_pH4_sanity = sanity(mod_pH4)
df_mod_cAIC_pH = data.frame(mod = "pH_thresh4",
                            condAIC = cAIC(mod_pH4))

# remove model to save space
rm(mod_pH4)

# > > 3.2.2: 8 weeks prior -----
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_pH8.rds"))) {
  mod_pH8 = read_rds(here("BES_mod_pH8.rds"))
} else {
  mod_pH8 =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2 +
             pH_thresh8,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # save model
  write_rds(mod_pH8, here("BES_mod_pH8.rds"))
  
}

# Record Sanity and start building out cAIC table
mod_pH8_sanity = sanity(mod_pH8)
df_mod_cAIC_pH = 
  df_mod_cAIC_pH |>
  add_row(mod = "pH_thresh8",
          condAIC = cAIC(mod_pH8))

# remove model to save space
rm(mod_pH8)

# > > 3.2.3: 12 weeks prior -----
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_pH12.rds"))) {
  mod_pH12 = read_rds(here("BES_mod_pH12.rds"))
} else {
  mod_pH12 =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2 +
             pH_thresh12,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # save model
  write_rds(mod_pH12, here("BES_mod_pH12.rds"))
  
}

# Record Sanity and start building out cAIC table
mod_pH12_sanity = sanity(mod_pH12)
df_mod_cAIC_pH = 
  df_mod_cAIC_pH |>
  add_row(mod = "pH_thresh12",
          condAIC = cAIC(mod_pH12))

# remove model to save space
rm(mod_pH12)

# > > 3.2.4: 16 weeks prior -----
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_pH16.rds"))) {
  mod_pH16 = read_rds(here("BES_mod_pH16.rds"))
} else {
  mod_pH16 =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2 +
             pH_thresh16,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # save model
  write_rds(mod_pH16, here("BES_mod_pH16.rds"))
  
}

# Record Sanity and start building out cAIC table
mod_pH16_sanity = sanity(mod_pH16)
df_mod_cAIC_pH = 
  df_mod_cAIC_pH |>
  add_row(mod = "pH_thresh16",
          condAIC = cAIC(mod_pH16))

# remove model to save space
rm(mod_pH16)

# > > 3.2.5: 20 weeks prior -----
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_pH20.rds"))) {
  mod_pH20 = read_rds(here("BES_mod_pH20.rds"))
} else {
  mod_pH20 =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2 +
             pH_thresh20,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # save model
  write_rds(mod_pH20, here("BES_mod_pH20.rds"))
  
}

# Record Sanity and start building out cAIC table
mod_pH20_sanity = sanity(mod_pH20)
df_mod_cAIC_pH = 
  df_mod_cAIC_pH |>
  add_row(mod = "pH_thresh20",
          condAIC = cAIC(mod_pH20))

# remove model to save space
rm(mod_pH20)

# > > 3.2.6: 24 weeks prior -----
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_pH24.rds"))) {
  mod_pH24 = read_rds(here("BES_mod_pH24.rds"))
} else {
  mod_pH24 =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2 +
             pH_thresh24,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # save model
  write_rds(mod_pH24, here("BES_mod_pH24.rds"))
  
}

# Record Sanity and start building out cAIC table
mod_pH24_sanity = sanity(mod_pH24)
df_mod_cAIC_pH = 
  df_mod_cAIC_pH |>
  add_row(mod = "pH_thresh24",
          condAIC = cAIC(mod_pH24))

# remove model to save space
rm(mod_pH24)

# > > 3.2.7: Compare pH models -----
df_mod_cAIC_pH = 
  df_mod_cAIC_pH |>
  mutate(dcondAIC = condAIC - min(condAIC))

write_csv(df_mod_cAIC_pH, here("df_mod_cAIC_pH.csv"))

var_pH_selected_window = 
  df_mod_cAIC_pH[which.min(df_mod_cAIC_pH$dcondAIC),]$mod

var_pH_selected_window

# pH 8 weeks preceding

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# 4: Environmental Variable Selection ------------------------------------------
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# > 4.1: Pull cAICs from environmental window selection ----
# No need to rerun models of single environmental variables
df_mod_cAIC_env = 
  rbind(df_mod_cAIC_temp[df_mod_cAIC_temp$dcondAIC == 0,],
        df_mod_cAIC_pH[df_mod_cAIC_pH$dcondAIC == 0,])

# set to NA, as these are not necessarily the best models for this set
df_mod_cAIC_env$dcondAIC = NA 

# > 4.2: Demographic Model -----
# Base model, demographic variables only
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_demo.rds"))) {
  mod_demo = read_rds(here("BES_mod_demo.rds"))
} else {
  mod_demo =
    sdmTMB(BES ~
             CPUE_std +
             CPUE_ind_std +
             size_std +
             shell2,
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # Save model
  write_rds(mod_demo, here("BES_mod_demo.rds"))
  
}

# Update cAIC table
mod_demo_sanity = sanity(mod_demo)
df_mod_cAIC_env = 
  df_mod_cAIC_env |>
  add_row(mod = "demo",
          condAIC = cAIC(mod_demo))

# Remove model for space
rm(mod_demo) 

# > 4.3: Combined Model -----
# Avoid refitting model if it already exists
if(file.exists(here("BES_mod_temp_pH.rds"))) {
  mod_temp_pH = read_rds(here("BES_mod_temp_pH.rds"))
} else {
  mod_temp_pH =
    sdmTMB(formula = 
             as.formula(paste("BES ~ CPUE_std + CPUE_ind_std + size_std + shell2",
                              var_temp_selected_window,
                              var_pH_selected_window,
                              sep = " + ")),
           data = dfm_QO1821,
           mesh = obj_mesh,
           family = binomial(link = "logit"),
           time = "fishery_ind",
           spatiotemporal = "iid")
  
  # Save model
  write_rds(mod_temp_pH, here("BES_mod_temp_pH.rds"))
  
}

# Update cAIC table
mod_temp_pH_sanity = sanity(mod_temp_pH)
df_mod_cAIC_env = 
  df_mod_cAIC_env |>
  add_row(mod = "temp_pH",
          condAIC = cAIC(mod_temp_pH))

# Remove model to save space
rm(mod_temp_pH) 

# > 4.4: Compare Environmental Variable models -----
df_mod_cAIC_env = 
  df_mod_cAIC_env |>
  mutate(dcondAIC = condAIC - min(condAIC))

write_csv(df_mod_cAIC_env, here("df_mod_cAIC_env.csv"))

var_env_selected = 
  df_mod_cAIC_env[which.min(df_mod_cAIC_env$dcondAIC),]$mod

var_env_selected

# full model selected

# > 4.5: Combine all cAIC tables -----
df_mod_cAIC_all = 
  rbind(
    df_mod_cAIC_temp,
    df_mod_cAIC_pH,
    df_mod_cAIC_env
  ) |>
  mutate(dcondAIC = condAIC - min(condAIC)) |>
  unique()

write.csv(df_mod_cAIC_all, here("df_mod_cAIC_all.csv"))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# NOTE: You can resume here after skipping model (re)fitting and comparison ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# 5: Final model ---------------------------------------------------------------
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# > 5.1: Load final model -----
# Full model was selected

# (Un)comment below to toggle declaring the selected model as the full model
# (This is if model (re)fitting and comparison is skipped)
var_env_selected = "temp_pH"

# Load final model
mod_fin = read_rds(here(paste0("BES_mod_",
                               var_env_selected,
                               ".rds")))

# check out model fit
summary(mod_fin)

# double check sanity 
sanity(mod_fin)

# > 5.2: Model Checking and Diagnostics -----
# > > 5.2.1: Residual Diagnostics -----
# QQ Plot
dfm_QO1821$resid = residuals(mod_fin)
qqnorm(dfm_QO1821$resid)
qqline(dfm_QO1821$resid)

# Spatiotemporal mean residuals
ggplot() +
  geom_point(aes(x = gridXCkm,
                 y = gridYCkm,
                 color = resid_M),
             size = 3,
             data = 
               dfm_QO1821 |>
               group_by(gridXCkm,
                        gridYCkm,
                        fishery_ind) |>
               summarize(resid_M = mean(resid))) +
  scale_color_gradient2(high = "red",
                        low = "blue",
                        mid = "grey50") +
  coord_fixed() +
  labs(x = "X (km)",
       y = "Y (km)",
       color = "mean residual") +
  theme_bw() +
  facet_grid(. ~ fishery_ind)

# > > 5.2.2: ROC-AUC Evaluation -----
dfm_QO1821$est = predict(mod_fin)$est

roc(BES ~ est,
    dfm_QO1821,
    auc = TRUE,
    plot = TRUE)

# > 5.3: Visualizing Estimates of Final Model ----- 
# > > 5.3.1: Extract fixed effects -----
df_mod_fin_fixed = tidy(mod_fin, effects = "fixed")
df_mod_fin_fixed$term = factor(df_mod_fin_fixed$term,
                               levels = c("(Intercept)",
                                          "CPUE_std",
                                          "CPUE_ind_std",
                                          "size_std",
                                          "shell2old",
                                          "shell2v.old",
                                          "temp_mean4_std",
                                          "pH_thresh8"))

# exponentiate estimates for odds ratios
df_mod_fin_fixed = 
  df_mod_fin_fixed |>
  mutate(estimate_exp = exp(estimate),
         std.error_exp_upr = exp(estimate + std.error),
         std.error_exp_low = exp(estimate - std.error),
         CI95_exp_upr = exp(estimate + 1.96 * std.error),
         CI95_exp_low = exp(estimate - 1.96 * std.error))

# colors 
vec_BES_colors =
  c(
    "(Intercept)" = "grey20",
    "CPUE_std" = "black",
    "CPUE_ind_std" = "black",
    "size_std" = "black",
    "new" = "orange",
    "old" = "orange3",
    "v.old" = "orange4",
    "shell2old" = "orange3",
    "shell2v.old" = "orange4",
    "temp_mean4_std" = "firebrick",
    "pH_thresh8" = "magenta3"
  )

# create lookup for labeling terms
df_term_labs = 
  data.frame(term = c("(Intercept)",
                      "CPUE_std",
                      "CPUE_ind_std",
                      "size_std",
                      "shell2old",
                      "shell2v.old",
                      "temp_mean4_std",
                      "pH_thresh8"),
             term_lab = c("intercept",
                          "fishery CPUE, std.",
                          "local CPUE, std.",
                          "size, std.",
                          "shell: old",
                          "shell: very old",
                          "temp., std.",
                          "pH, thresh."))

# match look up
df_mod_fin_fixed = 
  df_mod_fin_fixed |>
  # for easier editing, omit term_lab if it already exists
  select(!any_of("term_lab")) |> 
  left_join(df_term_labs)

# order term labels
df_mod_fin_fixed$term_lab = 
  factor(df_mod_fin_fixed$term_lab,
         levels = c("intercept",
                    "fishery CPUE, std.",
                    "local CPUE, std.",
                    "size, std.",
                    "shell: old",
                    "shell: very old",
                    "temp., std.",
                    "pH, thresh."))

# graph odds ratios
ggplot(data = df_mod_fin_fixed) +
  geom_errorbar(aes(x = term_lab,
                    ymin = CI95_exp_low,
                    ymax = CI95_exp_upr,
                    color = term),
                width = 0.0) +
  geom_point(aes(x = term_lab,
                 y = estimate_exp,
                 color = term)) +
  geom_hline(yintercept = 1,
             linetype = 2) +
  scale_color_manual(values = vec_BES_colors,
                     guide = "none") +
  labs(y = "odds ratio\nexp(β)",
       x = "model term") +
  scale_x_discrete(limits = rev) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(family = "serif"))

# export graph
# ggsave("BES_plot_fixed.png",
#        device = "png",
#        dpi = 300,
#        height = 2,
#        width = 3,
#        units = "in")

# > > 5.3.2: Visualize conditional effects -----
# > > > 5.3.2.1: Temperature effect -----
df_pred_temp = 
  predict(mod_fin,
          newdata = data.frame(gridXCkm = rep(0, 20),
                               gridYCkm = rep(0, 20),
                               fishery_ind = 1,
                               temp_mean4_std = seq(min(dfm_QO1821$temp_mean4_std),
                                                    max(dfm_QO1821$temp_mean4_std),
                                                    length.out = 20),
                               pH_thresh8 = rep(0, 20),
                               O2_thresh16 = rep(0, 20),
                               CPUE_std = rep(0, 20),
                               CPUE_ind_std = rep(0, 20),
                               size_std = rep(0, 20),
                               shell2 = "old"),
          re_form = ~ 0,
          se_fit = TRUE)

val_temp_increment = 
  mean(diff(seq(min(dfm_QO1821$temp_mean4),
                max(dfm_QO1821$temp_mean4),
                length.out = 20))) / 2

df_N_temp = 
  dfm_QO1821 |>
  mutate(bin = cut(temp_mean4,
                   breaks = seq(min(temp_mean4) - 0.001, # tiny pad
                                max(temp_mean4),
                                length.out = 20),
                   labels = seq(min(temp_mean4) + val_temp_increment - 0.001,
                                max(temp_mean4) + val_temp_increment,
                                length.out = 19)),
         include.lowest	= TRUE) |>
  group_by(BES, bin) |>
  summarise(N = n()) |>
  mutate(prop = ifelse(BES == 0,
                       N / sum(N) / 2,
                       1 - N / sum(N) / 2),
         bin = as.numeric(as.character(bin)))

plot_temp = 
  ggplot() +
  geom_ribbon(aes(x = temp_mean4_std * 
                    vec_std_SD["temp_mean4"] +
                    vec_std_M["temp_mean4"],
                  ymax = plogis(est + 1.96 * est_se),
                  ymin = plogis(est - 1.96 * est_se)),
              fill = vec_BES_colors["temp_mean4_std"],
              alpha = 0.2,
              data = df_pred_temp) +
  geom_segment(aes(x = bin,
                   xend = bin,
                   y = BES * 0.08,
                   yend = prop * 0.08),
               linewidth = 3,
               alpha = 0.2,
               color = vec_BES_colors["temp_mean4_std"],
               data = df_N_temp) +
  geom_line(aes(x = temp_mean4_std * 
                  vec_std_SD["temp_mean4"] +
                  vec_std_M["temp_mean4"],
                y = plogis(est)),
            linewidth = 1,
            color = vec_BES_colors["temp_mean4_std"],
            data = df_pred_temp) +
  # geom_hline(yintercept = 0.08,
  #            linetype = 2,
  #            color = "grey50") +
  coord_cartesian(ylim = c(0, 0.08)) +
  labs(x = "temperature (°C)",
       y = "BES conditional\nprobability") +
  theme_bw() +
  theme(text = element_text(family = "serif"))
plot_temp

# > > > 5.3.2.2: pH effect -----
df_pred_pH = 
  predict(mod_fin,
          newdata = data.frame(gridXCkm = rep(0, 2),
                               gridYCkm = rep(0, 2),
                               fishery_ind = 1,
                               temp_mean4_std = rep(0, 2),
                               pH_thresh8 = c(0, 1),
                               O2_thresh16 = rep(0, 2),
                               CPUE_std = rep(0, 2),
                               CPUE_ind_std = rep(0, 2),
                               size_std = rep(0, 2),
                               shell2 = "old"),
          re_form = ~ 0,
          se_fit = TRUE)

df_pred_pH = 
  dfm_QO1821 |>
  group_by(pH_thresh8) |>
  summarize(N = length(pH_thresh8)) |>
  right_join(df_pred_pH)

df_pred_pH$pH = factor(c("normal",
                         "acidic"),
                       levels = c("normal",
                                  "acidic"))

plot_pH = 
  ggplot() +
  geom_errorbar(aes(x = pH,
                    ymax = plogis(est + 1.96 * est_se),
                    ymin = plogis(est - 1.96 * est_se)),
                width = 0.0,
                linewidth = 1,
                color = vec_BES_colors["pH_thresh8"],
                data = df_pred_pH) +
  geom_point(aes(x = pH,
                 y = plogis(est)),
             size = 3,
             alpha = 0.5,
             color = vec_BES_colors["pH_thresh8"],
             data = df_pred_pH) +
  geom_text(aes(x = pH,
                y = plogis(est + 1.96 * est_se) + 0.005,
                label = format(N, big.mark=",")),
            color = vec_BES_colors["pH_thresh8"],
            size = unit(3, "pt"),
            family = "serif",
            data = df_pred_pH)+
  # geom_hline(yintercept = 0.08,
  #            linetype = 2,
  #            color = "grey50") +
  coord_cartesian(ylim = c(0, 0.08)) +
  labs(x = "pH",
       y = "BES conditional\nprobability") +
  theme_bw() +
  theme(text = element_text(family = "serif"))
plot_pH

# > > > 5.3.2.3: Size effect -----

# Predict/plot all the way to max size?
# val_size_max = max(dfm_QO1821$size_std) # Yes
val_size_max = 150 # No, constrain

df_pred_size = 
  predict(mod_fin,
          newdata = data.frame(gridXCkm = rep(0, 20),
                               gridYCkm = rep(0, 20),
                               fishery_ind = 1,
                               temp_mean4_std = rep(0, 20),
                               pH_thresh8 = rep(0, 20),
                               O2_thresh16 = rep(0, 20),
                               CPUE_std = rep(0, 20),
                               CPUE_ind_std = rep(0, 20),
                               size_std = seq(min(dfm_QO1821$size_std),
                                              (val_size_max - 98.74568) / 
                                                10.67769,
                                              length.out = 20),
                               shell2 = "old"),
          re_form = ~ 0,
          se_fit = TRUE)

# For reference:
#   > SD_size
# [1] 10.67769
# > M_size
# [1] 98.74568

val_size_increment = 
  mean(diff(seq(min(dfm_QO1821$size),
                val_size_max,
                length.out = 20))) / 2

df_N_size = 
  dfm_QO1821 |>
  mutate(bin = cut(size,
                   breaks = seq(min(size) - 0.001, # tiny pad
                                val_size_max,
                                length.out = 20),
                   labels = seq(min(size) + val_size_increment - 0.001,
                                val_size_max + val_size_increment,
                                length.out = 19)),
         include.lowest	= TRUE) |>
  group_by(BES, bin) |>
  summarise(N = n()) |>
  mutate(prop = ifelse(BES == 0,
                       N / sum(N) / 2,
                       1 - N / sum(N) / 2),
         bin = as.numeric(as.character(bin)))
# NAs here indicate crab outside of the selected range (larger crab)

plot_size = 
  ggplot() +
  geom_ribbon(aes(x = size_std * 
                    10.67769 +
                    98.74568,
                  ymax = plogis(est + 1.96 * est_se),
                  ymin = plogis(est - 1.96 * est_se)),
              alpha = 0.2,
              data = df_pred_size) +
  geom_segment(aes(x = bin,
                   xend = bin,
                   y = BES * 0.08,
                   yend = prop * 0.08),
               linewidth = 2,
               alpha = 0.2,
               data = df_N_size) +
  geom_line(aes(x = size_std * 
                  10.67769 +
                  98.74568,
                y = plogis(est)),
            linewidth = 1,
            data = df_pred_size) +
  # geom_hline(yintercept = 0.08,
  #            linetype = 2,
  #            color = "grey50") +
  coord_cartesian(ylim = c(0, 0.08),
                  xlim = c(50, val_size_max + val_size_increment)) + # match environmental
  # coord_cartesian(ylim = c(0, 0.32)) + # match crab size
  labs(x = "carapace width (mm)",
       y = "BES conditional\nprobability") +
  theme_bw() +
  theme(text = element_text(family = "serif"))
plot_size

# > > > 5.3.2.4: Local CPUE effect -----
df_pred_CPUE_ind = 
  predict(mod_fin,
          newdata = data.frame(gridXCkm = rep(0, 20),
                               gridYCkm = rep(0, 20),
                               fishery_ind = 1,
                               temp_mean4_std = rep(0, 20),
                               pH_thresh8 = rep(0, 20),
                               O2_thresh16 = rep(0, 20),
                               CPUE_std = rep(0, 20),
                               CPUE_ind_std = seq(min(dfm_QO1821$CPUE_ind_std),
                                                  max(dfm_QO1821$CPUE_ind_std),
                                                  length.out = 20),
                               size_std = rep(0, 20),
                               shell2 = "old"),
          re_form = ~ 0,
          se_fit = TRUE)

# Note: CPUE values for standardization:
# Fishery-level CPUE:
# > M_CPUE
# [1] 184.4186
# > SD_CPUE
# [1] 65.13364

# Local CPUE:
# > CPUE_ind_M
# [1] 339.7068
# > CPUE_ind_SD
# [1] 252.2821

val_CPUE_ind_increment = 
  mean(diff(seq(min(dfm_QO1821$CPUE_ind_std),
                max(dfm_QO1821$CPUE_ind_std),
                length.out = 20))) / 2

df_N_CPUE_ind = 
  dfm_QO1821 |>
  mutate(bin = cut(CPUE_ind_std,
                   breaks = seq(min(CPUE_ind_std)  - 0.001, # tiny pad
                                max(CPUE_ind_std),
                                length.out = 20),
                   labels = seq(min(CPUE_ind_std) + val_CPUE_ind_increment - 
                                  0.001,
                                max(CPUE_ind_std) + val_CPUE_ind_increment,
                                length.out = 19)),
         include.lowest	= TRUE) |>
  group_by(BES, bin) |>
  summarise(N = n()) |>
  mutate(prop = ifelse(BES == 0,
                       N / sum(N) / 2,
                       1 - N / sum(N) / 2),
         bin = as.numeric(as.character(bin)))

plot_CPUE_ind = 
  ggplot() +
  geom_ribbon(aes(x = CPUE_ind_std,
                  ymax = plogis(est + 1.96 * est_se),
                  ymin = plogis(est - 1.96 * est_se)),
              alpha = 0.2,
              data = df_pred_CPUE_ind) +
  geom_segment(aes(x = bin,
                   xend = bin,
                   y = BES * 0.08,
                   yend = prop * 0.08),
               linewidth = 2,
               alpha = 0.2,
               data = df_N_CPUE_ind) +
  geom_line(aes(x = CPUE_ind_std,
                y = plogis(est)),
            linewidth = 1,
            data = df_pred_CPUE_ind) +
  # geom_hline(yintercept = 0.08,
  #            linetype = 2,
  #            color = "grey50") +
  coord_cartesian(ylim = c(0, 0.08)) + # match environmental
  # coord_cartesian(ylim = c(0, 0.32)) + # match crab size
  labs(x = "local crab density, std.",
       y = "BES conditional\nprobability") +
  theme_bw() +
  theme(text = element_text(family = "serif"))
plot_CPUE_ind

# > > > 5.3.2.5: Shell Condition effect -----
df_pred_shell = 
  predict(mod_fin,
          newdata = data.frame(gridXCkm = rep(0, 3),
                               gridYCkm = rep(0, 3),
                               fishery_ind = 1,
                               temp_mean4_std = rep(0, 3),
                               pH_thresh8 = rep(0, 3),
                               O2_thresh16 = rep(0, 3),
                               CPUE_std = rep(0, 3),
                               CPUE_ind_std = rep(0, 3),
                               size_std = rep(0, 3),
                               shell2 = c("new",
                                          "old",
                                          "v.old")),
          re_form = ~ 0,
          se_fit = TRUE)

df_pred_shell = 
  dfm_QO1821 |>
  group_by(shell2) |>
  summarize(N = length(shell2)) |>
  right_join(df_pred_shell)

df_pred_shell$shell = factor(c("new",
                               "old",
                               "very old"),
                             levels = c("new",
                                        "old",
                                        "very old"))

plot_shell = 
  ggplot() +
  geom_errorbar(aes(x = shell,
                    ymax = plogis(est + 1.96 * est_se),
                    ymin = plogis(est - 1.96 * est_se),
                    color = shell2),
                width = 0.0,
                linewidth = 1,
                data = df_pred_shell) +
  geom_point(aes(x = shell,
                 y = plogis(est),
                 color = shell2),
             size = 3,
             alpha = 0.5,
             data = df_pred_shell) +
  geom_text(aes(x = shell,
                y = plogis(est + 1.96 * est_se) + 0.005,
                label = format(N, big.mark=","),
                color = shell2),
            size = unit(3, "pt"),
            family = "serif",
            data = df_pred_shell)+
  # geom_hline(yintercept = 0.08,
  #            linetype = 2,
  #            color = "grey50") +
  scale_color_manual(values = vec_BES_colors,
                     guide = "none") +
  coord_cartesian(ylim = c(0, 0.08)) + # match environmental
  # coord_cartesian(ylim = c(0, 0.32)) + # match crab size
  labs(x = "shell condition",
       y = "BES conditional\nprobability") +
  theme_bw() +
  theme(text = element_text(family = "serif"))
plot_shell

# > > > 5.3.2.6: Combine Figures -----
# Environmental
(plot_temp | plot_pH)  +
  plot_layout(axes = "collect_y",
              widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "serif"))

ggsave("BES_plot_cond_env.png",
       device = "png",
       dpi = 300,
       height = 3,
       width = 6,
       units = "in")

# Demographic
(plot_size | plot_shell | plot_CPUE_ind) +
  plot_layout(axes = "collect_y",
              widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "serif"))

ggsave("BES_plot_cond_demo.png",
       device = "png",
       dpi = 300,
       height = 3,
       width = 6,
       units = "in")

# All
(plot_size + plot_shell + plot_CPUE_ind + plot_layout(axes = "collect_y")) /
  (plot_temp + plot_pH + plot_layout(axes = "collect_y")) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "serif"),
        plot.tag.location = "panel",
        plot.tag.position = c(0.1, 0.90))

ggsave("BES_plot_cond_all.png",
       device = "png",
       dpi = 300,
       height = 6,
       width = 6,
       units = "in")

# > > > 5.3.2.7: Generate effect text -----
# Size
# Minimum legal vs market preferred
# 3.1 in = ~78 mm
# 4 in = ~101 mm
df_pred_misc = 
  predict(mod_fin,
          newdata = data.frame(gridXCkm = 0,
                               gridYCkm = 0,
                               fishery_ind = 1,
                               temp_mean4_std = 0,
                               pH_thresh8 = 0,
                               O2_thresh16 = 0,
                               CPUE_std = 0,
                               CPUE_ind_std = 0,
                               size_std = (c(78, 101) - 98.74568) / # standardize
                                 10.67769,
                               shell2 = "old"),
          re_form = ~ 0) |>
  mutate(prob = plogis(est))
# increase from 78 to 101 mm
(df_pred_misc$prob[2] / df_pred_misc$prob[1] - 1) * 100 

# Shell Condition
df_pred_misc = 
  predict(mod_fin,
          newdata = data.frame(gridXCkm = 0,
                               gridYCkm = 0,
                               fishery_ind = 1,
                               temp_mean4_std =
                                 (3 - vec_std_M["temp_mean4"]) / # standardize
                                 vec_std_SD["temp_mean4"],
                               pH_thresh8 = 0,
                               O2_thresh16 = 0,
                               CPUE_std = 0,
                               CPUE_ind_std = 0,
                               size_std = 0,
                               shell2 = c("new",
                                          "old",
                                          "v.old")),
          re_form = ~ 0) |>
  mutate(prob = plogis(est))
# new vs old shell
(df_pred_misc$prob[2] / df_pred_misc$prob[1] - 1) * 100 
# new vs very old shell
(df_pred_misc$prob[3] / df_pred_misc$prob[1] - 1) * 100 

# Temperature:
df_pred_misc = 
  predict(mod_fin,
          newdata = data.frame(gridXCkm = 0,
                               gridYCkm = 0,
                               fishery_ind = 1,
                               temp_mean4_std =
                                 (c(2, 3, 4) - vec_std_M["temp_mean4"]) / # standardize
                                 vec_std_SD["temp_mean4"],
                               pH_thresh8 = 0,
                               O2_thresh16 = 0,
                               CPUE_std = 0,
                               CPUE_ind_std = 0,
                               size_std = (101 - 98.74568) / # standardize
                                 10.67769,
                               shell2 = "old"),
          re_form = ~ 0) |>
  mutate(prob = plogis(est))
# Increase temp 2 to 3
(df_pred_misc$prob[2] / df_pred_misc$prob[1] - 1) * 100 
# Increase temp 3 to 4
(df_pred_misc$prob[3] / df_pred_misc$prob[2] - 1) * 100 

# pH:
df_pred_misc = 
  predict(mod_fin,
          newdata = data.frame(gridXCkm = 0,
                               gridYCkm = 0,
                               fishery_ind = 1,
                               temp_mean4_std = 0,
                               pH_thresh8 = c(0, 1),
                               O2_thresh16 = 0,
                               CPUE_std = 0,
                               CPUE_ind_std = 0,
                               size_std = (101 - 98.74568) / # standardize
                                 10.67769,
                               shell2 = "old"),
          re_form = ~ 0) |>
  mutate(prob = plogis(est))
# Normal pH to acidic
(df_pred_misc$prob[2] / df_pred_misc$prob[1] - 1) * 100 

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# 6: Maps ----------------------------------------------------------------------
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# Labels for plotting
fishery_labels = c(
  QO18 = "2018-2019",
  QO19 = "2019-2020",
  QO20 = "2020-2021",
  QO21 = "2021-2022"
)

# Aggregate observations by fishery season
dfm_QO1821 |>
  group_by(fishery) |>
  summarise(prev = mean(BES, na.rm = TRUE),
            N = length(BES),
            SE = sqrt(prev * (1 - prev) / N))

# Aggregate observations to statistical area
dfm_QO1821_prev = dfm_QO1821 |>
  group_by(statarea, fishery) |>
  summarise(prev = mean(BES, na.rm = TRUE),
            N = length(BES))
dfm_QO1821_prev = left_join(dfm_QO1821_prev,
                            df_statarea)

# Aggregate to the grid cell
dfm_QO1821_prev_gridind =
  dfm_QO1821 |>
  group_by(gridind,
           fishery,
           gridX3571,
           gridY3571) |>
  summarize(prev = mean(BES),
            N = length(BES))

# How many grid cell-season combos had nonzero prevalence
length(dfm_QO1821_prev_gridind[dfm_QO1821_prev_gridind$prev > 0,]$N) /
  length(dfm_QO1821_prev_gridind$N)

# Not a map, but helpful...
plot_hist_prev = 
  ggplot() +
  geom_histogram(aes(x = prev * 100),
                 fill = "firebrick4",
                 binwidth = 2.5,
                 data = dfm_QO1821_prev_gridind) +
  labs(x = "BES Prevalence (%)",
       y = "# grid cells") +
  scale_y_continuous(expand = expansion(mult = 0,
                                        add = c(1, 1))) +
  coord_cartesian(xlim = c(0, 40)) +
  theme_classic() +
  theme(text = element_text(family = "serif"),
        strip.text = element_text(size = unit (8, "pt")),
        axis.title = element_text(size = unit(8, "pt")),
        axis.text = element_text(size = unit(6, "pt")),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5)) +
  facet_grid(. ~ fishery, ,
             labeller = labeller(fishery = fishery_labels))

plot_hist_prev

# export graph
# ggsave("BES_plot_hist.png",
#        device = "png",
#        dpi = 300,
#        height = 2,
#        width = 3,
#        units = "in")

# > 6.1: Observed Prevalence Map -----

# What proportion of stations are plotted?
length(dfm_QO1821_prev[dfm_QO1821_prev$plot,]$statarea) /
  length(dfm_QO1821_prev$statarea)

# Plot with gridind - this does not get released for confidentiality
ggplot() +
  geom_point(aes(x = gridX3571,
                 y = gridY3571,
                 color = prev),
             size = 2,
             data = dfm_QO1821_prev_gridind) +
  scale_color_viridis_c(option = "inferno",
                        direction = -1) +
  coord_fixed() +
  labs(x = "X (km)",
       y = "Y (km)",
       color = "mean prevalence") +
  theme_bw() +
  facet_grid(. ~ fishery)

# Plot for manuscript: 
plot_map_prev = 
  basemap(limits = c(177, -162, 53, 63),
          lon.interval = 5,
          lat.interval = 2.5,
          # crs = 3571,
          rotate = TRUE,
          bathy.style = "raster_continuous_blues",
          # bathy.alpha = 0.5,
          legends = TRUE,
          land.col = "black") +
  geom_spatial_point(aes(x = saX3571,
                         y = saY3571),
                     shape = 19,
                     color = "black",
                     size = 2,
                     alpha = 1,
                     crs = 3571,
                     data = dfm_QO1821_prev[sample(nrow(dfm_QO1821_prev)) &
                                              dfm_QO1821_prev$plot,]) +
  geom_spatial_point(aes(x = saX3571,
                         y = saY3571,
                         # size = prev,
                         color = prev),
                     # shape = 1,
                     # shape = 13,
                     shape = 19,
                     size = 1,
                     alpha = 1,
                     crs = 3571,
                     data = dfm_QO1821_prev[sample(nrow(dfm_QO1821_prev)) &
                                              dfm_QO1821_prev$plot,]) +
  geom_spatial_text(aes(x = 177,
                        y = 61.5,
                        label = "A"),
                    color = "black",
                    fontface = "bold",
                    family = "serif",
                    crs = 4326,
                    data = data.frame(fishery = "QO18")) +
  annotation_north_arrow(location = "tr",
                         which_north = "true",
                         style = north_arrow_orienteering(line_col = "grey75",
                                                          fill = c("white",
                                                                   "grey75"),
                                                          text_col = "grey75",
                                                          text_size = 6),
                         height = unit(0.5, "cm"),
                         width = unit(0.5, "cm"),
                         data = data.frame(fishery = "QO18")) +
  annotation_scale(location = "bl",
                   unit_category = "metric",
                   style = "ticks",
                   line_col = "grey25",
                   text_col = "grey25",
                   text_cex = 0.5,
                   data = data.frame(fishery = "QO18")) +
  scale_color_viridis_c(option = "inferno",
                        guide = "legend",
                        direction = -1,
                        begin = 0.0,
                        end = 0.95,
                        limits = c(0, 0.155)) +
  scale_size(range = c(1, 5)) +
  labs(color = "BES\nprevalence") +
  # coord_sf(crs = st_crs(3571),
  #          xlim = c(min(df_mesh_points$gridX3571),
  #                   max(df_mesh_points$gridX3571)),
  #          ylim = c(min(df_mesh_points$gridY3571),
  #                   max(df_mesh_points$gridY3571))) +
  theme(text = element_text(family = "serif"),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(size = unit(6, "pt")),
        axis.text.y = element_text(size = unit(6, "pt")),
        axis.title = element_blank(),
        legend.title = element_text(size = unit(8, "pt")),
        legend.key.size = unit(6,"pt"),
        legend.text = element_text(size = unit(6, "pt"))) +
  guides(color = "colorbar",
         fill = "none") +
  # ggtitle("A) Observed BES Prevalence") +
  facet_grid(. ~ fishery,
             labeller = labeller(fishery = fishery_labels))

plot_map_prev

# plot_map_prev + plot_hist_prev +
#   plot_layout(nrow = 2)
# 
# # export graph
# ggsave("BES_plot_map_hist.png",
#        device = "png",
#        dpi = 600,
#        height = 3.5,
#        width = 6.5,
#        units = "in")


# > 6.2: Temp Map -----
plot_map_temp = basemap(limits = c(177, -162, 53, 63),
                        lon.interval = 5,
                        lat.interval = 2.5,
                        rotate = TRUE,
                        bathymetry = FALSE,
                        # bathy.style = "poly_blues",
                        # bathy.alpha = 0.5,
                        legends = FALSE,
                        land.col = "black") +
  geom_spatial_point(aes(x = gridX3571,
                         y = gridY3571,
                         color = temp_mean4),
                     shape = 16,
                     size = 0.1,
                     crs = 3571,
                     data = df_envdat) +
  geom_spatial_text(aes(x = 177,
                        y = 61.5,
                        label = "B"),
                    color = "black",
                    fontface = "bold",
                    family = "serif",
                    crs = 4326,
                    data = data.frame(fishery = "QO18")) +
  scale_color_gradient2(high = "firebrick",
                        low = "dodgerblue3",
                        mid = "grey90",
                        midpoint = 2) +
  labs(color = "Temperature\n(°C)") +
  # coord_sf(crs = st_crs(3571),
  #          xlim = c(min(df_mesh_points$gridX3571),
  #                   max(df_mesh_points$gridX3571)),
  #          ylim = c(min(df_mesh_points$gridY3571),
  #                   max(df_mesh_points$gridY3571))) +
  theme(text = element_text(family = "serif"),
        strip.text.x = element_blank(),
        strip.background.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = unit(6, "pt")),
        axis.title = element_blank(),
        legend.title = element_text(size = unit(8, "pt")),
        legend.key.size = unit(6,"pt"),
        legend.text = element_text(size = unit(6, "pt"))) +
  facet_grid(. ~ fishery,
             labeller = labeller(fishery = fishery_labels))

plot_map_temp

ggplot(data = df_envdat) +
  geom_histogram(aes(x = temp_mean4)) +
  facet_grid(fishery ~ .)
  
# > 6.3: pH_thresh Map -----
plot_map_pH = basemap(limits = c(177, -162, 53, 63),
                      lon.interval = 5,
                      lat.interval = 2.5,
                      rotate = TRUE,
                      bathymetry = FALSE,
                      # bathy.style = "poly_blues",
                      # bathy.alpha = 0.5,
                      legends = FALSE,
                      land.col = "black") +
  geom_spatial_point(aes(x = gridX3571,
                         y = gridY3571,
                         color = factor(pH_thresh8)),
                     shape = 16,
                     size = 0.1,
                     crs = 3571,
                     data = df_envdat) +
  geom_spatial_text(aes(x = 177,
                        y = 61.5,
                        label = "C"),
                    color = "black",
                    fontface = "bold",
                    family = "serif",
                    crs = 4326,
                    data = data.frame(fishery = "QO18")) +
  scale_color_manual(values = c("magenta", "magenta4"),
                     labels = c("normal", "low pH")) +
  labs(color = "pH") +
  # coord_sf(crs = st_crs(3571),
  #          xlim = c(min(df_mesh_points$gridX3571),
  #                   max(df_mesh_points$gridX3571)),
  #          ylim = c(min(df_mesh_points$gridY3571),
  #                   max(df_mesh_points$gridY3571))) +
  theme(text = element_text(family = "serif"),
        strip.text.x = element_blank(),
        strip.background.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text = element_text(size = unit(6, "pt")),
        axis.title = element_blank(),
        legend.title = element_text(size = unit(8, "pt")),
        legend.key.size = unit(3,"pt"),
        legend.text = element_text(size = unit(6, "pt"))) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  facet_grid(. ~ fishery) 

plot_map_pH

# Combine maps
plot_map_prev + plot_map_temp + plot_map_pH +
  plot_layout(ncol = 1)

# export graph
ggsave("BES_plot_maps.png",
       device = "png",
       dpi = 300,
       height = 4.5,
       width = 6.5,
       units = "in")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# 7: Miscellaneous -------------------------------------------------------------
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# > 7.1: Does shell condition relate to size? ----
ggplot(data = dfm_QO1821) +
  geom_density(aes(x = size,
                   color = shell2,
                   fill = shell2),
               alpha = 0.5) +
  geom_boxplot(aes(x = size,
                   fill = shell2),
               orientation = "y",
               width = 0.01) +
  scale_color_manual(values = c("orange",
                                "orange3",
                                "orange4")) +
  scale_fill_manual(values = c("orange",
                               "orange3",
                               "orange4")) +
  labs(x = "carapace width (mm)") +
  guides(color = "none",
         fill = "none") +
  theme_bw() +
  theme(text = element_text(family = "serif")) +
  facet_grid(shell2 ~ ., 
             labeller = labeller(shell2 = c(new = "new",
                                            old = "old",
                                            v.old = "very old"))) 

# export graph
# ggsave("BES_plot_shellsize.png",
#        device = "png",
#        dpi = 300,
#        height = 3,
#        width = 4,
#        units = "in")


# Nope, not really.

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# END OF SCRIPT -----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

