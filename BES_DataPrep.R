# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

#  Black Eye Syndrome Spatiotemporal Analyses ----------------------------------
#  code by Reyn Yoshioka

# Data Prep Script -----
# Exploring potential environmental drivers of BES based on ADF&G crab observer
# data and Bering10K model hindcasts.

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# 1: Getting Started -----------------------------------------------------------
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# 1: Getting Started -----------------------------------------------------------
# > 1.1: Clear everything ----- 
# clear workspace
# rm(list=ls())
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

library(missForest) # Random Forest imputation of missing values
library(doParallel) # For registering parallel backend for missForest
library(doRNG) # set seed for reproducible parallel core use
library(foreach) # parallel for loops
library(RcppRoll) # rolling functions
library(ncdf4) # for netCDF files 
library(survival) # for neardate, nearest date matching
library(stringr) # character strings
library(patchwork) # arrange plots

# > > 1.2.1: here::i_am() -----
# setting working directory
here::i_am("SpatiotemporalBES/BES_DataPrep.R")

# > 1.3: Read in data ------
# > > 1.3.1: Observer data -----
df_QO18 = read.csv(here("QO18_932-Crab.csv"), na.strings = c("NA",""))
df_QO19 = read.csv(here("QO19_932-Crab.csv"), na.strings = c("NA",""))
df_QO20 = read.csv(here("QO20_932-Crab.csv"), na.strings = c("NA",""))
df_QO21 = read.csv(here("QO21_932-Crab.csv"), na.strings = c("NA",""))

df_EBS = read.csv(here("NOAA_FOSS_SnowCrab.csv"))

df_potsum = read.csv(here("SNOWCRAB-1990-2021_potsum.csv"))

df_retcatch = read.csv(here("retained_catch_all_fisheries.csv"))

# Combine fishery seasons
df_QO1821 = rbind(df_QO18,
                  df_QO19,
                  df_QO20,
                  df_QO21)

# > 1.4: Data management and cleaning ------------------------------------------

# Properly define NAs
df_QO1821$mi_lon = ifelse(df_QO1821$mi_lon != 0.0015,
                          df_QO1821$mi_lon,
                          NA)

df_QO1821$mi_lat = ifelse(df_QO1821$mi_lat != -0.0015,
                          df_QO1821$mi_lat,
                          NA)

df_QO1821$depth = ifelse(df_QO1821$depth != -9, # -9 is NA
                         df_QO1821$depth,
                         NA)

df_QO1821$legal = ifelse(df_QO1821$legal != -7, # -7 is NA
                         df_QO1821$legal,
                         NA)

df_QO1821$clutch = ifelse(df_QO1821$clutch != -7, # -7 is NA
                          df_QO1821$clutch,
                          NA)

df_QO1821$eggdev = ifelse(df_QO1821$eggdev != -7, # -7 is NA
                          df_QO1821$eggdev,
                          NA)

df_QO1821$clutchcon = ifelse(df_QO1821$clutchcon != -7, # -7 is NA
                             df_QO1821$clutchcon,
                             NA)

df_QO1821$maturity = ifelse(df_QO1821$maturity != -7, # -7 is NA
                            df_QO1821$maturity,
                            NA)

df_QO1821$sex = ifelse(df_QO1821$sex != 0, # 0 indicates an undetermined sex
                       df_QO1821$sex,
                       NA)

# Draw BES status from parasite code
df_QO1821$BES = ifelse(grepl("13", df_QO1821$parasite),
                       1,
                       0)

# redefine variables (attach labels to indices)
# create keys for matching
key_sex = data.frame(ind = c("1",
                             "2"),
                     sex = c("male",
                             "female"))
key_shell = data.frame(ind = c("0",
                               "1",
                               "2",
                               "3",
                               "4",
                               "5",
                               "9"),
                       shell = c("premolt.molting",
                                 "soft",
                                 "new",
                                 "old",
                                 "v.old",
                                 "v.v.old",
                                 "new.pliable"))
# apply keys to create new variables
df_QO1821$sex2 = key_sex$sex[match(as.character(df_QO1821$sex),
                                   key_sex$ind)]
df_QO1821$shell2 = key_shell$shell[match(as.character(df_QO1821$shell),
                                         key_shell$ind)]
# relevel shell2 to put in order (9 or new.pliable was out of place)
df_QO1821$shell2 = factor(df_QO1821$shell2,
                          levels = c("premolt.molting",
                                     "soft",
                                     "new.pliable",
                                     "new",
                                     "old",
                                     "v.old",
                                     "v.v.old"))

# Create a pot id
df_QO1821$pot = paste("pot",
                      df_QO1821$adfg,
                      df_QO1821$trip,
                      df_QO1821$spn,
                      sep = ".")

# filter df_potsum to only include snow crab (QO...)
df_potsum = df_potsum %>%
  filter(grepl("QO", fishery))

# filter df_potsum to only include relevant fishery years
df_potsum = df_potsum[df_potsum$fishery %in% c("QO18",
                                               "QO19",
                                               "QO20",
                                               "QO21"),]

# filter df_retcatch to only include snow crab (fishery = "snow")
df_retcatch = df_retcatch[df_retcatch$fishery == "snow" &
                            df_retcatch$dir_inc == "directed", ]

# Create a pot id 
df_potsum$pot = paste("pot",
                      df_potsum$adfg,
                      df_potsum$trip,
                      df_potsum$spn,
                      sep = ".")

# Create CPUE in df_retcatch
df_retcatch$CPUE = df_retcatch$ret_cat_crabs / df_retcatch$Fish_effort

# Standardize CPUE in df_retcatch
df_retcatch$CPUE_std = (df_retcatch$CPUE - mean(df_retcatch$CPUE)) / 
  sd(df_retcatch$CPUE)
hist(df_retcatch$CPUE_std)
M_CPUE = mean(df_retcatch$CPUE)
SD_CPUE = sd(df_retcatch$CPUE)

# Determine and create filter for statistical areas
# only show statistical areas where three or more vessels fished in a year
df_statarea = df_QO1821 |>
  group_by(fishery, statarea) |>
  summarize(vessels = length(unique(adfg)))
df_statarea$plot = df_statarea$vessels >= 3

# Check whether any records are in the irregular statistical areas around
# the Pribilof Islands
df_statarea[df_statarea$statarea == 705701,]$statarea
df_statarea[df_statarea$statarea == 705703,]$statarea
df_statarea[df_statarea$statarea == 695701,]$statarea
df_statarea[df_statarea$statarea == 695631,]$statarea
df_statarea[df_statarea$statarea == 695632,]$statarea

# > 1.5: Handling and imputation of missing data -----
# Trim df of all unnecessary variables and then impute
# See: https://rpubs.com/lmorgan95/MissForest
# Also see: https://www.kaggle.com/code/lmorgan95/missforest-the-best-imputation-algorithm/report
# Do impute: sex, size, shell condition
# Do not impute: location (mi_lon, mi_lat)

# > > 1.5.1: Select columns to retain in subset ("s" in "dfs") -----
dfs_QO1821 = df_QO1821[,c("fishery",
                          "sampdate",
                          "statarea",
                          "mi_lon",
                          "mi_lat",
                          "soaktime",
                          "sex2",
                          "size",
                          "shell2",
                          "pot",
                          "BES")]

summary(dfs_QO1821)
# > > 1.5.2: Drop missing locations -----
length(dfs_QO1821[is.na(dfs_QO1821$mi_lon),]$fishery) / 
  length(dfs_QO1821$fishery) * 100
# 662 missing location data, = 0.09%

dfs_QO1821 = dfs_QO1821[!is.na(dfs_QO1821$mi_lon),]
summary(dfs_QO1821)

# > > 1.5.3: Make sure classes are specified correctly -----
dfs_QO1821$fishery = as.factor(dfs_QO1821$fishery)
dfs_QO1821$sex2 = as.factor(dfs_QO1821$sex2)
dfs_QO1821$BES = as.factor(dfs_QO1821$BES)
summary(dfs_QO1821)

# Convert sampdate to a numeric
dfs_QO1821$sampdate = as.POSIXct(dfs_QO1821$sampdate,
                                 format = "%Y-%m-%d",
                                 origin = "1970-01-01 00:00:00",
                                 tz = "America/Anchorage")

dfs_QO1821$sampdate2 = as.numeric(dfs_QO1821$sampdate)

# > > 1.5.4: Make labels for which values are imputed ----
dfs_QO1821$imp_sex2 = ifelse(is.na(dfs_QO1821$sex2),
                             "sex2",
                             "")
dfs_QO1821$imp_shell2 = ifelse(is.na(dfs_QO1821$shell2),
                               "shell2",
                               "")
dfs_QO1821$imp_size = ifelse(is.na(dfs_QO1821$size),
                             "size",
                             "")
dfs_QO1821$imp = paste(dfs_QO1821$imp_sex2,
                       dfs_QO1821$imp_size,
                       dfs_QO1821$imp_shell2,
                       sep = "")
dfs_QO1821[dfs_QO1821$imp == "",]$imp = "none"
dfs_QO1821$imp_sex2 = NULL
dfs_QO1821$imp_shell2 = NULL
dfs_QO1821$imp_size = NULL

# Tally missing/imputed crab
dfs_QO1821 |>
  summarize(total = length(sex2),
            sex2 = length(imp[grepl(imp, pattern = "sex2")]),
            shell2 = length(imp[grepl(imp, pattern = "shell2")]),
            size = length(imp[grepl(imp, pattern = "size")]),
            shell2size = length(imp[grepl(imp, pattern = "shell2") &
                                      grepl(imp, pattern = "size")])) |>
  mutate(shell2 = shell2 - shell2size,
         size = size - shell2size,
         perc_sex2 = sex2 / total * 100,
         perc_shell2 = shell2 / total * 100,
         perc_size = size / total * 100,
         perc_shell2size = shell2size / total * 100)

# > > 1.5.5: Data to be separated and reattached -----
# Separate, retain, and remove BES to be reattached to imputed data later
dfs_QO1821SBES = dfs_QO1821$BES # Vector of imp data
dfs_QO1821$BES = NULL # Remove imp from df to be imputed

# Separate, retain, and remove imp to be reattached to imputed data later
dfs_QO1821Simp = dfs_QO1821$imp # Vector of imp data
dfs_QO1821$imp = NULL # Remove imp from df to be imputed

# Separate, retain, and remove pot to be reattached to imputed data later
dfs_QO1821Spot = dfs_QO1821$pot # Vector of pot data
dfs_QO1821$pot = NULL # Remove pot from df to be imputed

# Separate, retain, and remove sampdate to be reattached to imputed data later
dfs_QO1821Ssampdate = dfs_QO1821$sampdate # Vector of sampdate data
dfs_QO1821$sampdate = NULL # Remove sampdate from df to be imputed

# Separate, retain, and remove statarea to be reattached to imputed data later
dfs_QO1821Sstatarea = dfs_QO1821$statarea # Vector of statarea data
dfs_QO1821$statarea = NULL # Remove statarea from df to be imputed

# # > > 1.5.6: Imputation -----
# # set up parallel
# doParallel::registerDoParallel(cores = 8)
# # Originally, we set seed arbitrarily to 36, but we leave it blank to be
# # actually random:
# doRNG::registerDoRNG()
# 
# # missForest
# # From a previous run, we found 6 iterations to be sufficient
# print(Sys.time())
# # 09:36:53 EDT
# dfs_QO1821.imp = missForest(dfs_QO1821,
#                             parallelize = 'forests',
#                             variablewise = TRUE,
#                             ntree = 100,
#                             maxiter = 6,
#                             verbose = TRUE)
# print(Sys.time())
# # 10:41:05 EDT
# 
# dfs_QO1821_imp = dfs_QO1821.imp$ximp
# 
# dfs_QO1821.imp$OOBerror
# 
# # Save data (imputation takes a while)
# save.image(paste0("ADFG_CrabData_",
#            format(Sys.time(),
#                   "%Y%m%d_%H%M"),
#            ".RData"))
# 
# write.csv(dfs_QO1821_imp,
#           paste0("dfs_QO1821_imp_",
#                  format(Sys.time(),
#                         "%Y%m%d_%H%M"),
#                  ".csv"))
# 
# write.csv(dfs_QO1821.imp$OOBerror,
#           paste0("dfs_QO1821_imp_OOBerr_",
#                  format(Sys.time(),
#                         "%Y%m%d_%H%M"),
#                  ".csv"))

# > > 1.5.7: Read imputed data back in ("i" in "dfi") -----
dfi_QO1821 = read.csv(here("dfs_QO1821_imp_20240426_1042.csv"))
dfi_QO1821$X = NULL

dfi_QO1821$shell2 = factor(dfi_QO1821$shell2,
                           levels = c("premolt.molting",
                                      "soft",
                                      "new.pliable",
                                      "new",
                                      "old",
                                      "v.old",
                                      "v.v.old"))

dfi_QO1821$sex2 = factor(dfi_QO1821$sex2)

# > > 1.5.8: Reattach separated data -----
# Reattach BES
dfi_QO1821$BES = dfs_QO1821SBES

# Reattach imp identification
dfi_QO1821$imp = dfs_QO1821Simp

# Remove dfs_QO1821Simp 
rm(dfs_QO1821Simp)

# Reattach pot identification
dfi_QO1821$pot = dfs_QO1821Spot

# Remove dfs_QO1821Spot
rm(dfs_QO1821Spot)

# Reattach sampdate
dfi_QO1821$sampdate = dfs_QO1821Ssampdate

# Remove dfs_QO1821Ssampdate 
rm(dfs_QO1821Ssampdate)

# Reattach statarea
dfi_QO1821$statarea = dfs_QO1821Sstatarea

# Remove dfs_QO1821Sstatarea 
rm(dfs_QO1821Sstatarea)

# > > 1.5.9: Index Seasons (QO18 --> 1, etc.) -----
dfi_QO1821$fishery_ind = ifelse(dfi_QO1821$fishery == "QO18",
                                1,
                                ifelse(dfi_QO1821$fishery == "QO19",
                                       2,
                                       ifelse(dfi_QO1821$fishery == "QO20",
                                              3,
                                              4)))

# > > 1.5.10: Limit data to males only -----
length(dfi_QO1821[dfi_QO1821$sex2 == "female",]$sex2) /
  length(dfi_QO1821$sex2) * 100
# only 0.3% of the crab are female
dfm_QO1821 = dfi_QO1821[dfi_QO1821$sex2 == "male",]

# > > 1.5.11: Standardize size -----
M_size = mean(dfm_QO1821$size)
# 98.66659
SD_size = sd(dfm_QO1821$size)
# 10.43882

dfm_QO1821$size_std = (dfm_QO1821$size - M_size) / SD_size

# Make a key for trips to vessel
vessel_trip = unique(paste0("V",
                            paste(df_QO1821$adfg,
                                  df_QO1821$trip,
                                  sep = "-")))
trip_key = data.frame(vessel = str_split_i(vessel_trip,
                                           "-",
                                           1),
                      trip = str_split_i(vessel_trip,
                                         "-",
                                         2))

write.csv(file = "trip_key.csv",
          trip_key)

# 2: Spatiotemporal structure and ROMS data matching -----
# > 2.1: Initial coordinate conversion for spatial structure -----
# Convert coordinates to WGS 84 / North Pole LAEA Bering Sea (3571) for meters,
# assuming WGS 84 (4326) in decimal degrees as original coordinate system
# This means that spatial measures will be in meters

# > > 2.1.1: Set coordinates as simple feature points -----
dfm_QO1821 = dfm_QO1821 %>%
  rowwise() %>%
  mutate(coord3571 = list(st_point(c(mi_lon,
                                     mi_lat))))

# > > 2.1.2: CRS Conversion from 4326 to 3571 -----
# (This takes a bit of time with all the crab...)
dfm_QO1821 = dfm_QO1821 %>%
  rowwise() %>%
  mutate(coord3571 = coord3571 %>%
           st_sfc(crs = 4326) %>%
           st_transform(crs = 3571) %>%
           st_coordinates() %>%
           list(),
         X3571 = unlist(coord3571)[1],
         Y3571 = unlist(coord3571)[2])

# > > 2.1.3: View sample locations -----
basemap(limits = c(-180, -160, 55, 65),
        rotate = TRUE) +
  geom_spatial_point(aes(x = X3571,
                         y = Y3571,
                         # size = BES,
                         color = BES),
                     alpha = 0.1,
                     data = dfm_QO1821,
                     crs = 3571) +
  scale_color_viridis_d(direction = -1) +
  theme_bw() +
  theme() +
  facet_wrap(. ~ fishery)

# > 2.2: Form analysis boundary ------
# mesh boundaries
obj_mesh_bounds = 
  fm_extensions(unique(cbind(dfm_QO1821$X3571,
                             dfm_QO1821$Y3571)),
                convex = -0.15)
# NOTE: We use a looser boundary so that when graphing, the data points are
# more ambiguous. The mesh used for analysis is much tighter.

poly_mesh_bounds = 
  obj_mesh_bounds[[1]] |>
  st_sfc(crs = 3571) |>
  st_transform(crs = 4326) # convert to 4326 so later matching is easier

# > 2.3: Prepare df_potsum and df_retcatch -----
# > >  2.3.1: Merge df_retcatch to df_potsum -----
df_retcatch$fishery2 = df_retcatch$fishery
df_retcatch$fishery = paste0("QO",
                             str_sub(df_retcatch$year, -2))
df_potsum = left_join(df_potsum,
                      df_retcatch[,
                                  c("fishery",
                                    "CPUE",
                                    "CPUE_std")],
                      by = join_by(fishery))

df_potsum = df_potsum %>%
  rowwise() %>%
  mutate(catch = sum(legal_ret,
                     legal_nr,
                     legal_ur,
                     sublegal),
         p_retcatch = legal_ret / catch)

# > > 2.3.2: Coordinate conversion for df_potsum -----
df_potsum = df_potsum %>%
  rowwise() %>%
  mutate(coord3571 = list(st_point(c(longitude,
                                     latitude))))

# > > > 2.3.2.1: CRS Conversion from 4326 to 3571 -----
# (This takes a bit of time with all the pots...)
df_potsum = df_potsum %>%
  rowwise() %>%
  mutate(coord4326 = coord3571 %>%
           st_sfc(crs = 4326) %>%
           st_coordinates() %>%
           list(),
         X4326 = unlist(coord4326)[1],
         Y4326 = unlist(coord4326)[2],
         coord3571 = coord3571 %>%
           st_sfc(crs = 4326) %>%
           st_transform(crs = 3571) %>%
           st_coordinates() %>%
           list(),
         X3571 = unlist(coord3571)[1],
         Y3571 = unlist(coord3571)[2])

# > > 2.3.3: Time for df_potsum -----
df_potsum$sampdate = as.POSIXct(df_potsum$sampdate,
                                format = "%m-%d-%Y",
                                origin = "1970-01-01 00:00:00",
                                tz = "America/Anchorage")

df_potsum$sampdate2 = as.numeric(df_potsum$sampdate)

# 3: Read in and prepare ROMS data -----
# > 3.1: Define grid ----
# NOTE: This is not the best way, but we need to access data in a different
# directory:
B10Kgrid = 
  nc_open(
    paste0(here() |> str_remove("/2_BES_ADFG_CrabData"), # remove current
           "/X3_Bering10K/", # tack on where the data live
           "Bering10K_extended_grid.nc")) # the file

B10Kgrid_latrho = ncvar_get(B10Kgrid, "lat_rho")
B10Kgrid_lonrho = ncvar_get(B10Kgrid, "lon_rho")

B10Kgrid_lat = as.vector(B10Kgrid_latrho)
B10Kgrid_lon = as.vector(B10Kgrid_lonrho)

basegrid = data.frame(lon = B10Kgrid_lon,
                      lat = B10Kgrid_lat)
basegrid = basegrid %>%
  rowwise() %>%
  mutate(coord4326 = list(st_point(c(lon,
                                     lat))))

basegrid$coord4326 = st_sfc(basegrid$coord4326,
                            crs = 4326)

basegrid$index = lengths(st_intersects(basegrid$coord4326,
                                       poly_mesh_bounds,
                                       sparse = TRUE)) > 0

plot_mod_region = basemap(limits = c(-180, -160, 54, 65),
                          rotate = TRUE) +
  geom_spatial_point(aes(x = lon,
                         y = lat),
                     crs = 4326,
                     data = basegrid[basegrid$index,])
plot_mod_region

length(basegrid$index)
basegrid$gridind = 1:length(basegrid$index)

# Set the basegrid clipped to the prediction bounds and prepare for coords
# This is used later to clip ROMS data
basegrid_masked = basegrid[basegrid$index,]

basegrid_masked = basegrid_masked %>%
  rowwise() %>%
  mutate(coord3571 = list(st_point(c(lon,
                                     lat))))

# convert CRS and extract the 3571 coordinates
basegrid_masked = basegrid_masked %>%
  rowwise() %>%
  mutate(coord3571 = coord3571 %>%
           st_sfc(crs = 4326) %>%
           st_transform(crs = 3571) %>%
           st_coordinates() %>%
           list(),
         gridX3571 = unlist(coord3571)[1],
         gridY3571 = unlist(coord3571)[2])

# > > 3.1.1: Assign statistical areas to grid ----
statarea_grid = matrix(c(-180, -156, -156, -180, 
                         62, 62, 53.5, 53.5),
                       ncol = 2)
statarea_grid = st_multipoint(statarea_grid, dim = "XY")
plot(statarea_grid)

# obtain centers for later plotting
statarea_grid_cent = st_make_grid(statarea_grid,
                                  cellsize = c(1, 0.5),
                                  what = "centers",
                                  crs = 4326)
# also convert to 3571
statarea_grid_cent2 = statarea_grid_cent |>
  st_transform(crs = 3571)

# polygons for grid - statarea matching
statarea_grid = st_make_grid(statarea_grid,
                             cellsize = c(1, 0.5),
                             what = "polygons",
                             crs = 4326)
# temporary IDs
statarea_grid = cbind(statarea_grid, data.frame(ID = paste0("statarea.",
                                                            1:length(statarea_grid))))
plot(statarea_grid$geometry)
plot(statarea_grid_cent)

# Extract lower right coordinates, which are the basis of the naming
statarea_grid$xlr = lapply(1:length(statarea_grid$geometry), function(x) statarea_grid[[1]][[x]][[1]][2,1])
statarea_grid$ylr = lapply(1:length(statarea_grid$geometry), function(x) statarea_grid[[1]][[x]][[1]][2,2])
statarea_grid$sa1 = gsub("-1", "", statarea_grid$xlr)
statarea_grid$sa2 = as.character(floor(as.numeric(statarea_grid$ylr)))
statarea_grid$sa3 = ifelse(as.numeric(statarea_grid$ylr) - 
                             floor(as.numeric(statarea_grid$ylr)) == 0.5,
                           "30",
                           "00")
statarea_grid$statarea = paste0(statarea_grid$sa1,
                                statarea_grid$sa2,
                                statarea_grid$sa3)
# remove unnecessary variables
statarea_grid$ID = NULL
statarea_grid$xlr = NULL
statarea_grid$ylr = NULL
statarea_grid$sa1 = NULL
statarea_grid$sa2 = NULL
statarea_grid$sa3 = NULL

# attach statistical area centers to the key dataframe (df_statarea)
statarea_grid_cent = cbind(statarea_grid_cent,
                           statarea_grid_cent2,
                           statarea_grid$statarea)
colnames(statarea_grid_cent) = c("coord4326",
                                 "coord3571",
                                 "statarea")

statarea_grid_cent = statarea_grid_cent |>
  data.frame() |>
  rowwise() |>
  mutate(saX4326 = unlist(coord4326)[1],
         saY4326 = unlist(coord4326)[2],
         saX3571 = unlist(coord3571)[1],
         saY3571 = unlist(coord3571)[2])
statarea_grid_cent$statarea = as.numeric(as.character(statarea_grid_cent$statarea))

df_statarea = left_join(df_statarea,
                        statarea_grid_cent,
                        by = join_by(statarea))
df_statarea$coord4326 = NULL
df_statarea$coord3571 = NULL

# export
write.csv(df_statarea, "df_statarea.csv")

# match the statistical areas to basegrid_masked (used later for df_envdat)
basegrid_masked = st_join(st_sf(basegrid_masked),
                          st_sf(statarea_grid),
                          join = st_intersects)

basegrid_masked_lab = basegrid_masked |>
  group_by(statarea) |>
  summarize(gridX3571 = mean(gridX3571),
            gridY3571 = mean(gridY3571))

# plot it out
ggplot(aes(x = gridX3571,
           y = gridY3571,
           fill = statarea,
           color = (as.numeric(statarea)/10) %% 2),
       data = basegrid_masked) +
  geom_point(shape = 21,
             size = 2,
             alpha = 0.5) +
  geom_text(aes(x = gridX3571,
                y = gridY3571,
                label = statarea),
            color = "red",
            data = basegrid_masked_lab) +
  theme_bw() +
  theme(legend.position = "none")

# > 3.2: Match observer and ROMS data in space ----
# Set CRS to 4326 to prepare for matching
dfm_QO1821 = dfm_QO1821 %>%
  rowwise() %>%
  mutate(coord4326 = list(st_point(c(mi_lon,
                                     mi_lat))))
dfm_QO1821$coord4326 = st_sfc(dfm_QO1821$coord4326,
                              crs = 4326)

# Provides index of grid_masked that is the closest to the crab point
banana = st_nearest_feature(dfm_QO1821$coord4326,
                            basegrid_masked$coord4326)

# View matching of ROMS grid and observer data
basemap(limits = c(-180, -160, 54, 65),
        rotate = TRUE) +
  geom_spatial_point(aes(x = lon,
                         y = lat),
                     color = "grey",
                     crs = 4326,
                     data = basegrid_masked) +
  geom_spatial_point(aes(x = mi_lon,
                         y = mi_lat),
                     color = "blue",
                     crs = 4326,
                     data = dfm_QO1821 |>
                       select(mi_lon,
                              mi_lat) |>
                       unique()) +
  geom_spatial_point(aes(x = lon,
                         y = lat),
                     color = "red",
                     shape = 13,
                     size = 3,
                     crs = 4326,
                     data = basegrid_masked[banana,]) +
  geom_spatial_point(aes(x = saX4326,
                         y = saY4326),
                     color = "green",
                     crs = 4326,
                     data = df_statarea |>
                       filter(plot)) 

# attach basegrid data to dfm_QO1821
dfm_QO1821$gridind = basegrid_masked[banana,]$gridind
dfm_QO1821$lon_rho = basegrid_masked[banana,]$lon
dfm_QO1821$lat_rho = basegrid_masked[banana,]$lat
dfm_QO1821$gridX3571 = basegrid_masked[banana,]$gridX3571
dfm_QO1821$gridY3571 = basegrid_masked[banana,]$gridY3571

# > 3.3: Match df_potsum and ROMS data in space -----
df_potsum = df_potsum %>%
  rowwise() %>%
  mutate(coord4326 = list(st_point(c(longitude,
                                     latitude))))
df_potsum = df_potsum %>%
  rowwise() %>%
  mutate(coord4326 = coord4326 %>%
           st_sfc(crs = 4326))

# Using basegrid here instead of basegrid_masked allows us to filter out the 
# locations outside of the modeling range
banana = st_nearest_feature(df_potsum$coord4326,
                            basegrid$coord4326)

# attach basegrid data to dfm_QO1821
df_potsum$gridind = basegrid[banana,]$gridind
df_potsum$lon_rho = basegrid[banana,]$lon
df_potsum$lat_rho = basegrid[banana,]$lat
df_potsum$index = basegrid[banana,]$index

# filter df_potsum to within prediction grid
df_potsum = df_potsum[df_potsum$index,]

# 4: Environmental Variables -----
# > 4.1: Extract netcdf4 data -----
# Data in 5-year packets
# Temperature
temp_1519 = 
  nc_open(
    paste0(here() |> str_remove("/2_BES_ADFG_CrabData"), # remove current
           "/X3_Bering10K/", # tack on where the data live
           "B10K-K20P19_CORECFS_2015-2019_average_temp_bottom5m.nc")) # the file
temp_2024 = 
  nc_open(
    paste0(here() |> str_remove("/2_BES_ADFG_CrabData"), # remove current
           "/X3_Bering10K/", # tack on where the data live
           "B10K-K20P19_CORECFS_2020-2024_average_temp_bottom5m.nc")) # the file

# Dissolved oxygen
O2_1519 = 
  nc_open(
    paste0(here() |> str_remove("/2_BES_ADFG_CrabData"), # remove current
           "/X3_Bering10K/", # tack on where the data live
           "B10K-K20P19_CORECFS_2015-2019_average_oxygen_bottom5m.nc")) # the file
O2_2024 = 
  nc_open(
    paste0(here() |> str_remove("/2_BES_ADFG_CrabData"), # remove current
           "/X3_Bering10K/", # tack on where the data live
           "B10K-K20P19_CORECFS_2020-2024_average_oxygen_bottom5m.nc")) # the file

# pH
pH_1519 = 
  nc_open(
    paste0(here() |> str_remove("/2_BES_ADFG_CrabData"), # remove current
           "/X3_Bering10K/", # tack on where the data live
           "B10K-K20P19_CORECFS_2015-2019_average_pH_bottom5m.nc")) # the file
pH_2024 = 
  nc_open(
    paste0(here() |> str_remove("/2_BES_ADFG_CrabData"), # remove current
           "/X3_Bering10K/", # tack on where the data live
           "B10K-K20P19_CORECFS_2020-2024_average_pH_bottom5m.nc")) # the file

# > 4.2: Extract data from netcdf4 -----
# Temperature
# 2015-2019
temp_1519_array = ncvar_get(temp_1519,"temp")
fillval = ncatt_get(temp_1519, "temp", "_FillValue")
temp_1519_latrho = ncvar_get(temp_1519, "lat_rho")
dim(temp_1519_latrho)
temp_1519_lonrho = ncvar_get(temp_1519, "lon_rho")
dim(temp_1519_lonrho)
temp_1519_time = ncvar_get(temp_1519, "ocean_time")
dim(temp_1519_array)
temp_1519_array[temp_1519_array == fillval$value] = NA
temp_1519_lat = as.vector(temp_1519_latrho)
temp_1519_lon = as.vector(temp_1519_lonrho)

# 2020-2024 (actually to 2022)
temp_2024_array = ncvar_get(temp_2024,"temp")
fillval = ncatt_get(temp_2024, "temp", "_FillValue")
temp_2024_latrho = ncvar_get(temp_2024, "lat_rho")
dim(temp_2024_latrho)
temp_2024_lonrho = ncvar_get(temp_2024, "lon_rho")
dim(temp_2024_lonrho)
temp_2024_time = ncvar_get(temp_2024, "ocean_time")
dim(temp_2024_array)
temp_2024_array[temp_2024_array == fillval$value] = NA
temp_2024_lat = as.vector(temp_2024_latrho)
temp_2024_lon = as.vector(temp_2024_lonrho)

# Dissolved oxygen
# 2015-2019
O2_1519_array = ncvar_get(O2_1519,"oxygen")
fillval = ncatt_get(O2_1519, "oxygen", "_FillValue")
O2_1519_latrho = ncvar_get(O2_1519, "lat_rho")
dim(O2_1519_latrho)
O2_1519_lonrho = ncvar_get(O2_1519, "lon_rho")
dim(O2_1519_lonrho)
O2_1519_time = ncvar_get(O2_1519, "ocean_time")
dim(O2_1519_array)
O2_1519_array[O2_1519_array == fillval$value] = NA
O2_1519_lat = as.vector(O2_1519_latrho)
O2_1519_lon = as.vector(O2_1519_lonrho)
# 2020-2024
O2_2024_array = ncvar_get(O2_2024,"oxygen")
fillval = ncatt_get(O2_2024, "oxygen", "_FillValue")
O2_2024_latrho = ncvar_get(O2_2024, "lat_rho")
dim(O2_2024_latrho)
O2_2024_lonrho = ncvar_get(O2_2024, "lon_rho")
dim(O2_2024_lonrho)
O2_2024_time = ncvar_get(O2_2024, "ocean_time")
dim(O2_2024_array)
O2_2024_array[O2_2024_array == fillval$value] = NA
O2_2024_lat = as.vector(O2_2024_latrho)
O2_2024_lon = as.vector(O2_2024_lonrho)

# pH
# 2015-2019
pH_1519_array = ncvar_get(pH_1519,"pH")
fillval = ncatt_get(pH_1519, "pH", "_FillValue")
pH_1519_latrho = ncvar_get(pH_1519, "lat_rho")
dim(pH_1519_latrho)
pH_1519_lonrho = ncvar_get(pH_1519, "lon_rho")
dim(pH_1519_lonrho)
pH_1519_time = ncvar_get(pH_1519, "ocean_time")
dim(pH_1519_array)
pH_1519_array[pH_1519_array == fillval$value] = NA
pH_1519_lat = as.vector(pH_1519_latrho)
pH_1519_lon = as.vector(pH_1519_lonrho)
# 2020-2024
pH_2024_array = ncvar_get(pH_2024,"pH")
fillval = ncatt_get(pH_2024, "pH", "_FillValue")
pH_2024_latrho = ncvar_get(pH_2024, "lat_rho")
dim(pH_2024_latrho)
pH_2024_lonrho = ncvar_get(pH_2024, "lon_rho")
dim(pH_2024_lonrho)
pH_2024_time = ncvar_get(pH_2024, "ocean_time")
dim(pH_2024_array)
pH_2024_array[pH_2024_array == fillval$value] = NA
pH_2024_lat = as.vector(pH_2024_latrho)
pH_2024_lon = as.vector(pH_2024_lonrho)

# > 4.3: Extract data matching grid -----
# Temperature 
# Create temporary data frame called "banana"
banana = data.frame(lat = NA,
                    lon = NA,
                    temp = NA,
                    ind = NA,
                    gridind = NA,
                    t = NA)

# Loop through temp data and extract based on basegrid$index, which defines
# the overlap between predBounds and the ROMS grid
# Note that the temperature data will hold the "ind" for time index
for (i in 1:length(temp_1519_time)) {
  banana2 = data.frame(lat = temp_1519_lat[basegrid$index],
                       lon = temp_1519_lon[basegrid$index],
                       temp = as.vector(temp_1519_array[, , i])[basegrid$index],
                       ind = i,
                       gridind = basegrid[basegrid$index,]$gridind,
                       t = temp_1519_time[i])
  banana = rbind(banana, banana2)
}

indlength = length(unique(banana[!is.na(banana$ind),]$ind))

for (i in 1:length(temp_2024_time)) {
  banana2 = data.frame(lat = temp_2024_lat[basegrid$index],
                       lon = temp_2024_lon[basegrid$index],
                       temp = as.vector(temp_2024_array[, , i])[basegrid$index],
                       ind = indlength + i,
                       gridind = basegrid[basegrid$index,]$gridind,
                       t = temp_2024_time[i])
  banana = rbind(banana, banana2)
}

# indlength is the number of time points
indlength = length(unique(banana[!is.na(banana$ind),]$ind))

# Remove rows with missing temperature data
banana = banana[!is.na(banana$temp),]

# Reassign banana to B10K_btemp and remove banana
B10K_btemp = banana
rm(banana)

# Oxygen
# Create temporary data frame called "banana"
banana = data.frame(lat = NA,
                    lon = NA,
                    O2 = NA,
                    gridind = NA,
                    t = NA)

# Loop through temp data and extract based on basegrid$index, which defines
# the overlap between predBounds and the ROMS grid
for (i in 1:length(O2_1519_time)) {
  banana2 = data.frame(lat = O2_1519_lat[basegrid$index],
                       lon = O2_1519_lon[basegrid$index],
                       O2 = as.vector(O2_1519_array[, , i])[basegrid$index],
                       gridind = basegrid[basegrid$index,]$gridind,
                       t = O2_1519_time[i])
  banana = rbind(banana, banana2)
}

for (i in 1:length(O2_2024_time)) {
  banana2 = data.frame(lat = O2_2024_lat[basegrid$index],
                       lon = O2_2024_lon[basegrid$index],
                       O2 = as.vector(O2_2024_array[, , i])[basegrid$index],
                       gridind = basegrid[basegrid$index,]$gridind,
                       t = O2_2024_time[i])
  banana = rbind(banana, banana2)
}

# Remove rows with missing O2 data
banana = banana[!is.na(banana$O2),]

B10K_bO2 = banana

# Convert O2 in mmol m^-3 to mg/L
B10K_bO2$O2_mmolpm3 = B10K_bO2$O2
B10K_bO2$O2 = B10K_bO2$O2_mmolpm3 / 1000 * 15.999 * 2
#                               / L/m^-3 * mg/mmol  

# pH
# Create temporary data frame called "banana"
banana = data.frame(lat = NA,
                    lon = NA,
                    pH = NA,
                    gridind = NA,
                    t = NA)

# Loop through temp data and extract based on basegrid$index, which defines
# the overlap between predBounds and the ROMS grid
for (i in 1:length(pH_1519_time)) {
  banana2 = data.frame(lat = pH_1519_lat[basegrid$index],
                       lon = pH_1519_lon[basegrid$index],
                       pH = as.vector(pH_1519_array[, , i])[basegrid$index],
                       gridind = basegrid[basegrid$index,]$gridind,
                       t = pH_1519_time[i])
  banana = rbind(banana, banana2)
}

for (i in 1:length(pH_2024_time)) {
  banana2 = data.frame(lat = pH_2024_lat[basegrid$index],
                       lon = pH_2024_lon[basegrid$index],
                       pH = as.vector(pH_2024_array[, , i])[basegrid$index],
                       gridind = basegrid[basegrid$index,]$gridind,
                       t = pH_2024_time[i])
  banana = rbind(banana, banana2)
}

# Remove rows with missing pH data
banana = banana[!is.na(banana$pH),]

B10K_bpH = banana

# > 4.4: Join environmental data -----
B10K_benv = left_join(B10K_btemp,
                      B10K_bO2,
                      by = join_by(lat, lon, gridind, t))
B10K_benv = left_join(B10K_benv,
                      B10K_bpH,
                      by = join_by(lat, lon, gridind, t))

# Convert time from seconds since 1900 Jan 01 to POSIXct
B10K_benv$t2 = as.POSIXct(B10K_benv$t,
                          origin = "1900-01-01 00:00:00",
                          tz = "UTC")

# > 4.5: Rolling averages ----
banana = B10K_benv %>%
  group_by(gridind) %>%
  arrange(ind) %>%
  mutate(temp_mean24 = roll_mean(x = lag(temp,
                                         n = 1), 
                                 n = 24,
                                 align = "right",
                                 fill = NA,
                                 na.rm = TRUE),
         temp_mean20 = roll_mean(x = lag(temp,
                                         n = 1), 
                                 n = 20,
                                 align = "right",
                                 fill = NA,
                                 na.rm = TRUE),
         temp_mean16 = roll_mean(x = lag(temp,
                                         n = 1), 
                                 n = 16,
                                 align = "right",
                                 fill = NA,
                                 na.rm = TRUE),
         temp_mean12 = roll_mean(x = lag(temp,
                                         n = 1), 
                                 n = 12,
                                 align = "right",
                                 fill = NA,
                                 na.rm = TRUE),
         temp_mean8 = roll_mean(x = lag(temp,
                                        n = 1), 
                                n = 8,
                                align = "right",
                                fill = NA,
                                na.rm = TRUE),
         temp_mean4 = roll_mean(x = lag(temp,
                                        n = 1), 
                                n = 4,
                                align = "right",
                                fill = NA,
                                na.rm = TRUE),
         # "T" is for Time (eventually week above threshold)
         O2_threshT = ifelse(O2 > 3.21, # If O2 above threshold
                             0, # assign 0 (normoxic)
                             1), # assign 1 (hypoxic)
         # "U" is for Unit (eventually deficit O2 below threshold)
         O2_threshU = ifelse(O2 > 3.21, # If O2 above threshold
                             0, # assign 0 (normoxic)
                             3.21 - O2), # calculate difference
         # Averaged O2
         O2_mean24 = roll_mean(x = lag(O2,
                                       n = 1), 
                               n = 24,
                               align = "right",
                               fill = NA,
                               na.rm = TRUE),
         O2_mean20 = roll_mean(x = lag(O2,
                                       n = 1), 
                               n = 20,
                               align = "right",
                               fill = NA,
                               na.rm = TRUE),
         O2_mean16 = roll_mean(x = lag(O2,
                                       n = 1), 
                               n = 16,
                               align = "right",
                               fill = NA,
                               na.rm = TRUE),
         O2_mean12 = roll_mean(x = lag(O2,
                                       n = 1), 
                               n = 12,
                               align = "right",
                               fill = NA,
                               na.rm = TRUE),
         O2_mean8 = roll_mean(x = lag(O2,
                                      n = 1), 
                              n = 8,
                              align = "right",
                              fill = NA,
                              na.rm = TRUE),
         O2_mean4 = roll_mean(x = lag(O2,
                                      n = 1), 
                              n = 4,
                              align = "right",
                              fill = NA,
                              na.rm = TRUE),
         # Proportion of weeks below O2 threshold
         O2_threshT24 = roll_mean(x = lag(O2_threshT,
                                          n = 1), 
                                  n = 24,
                                  align = "right",
                                  fill = NA,
                                  na.rm = TRUE),
         O2_threshT20 = roll_mean(x = lag(O2_threshT,
                                          n = 1), 
                                  n = 20,
                                  align = "right",
                                  fill = NA,
                                  na.rm = TRUE),
         O2_threshT16 = roll_mean(x = lag(O2_threshT,
                                          n = 1), 
                                  n = 16,
                                  align = "right",
                                  fill = NA,
                                  na.rm = TRUE),
         O2_threshT12 = roll_mean(x = lag(O2_threshT,
                                          n = 1), 
                                  n = 12,
                                  align = "right",
                                  fill = NA,
                                  na.rm = TRUE),
         O2_threshT8 = roll_mean(x = lag(O2_threshT,
                                         n = 1), 
                                 n = 8,
                                 align = "right",
                                 fill = NA,
                                 na.rm = TRUE),
         O2_threshT4 = roll_mean(x = lag(O2_threshT,
                                         n = 1), 
                                 n = 4,
                                 align = "right",
                                 fill = NA,
                                 na.rm = TRUE),
         # Accumulated O2 deficit
         O2_threshU24 = roll_sum(x = lag(O2_threshU,
                                         n = 1), 
                                 n = 24,
                                 align = "right",
                                 fill = NA,
                                 na.rm = TRUE),
         O2_threshU20 = roll_sum(x = lag(O2_threshU,
                                         n = 1), 
                                 n = 20,
                                 align = "right",
                                 fill = NA,
                                 na.rm = TRUE),
         O2_threshU16 = roll_sum(x = lag(O2_threshU,
                                         n = 1), 
                                 n = 16,
                                 align = "right",
                                 fill = NA,
                                 na.rm = TRUE),
         O2_threshU12 = roll_sum(x = lag(O2_threshU,
                                         n = 1), 
                                 n = 12,
                                 align = "right",
                                 fill = NA,
                                 na.rm = TRUE),
         O2_threshU8 = roll_sum(x = lag(O2_threshU,
                                        n = 1), 
                                n = 8,
                                align = "right",
                                fill = NA,
                                na.rm = TRUE),
         O2_threshU4 = roll_sum(x = lag(O2_threshU,
                                        n = 1), 
                                n = 4,
                                align = "right",
                                fill = NA,
                                na.rm = TRUE),
         Hydronium = 10 ^ (-pH), # pH is log10 scaled, take it back to [H3O+]
         # "T" is for Time (eventually week above threshold)
         pH_threshT = ifelse(pH > 7.75, # If pH above threshold, 
                             0, # assign 0 (normal)
                             1), # assign 1 (acidic)
         # "U" is for Unit (eventually excess pH below threshold)
         pH_threshU = ifelse(pH > 7.75, # If pH threshold is surpassed...
                             0, # assign 0
                             Hydronium - (10 ^ (-7.75))), # Calc excess [H3O+]
         # Averaged pH, first put on [H3O+] scale
         pH_mean24 = -log10(roll_mean(x = lag(Hydronium, 
                                              n = 1), 
                                      n = 24,
                                      align = "right",
                                      fill = NA,
                                      na.rm = TRUE)),
         pH_mean20 = -log10(roll_mean(x = lag(Hydronium,
                                              n = 1), 
                                      n = 20,
                                      align = "right",
                                      fill = NA,
                                      na.rm = TRUE)),
         pH_mean16 = -log10(roll_mean(x = lag(Hydronium,
                                              n = 1), 
                                      n = 16,
                                      align = "right",
                                      fill = NA,
                                      na.rm = TRUE)),
         pH_mean12 = -log10(roll_mean(x = lag(Hydronium,
                                              n = 1), 
                                      n = 12,
                                      align = "right",
                                      fill = NA,
                                      na.rm = TRUE)),
         pH_mean8 = -log10(roll_mean(x = lag(Hydronium,
                                             n = 1), 
                                     n = 8,
                                     align = "right",
                                     fill = NA,
                                     na.rm = TRUE)),
         pH_mean4 = -log10(roll_mean(x = lag(Hydronium,
                                             n = 1), 
                                     n = 4,
                                     align = "right",
                                     fill = NA,
                                     na.rm = TRUE)),
         # Proportion of weeks below pH threshold
         pH_threshT24 = roll_mean(x = lag(pH_threshT,
                                          n = 1), 
                                  n = 24,
                                  align = "right",
                                  fill = NA,
                                  na.rm = TRUE),
         pH_threshT20 = roll_mean(x = lag(pH_threshT,
                                          n = 1), 
                                  n = 20,
                                  align = "right",
                                  fill = NA,
                                  na.rm = TRUE),
         pH_threshT16 = roll_mean(x = lag(pH_threshT,
                                          n = 1), 
                                  n = 16,
                                  align = "right",
                                  fill = NA,
                                  na.rm = TRUE),
         pH_threshT12 = roll_mean(x = lag(pH_threshT,
                                          n = 1), 
                                  n = 12,
                                  align = "right",
                                  fill = NA,
                                  na.rm = TRUE),
         pH_threshT8 = roll_mean(x = lag(pH_threshT,
                                         n = 1), 
                                 n = 8,
                                 align = "right",
                                 fill = NA,
                                 na.rm = TRUE),
         pH_threshT4 = roll_mean(x = lag(pH_threshT,
                                         n = 1), 
                                 n = 4,
                                 align = "right",
                                 fill = NA,
                                 na.rm = TRUE),
         # Accumulated Excess [H3O+] converted to pH scale
         # Add [H3O+] at pH 7.75 to each
         pH_threshU24 = -log10((10 ^ (-7.75)) +
                                 roll_sum(x = lag(pH_threshU, 
                                                  n = 1), 
                                          n = 24,
                                          align = "right",
                                          fill = NA,
                                          na.rm = TRUE)),
         pH_threshU20 = -log10((10 ^ (-7.75)) +
                                 roll_sum(x = lag(pH_threshU,
                                                  n = 1), 
                                          n = 20,
                                          align = "right",
                                          fill = NA,
                                          na.rm = TRUE)),
         pH_threshU16 = -log10((10 ^ (-7.75)) +
                                 roll_sum(x = lag(pH_threshU,
                                                  n = 1), 
                                          n = 16,
                                          align = "right",
                                          fill = NA,
                                          na.rm = TRUE)),
         pH_threshU12 = -log10((10 ^ (-7.75)) +
                                 roll_sum(x = lag(pH_threshU,
                                                  n = 1), 
                                          n = 12,
                                          align = "right",
                                          fill = NA,
                                          na.rm = TRUE)),
         pH_threshU8 = -log10((10 ^ (-7.75)) +
                                roll_sum(x = lag(pH_threshU,
                                                 n = 1), 
                                         n = 8,
                                         align = "right",
                                         fill = NA,
                                         na.rm = TRUE)),
         pH_threshU4 = -log10((10 ^ (-7.75)) +
                                roll_sum(x = lag(pH_threshU,
                                                 n = 1), 
                                         n = 4,
                                         align = "right",
                                         fill = NA,
                                         na.rm = TRUE))
  )

B10K_benv = banana

# remove temporary 
rm(banana)

# write csv -----
write.csv(B10K_benv, "BES_B10K_benv.csv")

# > 4.6: Match sampdate and t2 and attach t2, t, and ind -----
banana = neardate(id1 = dfm_QO1821$gridind,
                  id2 = B10K_btemp$gridind,
                  y1 = dfm_QO1821$sampdate2,
                  y2 = as.numeric(
                    as.POSIXct(B10K_benv$t2,
                               tz = "utc")))

dfm_QO1821$ind = B10K_benv$ind[banana]
dfm_QO1821$oceantime = B10K_benv$t2[banana]

# Last sampdate with valid ind
max(dfm_QO1821[!is.na(dfm_QO1821$ind),]$sampdate)
# First sampdate without valid ind (WOOO FAILS BECAUSE WE HAVE ALL THE DATES)
min(dfm_QO1821[is.na(dfm_QO1821$ind),]$sampdate)

# > 4.7: Attach B10K_env to dfm_QO1821 -----
dfm_QO1821 = left_join(dfm_QO1821,
                       B10K_benv)

# > 4.8: Match sampdate for df_potsum and t2 and attach t2, t, and ind -----
banana = neardate(id1 = df_potsum$gridind,
                  id2 = B10K_btemp$gridind,
                  y1 = df_potsum$sampdate2,
                  y2 = as.numeric(
                    as.POSIXct(B10K_benv$t2,
                               tz = "utc")))

df_potsum$ind = B10K_benv$ind[banana]
df_potsum$oceantime = B10K_benv$t2[banana]

# 5: Create Environmental Data for Graphing -----
# > 5.1: Find central index for each fishery season -----
# Note that this is not the true median, which would be weighted by the
# number of crab caught for each index's sampdate
banana = dfm_QO1821 %>%
  group_by(fishery) %>%
  summarise(median_ind = round(median(unique(ind))), # median index
            mean_ind = round(mean(unique(ind)))) # mean index, same as median
banana$oceantime = unique(B10K_benv[B10K_benv$ind %in% banana$median_ind,]$t2)

# > 5.2: Create df_envdat -----
# Begin with the basegrid
df_envdat = basegrid_masked

# Repeat the basegrid 4 times, once per fishery season
df_envdat = df_envdat[rep(seq_len(nrow(df_envdat)), 4),]

# Assign median_ind for each fishery season
df_envdat$ind = c(rep(banana$median_ind[1], length(basegrid_masked$gridind)),
                  rep(banana$median_ind[2], length(basegrid_masked$gridind)),
                  rep(banana$median_ind[3], length(basegrid_masked$gridind)),
                  rep(banana$median_ind[4], length(basegrid_masked$gridind)))

# Attach oceantime to df_envdat by matching ind/median_ind
df_envdat = left_join(df_envdat,
                      banana,
                      by = c("ind" = "median_ind"))
df_envdat$mean_ind = NULL

# Remove temporary
rm(banana)

# Omit missing data-lacking locations from df_envdat
# These should be grid cells with land
length(unique(df_envdat$gridind))
length(unique(B10K_benv$gridind))
df_envdat = df_envdat[df_envdat$gridind %in% unique(B10K_benv$gridind),]
length(unique(df_envdat$gridind))

# > 5.3: Attach environmental data to df_envdat -----
df_envdat = left_join(df_envdat,
                      B10K_benv)

# Add fishery index
df_envdat$fishery_ind = ifelse(df_envdat$fishery == "QO18",
                               1,
                               ifelse(df_envdat$fishery == "QO19",
                                      2,
                                      ifelse(df_envdat$fishery == "QO20",
                                             3,
                                             4)))

# 6: Center grid to the median and convert to km -----
range(df_envdat$gridX3571)
range(df_envdat$gridY3571)

plot(df_envdat$gridX3571,
     df_envdat$gridY3571)

gridX3571med = median(unique(df_envdat$gridX3571))
gridY3571med  = median(unique(df_envdat$gridY3571))

dfm_QO1821$gridXCkm = (dfm_QO1821$gridX3571 - gridX3571med) / 1000
dfm_QO1821$gridYCkm = (dfm_QO1821$gridY3571 - gridY3571med) / 1000

df_envdat$gridXCkm = (df_envdat$gridX3571 - gridX3571med) / 1000
df_envdat$gridYCkm = (df_envdat$gridY3571 - gridY3571med) / 1000

plot(df_envdat$gridXCkm,
     df_envdat$gridYCkm)

range(df_envdat$gridXCkm)
range(df_envdat$gridYCkm)

# 7: Attach df_potsum_agg to main data and environmental data frames ----
# Provides CPUE at the annual, total-fishery level and effort (as pots/fishery)
# at the gridind-ind level.
# Also attaching "ncrab" to account for local abundance

# aggregate df_potsum to ind and fishery year and gridind
df_potsum_agg = df_potsum %>%
  group_by(gridind, ind, fishery) %>%
  summarise(CPUE = unique(CPUE),
            CPUE_std = unique(CPUE_std),
            ncrab_ind = sum(catch),
            npots_ind = length(pot))

# (time) index-level CPUE
df_potsum_agg$CPUE_ind = df_potsum_agg$ncrab_ind / df_potsum_agg$npots_ind

# standardize ind-level CPUE
CPUE_ind_M = mean(df_potsum_agg$CPUE_ind)
CPUE_ind_SD = sd(df_potsum_agg$CPUE_ind)

df_potsum_agg$CPUE_ind_std = 
  (df_potsum_agg$CPUE_ind - CPUE_ind_M) /
  CPUE_ind_SD

# dfm_QO1821 uses the data to the gridind level
# note that pots are counted at the gridind level, but the CPUE calculations
# ultimately come from df_retcatch (fishery-level)
dfm_QO1821 = left_join(dfm_QO1821,
                       df_potsum_agg)

# 8: Export Both Main and Environmental Data Frames -----
# > 8.1: Main Data -----
# Write a csv
banana = dfm_QO1821
banana$coord4326 = NULL
banana$coord3571 = NULL

write.csv(banana,
          paste0("BES_dfm_QO1821_",
                 format(Sys.time(),
                        "%Y%m%d_%H%M"),
                 ".csv"))

# > 8.2: Environmental Data -----
# Write a csv
banana = df_envdat
banana$coord4326 = NULL
banana$coord3571 = NULL

write.csv(banana,
          paste0("BES_df_envdat_",
                 format(Sys.time(),
                        "%Y%m%d_%H%M"),
                 ".csv"))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
# END OF SCRIPT ------
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####