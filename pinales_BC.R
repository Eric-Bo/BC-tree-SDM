
# set working directory to directory of script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#getwd()

# load packages
library(terra)
library(mecofun)
library(ggplot2)
library(gridExtra)
library(geodata)
library(bcmaps)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(biomod2)
library(tidyverse)
library(sf)

# define extent of BC
# source: https://gis.stackexchange.com/questions/154740/lat-and-long-extents-of-canadian-provinces
extent_environ <- ext(-139.06,-114.03,48.30,60.00)

#----------bioclimatic variables-------------------

# download current bioclimatic conditions from worldclim
#d <- worldclim_global(var="bio", res=2.5, path=getwd())

# if files are downloaded, they can be read directly
download_current <- rast(paste0("climate/wc2.1_2.5m/wc2.1_2.5m_bio_", 1:19, ".tif"))
# to study region
current_climate_bc <- crop(download_current, extent_environ)
# transform to dataframe
env_current <- as.data.frame(current_climate_bc, xy=TRUE)
names(env_current) <- c("x", "y", paste0("b", 1:19))
names(env_current)

# download climate bioclimatic data from worldclim for SSPRCP585 for ACCESS climate model
#r <- cmip6_world("ACCESS-CM2", "370", "2081-2100", var="bioc", path=getwd(), res=2.5)
# load worldclim data for future conditions if already downloaded
download_future <- rast(paste0("climate/wc2.1_2.5m/wc2.1_2.5m_bioc_ACCESS-CM2_ssp585_2081-2100.tif"))
# crop to study region
future_Climate_bc <- crop(download_future, extent_environ)
# make sure that current and future biolimate have same grid
future_Climate_bc <- resample(future_Climate_bc,current_climate_bc)
# transform to dataframe
env_future <- as.data.frame(future_Climate_bc, xy=TRUE)
names(env_future)  <- c("x", "y", paste0("b", 1:19))


# merge bioclimate data sets for current and future biolimate, then extract
# x and y coordinates; this makes sure that grids for current and future
# conditions are identical
coo_merge <- merge(env_current, env_future, by=c("x","y"))[,c("x","y")]
dim(coo_merge)

# merge bioclimate with coordinates such that current and future bioclimate
# grids are identical
env_current_merge <- merge(coo_merge, env_current, by=c("x","y"))
env_future_merge  <- merge(coo_merge, env_future, by=c("x","y"))
dim(env_current_merge)
dim(env_future_merge)

# save bioclimatic data
write.table(env_current_merge, "climate_grid_current.dat", col.names=T, row.names=F)
write.table(env_future_merge, "climate_grid_future585.dat", col.names=T, row.names=F)

#----------occurrence data------------

occurrence_raw <- data.frame(read.table("data/occurrence.txt", header=T, sep="\t", encoding="UTF-8", quote = "", na.strings=c("", "NA"), fill=T))
#needed columns: "datasetName", "year", "decimalLatitude", "decimalLongitude", "species"

occurrence <- occurrence_raw[, c("datasetName", "year", "decimalLatitude", "decimalLongitude", "species")]

#cleanup data 
#drop all rows with empty fields
occurrence <- na.omit(occurrence)
#drop all occurences before 2000
occurrence <- subset(occurrence, year>=2000)
#only use iNaturalist for now
occurrence <- subset(occurrence, datasetName == "iNaturalist research-grade observations")
#only use species with more than 100 observations
species_tab <- table(occurrence$species)

valid_species <- c()
for(i in 1:length(species_tab)){
  
  if(species_tab[i] > 100){
    valid_species[i] <- names(species_tab)[i]
  }
}

occurrence <- subset(occurrence, species %in% valid_species)
# 35 species meet the criteria
#length(valid_species)

# next few lines make the df into a wide frame thats easier to use
occ <- occurrence

# adds a presence column to make things easier
occ <- occ %>%
  mutate(presence = 1)

# drops points with multiple obs
occ <- occ %>%
  distinct(decimalLatitude, decimalLongitude, species, .keep_all = TRUE)

# widens the df
occ_wide <- occ  %>%
  pivot_wider(id_cols = c("decimalLatitude", "decimalLongitude"),
              names_from = species, 
              values_from = presence,
              values_fill = 0)

occurrence <- as.data.frame(occ_wide)

#also rename longitude and latitude
names(occurrence)[names(occurrence) == "decimalLatitude"] <- "y"
names(occurrence)[names(occurrence) == "decimalLongitude"] <- "x"

saveRDS(occurrence, file = "occurrence_iNat_100obs_y2000.Rds")

#---------test data----------
# TODO: I dont think this would be a smart idea
# use different dataset for the independent test data
occurrence_raw <- data.frame(read.table("data/occurrence.txt", header=T, sep="\t", encoding="UTF-8", quote = "", na.strings=c("", "NA"), fill=T))
#needed columns: "datasetName", "year", "decimalLatitude", "decimalLongitude", "species"

occ <- occurrence_raw[, c("datasetName", "year", "decimalLatitude", "decimalLongitude", "species")]

#cleanup data 
#drop all rows with empty fields
occ <- na.omit(occ)

occ_test <- subset(occ, datasetName == "Royal BC Museum - Botany Collection" | datasetName == "University of British Columbia Herbarium (UBC) - Vascular Plant Collection")
occ_test <- subset(occ_test, species %in% valid_species)
table(occ_test$species)

#--------build models------

# Modellierung mit Biomod2 Paket
# Erstmal testweise für nur eine Art
# dann in einer loop für alle

#loading data into a biomod object
occurrence <- readRDS(file = "occurrence_iNat_100obs_y2000.Rds")
myRespName <- "Abies amabilis"
myResp <- subset(occurrence, occurrence[[myRespName]] == 1)[,myRespName]
myRespXY <- subset(occurrence, occurrence[[myRespName]] == 1)[,c("x", "y")]
myExpl <- read.table("climate_grid_current.dat", header=T)
biomod_data <- BIOMOD_FormatingData(resp.var = myResp,
                                    expl.var = rast(myExpl,
                                                    crs="+proj=latlon +ellps=WGS84"),
                                    resp.xy = myRespXY,
                                    resp.name = myRespName,
                                    PA.nb.rep = 2,
                                    PA.nb.absences = 1000,
                                    PA.strategy = 'random',
                                    filter.raster = T)

#plot(biomod_data)


#before proceeding to forecasts the number of explaining variables should
#be minimized. using select07 from mecofun library

# get the pseudo absence data from the biomod obj
pa_data <- biomod_data@data.species
pa_data <- pa_data %>% replace(is.na(.), 0)
pa_coords <- biomod_data@coord

# Combine into a dataframe
pa_points <- data.frame(pa_coords, Presence_Absence = pa_data)

#get bioclim variables for the selected points
bio_sites <- data.frame(terra::extract(current_climate_bc,
                                cbind(pa_points$x,pa_points$y)))
names(bio_sites) <- paste0("b", 1:19)
var_sel <- select07(X=bio_sites[,paste0("b",1:19)], y=pa_points$Presence_Absence)

pred_sel <- var_sel$pred_sel
#only use the selected variables for modeling
myExpl_sel <- myExpl[,c("x", "y", pred_sel)]

biomod_data_sel <- BIOMOD_FormatingData(resp.var = myResp,
                                        expl.var = rast(myExpl_sel,
                                                        crs="+proj=latlon +ellps=WGS84"),
                                        resp.xy = myRespXY,
                                        resp.name = myRespName,
                                        PA.nb.rep = 2,
                                        PA.nb.absences = 1000,
                                        PA.strategy = 'random',
                                        filter.raster = T)

# models with advanced options:
# Define custom modeling options
#TODO. model options could be adjusted
# Modeling() needs the arg bm.options
# model_options <- bm_ModelingOptions(data.type = "binary",
#                                     models = c("GLM", "RF", "ANN", "MAXNET"),
#                                     bm.format = biomod_data_sel,
#                                     strategy = "default")

# Modifying GLM formula not necessary
# fllowing line shows that it already uses both linear and quadratic terms
#model_options@options$GLM.binary.stats.glm@args.values

# Modify ANN parameters (number of neurons and weight decay)
#model_options@options$ANN.binary.nnet.nnet@args.values$`_PA1_allRun`$size <- 99

# Check the updated options
#print(model_options)


#the data is then used for modeling

BiomodModelOut <- BIOMOD_Modeling(bm.format = biomod_data_sel,
                                  modeling.id="test",
                                 models =c("GLM", "RF",
                                           "MAXNET"),
                                 CV.strategy = 'random',
                                 CV.nb.rep = 2,
                                 CV.perc = 0.8,
                                 var.import = 3,
                                 metric.eval = c('TSS','ROC'),
                                 nb.cpu = 4)

# following two plots show evals of the different models
bm_PlotEvalMean(bm.out = BiomodModelOut)
# and the importance of explaining varaibles
#bm_PlotVarImpBoxplot(bm.out = BiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))

#projecting models with current conditions
# myBiomodProj <- BIOMOD_Projection(bm.mod = BiomodModelOut,
#                                   proj.name = 'Current',
#                                   new.env = rast(myExpl,
#                                                  crs="+proj=latlon +ellps=WGS84"),
#                                   models.chosen = 'all',
#                                   metric.binary = 'all',
#                                   metric.filter = 'all',
#                                   build.clamping.mask = TRUE)

#plot(myBiomodProj, str.grep = 'GLM')

#TODO. decide on a em.algo
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = BiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      em.algo = c('EMca', 'EMwmean'),
                                      metric.select = 'TSS',
                                      metric.select.thresh = c(0.7),
                                      var.import = 3,
                                      metric.eval = c("TSS", "ROC"),
                                      nb.cpu = 4 )
# the response curves of the different explaining variables
# myResponseCurve <- bm_PlotResponseCurves(bm.out = myBiomodEM,
#                                          models.chosen = get_built_models(myBiomodEM)[1],
#                                          fixed.var = 'mean')

#----------Forecasting-------------

#using the ensemble to predict the species distribution with current conditions
myBiomodEMProj_curr <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                 proj.name = 'CurrentEM',
                                                 new.env = rast(myExpl_sel,
                                                                crs="+proj=latlon +ellps=WGS84"),
                                                 models.chosen = 'all',
                                                 metric.binary = 'all',
                                                 metric.filter = 'all',
                                                 nb.cpu=4)
plot(myBiomodEMProj_curr)


#use the model to predict future distribution
myExpl_fut <- read.table("climate_grid_future.dat", header=T)
myExpl_sel_fut <- myExpl_fut[,c("x", "y", pred_sel)]
myBiomodEMProj_fut <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                 proj.name = 'FutureEM',
                                                 new.env = rast(myExpl_sel_fut,
                                                                crs="+proj=latlon +ellps=WGS84"),
                                                 models.chosen = 'all',
                                                 metric.binary = 'all',
                                                 metric.filter = 'all',
                                                 nb.cpu=4)
#show(myBiomodEMProj_fut)

#----------Forecasting from saved models-------------

# loading models from the saved models.out files

my_species <- str_replace(myRespName, " ", ".")
# use a seperate env bc load just dumps the obj into namespace
model_env <- new.env()
load(paste0(my_species, "/", my_species, ".test.ensemble.models.out"), envir = model_env)
loaded_mod <- model_env[[ls(model_env)[1]]]
# loaded em is just a vec of names
loaded_em <- BIOMOD_LoadModels(loaded_mod)

myBiomodEMProj_curr <- BIOMOD_EnsembleForecasting(bm.em = loaded_mod,
                                                  proj.name = 'CurrentEM',
                                                  new.env = rast(myExpl_sel,
                                                                 crs="+proj=latlon +ellps=WGS84"),
                                                  models.chosen = 'all',
                                                  metric.binary = 'all',
                                                  metric.filter = 'all',
                                                  nb.cpu=4)

myExpl_fut <- read.table("climate_grid_future.dat", header=T)
myExpl_sel_fut <- myExpl_fut[,c("x", "y", pred_sel)]
myBiomodEMProj_fut <- BIOMOD_EnsembleForecasting(bm.em = loaded_mod,
                                                 proj.name = 'FutureEM',
                                                 new.env = rast(myExpl_sel_fut,
                                                                crs="+proj=latlon +ellps=WGS84"),
                                                 models.chosen = 'all',
                                                 metric.binary = 'all',
                                                 metric.filter = 'all',
                                                 nb.cpu=4)
#------------plots (works)-----------
#TODO: add algo param so wie unten
# Load your ensemble forecast results
ensemble_current <- get_predictions(myBiomodEMProj_curr, full.name=myBiomodEMProj_curr@models.projected[2])
ensemble_future  <- get_predictions(myBiomodEMProj_fut, full.name=myBiomodEMProj_fut@models.projected[2])

# Convert to data frames for ggplot2
df_current <- as.data.frame(ensemble_current, xy = TRUE)
df_future  <- as.data.frame(ensemble_future, xy = TRUE)

# Rename columns for clarity
#myBiomodEMProj_curr@models.projected
colnames(df_current) <- c("Longitude", "Latitude", "Probability")
colnames(df_future)  <- c("Longitude", "Latitude", "Probability")

# Add scenario labels
df_current$Scenario <- "Current"
df_future$Scenario  <- "Future"

# Combine both datasets
df_all <- bind_rows(df_current, df_future)

# Get Canada boundaries
canada <- ne_states(country = "Canada", returnclass = "sf")

# Filter for British Columbia
bc_boundary <- canada %>% filter(name == "British Columbia")


ggplot() +
  geom_tile(data = df_all, aes(x = Longitude, y = Latitude, fill = Probability)) +  
  scale_fill_viridis_c(option = "C", name = "Suitability") +
  geom_sf(data = bc_boundary, fill = NA, color = "black", size = 5) +  # Add BC boundary
  facet_wrap(~Scenario) +
  coord_sf() +  # Ensure spatial alignment
  theme_minimal() +
  labs(title = "Species Distribution Model: British Columbia (Current vs Future)",
       x = "Longitude", y = "Latitude") +
  theme(legend.position = "bottom",
        text = element_text(size = 14))

# possible additional features:

# add the occurrence data:
# geom_point(data = myRespXY, aes(x=x, y=y), size=0.5, alpha = 0.7) +
# only works in this env bc i have myRespXY. would have to retrieve that later

#----------plots (discrete)-----------------

# Load your ensemble forecast results
# 1 = committee average
# 2 = weighted average
algo <- 1
ensemble_current <- get_predictions(myBiomodEMProj_curr, full.name=myBiomodEMProj_curr@models.projected[algo])
ensemble_future  <- get_predictions(myBiomodEMProj_fut, full.name=myBiomodEMProj_fut@models.projected[algo])


# Convert to data frames for ggplot2
df_current <- as.data.frame(ensemble_current, xy = TRUE)
df_future  <- as.data.frame(ensemble_future, xy = TRUE)

# Rename columns for clarity
colnames(df_current) <- c("Longitude", "Latitude", "Probability")
colnames(df_future)  <- c("Longitude", "Latitude", "Probability")

# Add scenario labels
df_current$Scenario <- "Current"
df_future$Scenario  <- "Future"

# Combine both datasets
df_all <- bind_rows(df_current, df_future)

# Get Canada boundaries
canada <- ne_states(country = "Canada", returnclass = "sf")

# Filter for British Columbia
bc_boundary <- canada %>% filter(name == "British Columbia")

# convert to discrete prediction
evals <- get_evaluations(loaded_mod)
cutoff <- evals[evals$full.name==myBiomodEMProj_curr@models.projected[algo] & 
                  evals$metric.eval=="TSS", ]$cutoff

df_all_bin <- df_all %>% mutate(Probability = ifelse(Probability>=cutoff,"Present","Absent"))

ggplot() +
  geom_tile(data = df_all_bin, aes(x = Longitude, y = Latitude, fill = Probability)) +  
  scale_fill_manual(values=c("Present"="green", "Absent"="grey"))+
  geom_sf(data = bc_boundary, fill = NA, color = "black", size = 5) +  # Add BC boundary
  facet_wrap(~Scenario) +
  coord_sf() +  # Ensure spatial alignment
  geom_point(data = myRespXY, aes(x=x, y=y, color="red"), size=0.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Species Distribution Model: British Columbia (Current vs Future)",
       x = "Longitude", y = "Latitude") +
  theme(legend.position = "bottom",
        text = element_text(size = 14))


#-----biogeoclimatic zones---------

# Create the map
# Load required libraries
library(bcmaps)
library(ggplot2)
library(sf)
library(dplyr)

# Load biogeoclimatic zones data using bec()
bc_zones <- bec()  # Load biogeoclimatic zones
bc_boundary <- bc_bound()  # Load BC boundary

bc_zones_simplified <- st_simplify(bc_zones, dTolerance = 1000)  # Increase tolerance to reduce resolution

bc_zones_grouped <- bc_zones_simplified %>%
  mutate(
    ZoneGroup = case_when(
      ZONE %in% c("CWH", "CWHvh1", "CWHvh2", "CWHws", "CDF", "CMA") ~ "Coastal",  # Coastal
      ZONE %in% c("IDF", "IDFxh", "IDFdk", "BAFA", "PP") ~ "Interior Dry",  # Interior Dry
      ZONE %in% c("BWBS") ~ "Boreal Forest",  # Boreal Forest
      ZONE %in% c("SBPS", "SBS") ~ "Interior Forest",  # Interior Forest
      ZONE %in% c("ESSF", "ICH", "MH", "MS") ~ "Mountain",  # Mountain
      TRUE ~ "Other"  # Default category for others
    )
  )

# Create the simplified map without borders between zones
p <- ggplot() +
  geom_sf(data = bc_zones_grouped, aes(fill = ZoneGroup), color = NA, size = 0.2) +  # Remove borders between zones
  geom_sf(data = bc_boundary, fill = NA, color = "black", size = 0.5) +  # Outline BC
  scale_fill_manual(values = c("Coastal" = "#1f78b4", 
                               "Interior Dry" = "#33a02c", 
                               "Boreal Forest" = "#e31a1c", 
                               "Interior Forest" = "#ff7f00", 
                               "Mountain" = "#6a3d9a", 
                               "Other" = "grey"),  # Custom color palette
                    name = "Biogeoclimatic Zone Group") +
  theme_minimal() +
  labs(title = "Biogeoclimatic Zones of British Columbia") +  # No caption
  theme(legend.position = "right",
        plot.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 10))
# Save the plot as a PDF file
ggsave("results/bc_bioclimatic_zones_grouped.pdf", plot = p, width = 10, height = 7, dpi = 300)

other_zones <- bc_zones_grouped %>% filter(ZoneGroup == "Other")
print(unique(other_zones$ZONE))
print(unique(bc_zones_simplified$ZONE))
