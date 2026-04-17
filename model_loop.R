# This script loops through all the relevant species and does the modeling
# steps:
# load occurence data created by the pinales_BC file
# load climate data
# loop through species:
  # 1. create biomod obj with BIOMOD_FormatingData
  # 2. use select07 from mecofun to select relevant bioclim variables
  # 3. reload into new biomod obj using BIOMOD_FormatingData
  # 4. create a number of models (GLM, MAXNET, RF) using the BIOMOD_Modeling function
  # 5. combine into a single ensemble model

# TODO: some species are missing ensemble models. prob bc not meeting threshold
# fix or toss?

# set working directory to directory of script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(terra)
library(mecofun)
library(dplyr)
library(biomod2)
library(tidyverse)

# study region
extent_environ <- ext(-139.06,-114.03,48.30,60.00)

# bioclim files
download_current <- rast(paste0("climate/wc2.1_2.5m/wc2.1_2.5m_bio_", 1:19, ".tif"))

current_climate_bc <- crop(download_current, extent_environ)

# load occurrence records
occurrence <- readRDS(file = "occurrence_iNat_100obs_y2000.Rds")

# species names
sp_names <- setdiff(names(occurrence), c("x", "y"))

counter <- 1

# looping through all the species
for(sp_name in sp_names){
  
  print(paste0("########## modeling species " , sp_name , " (" , counter , "/" ,
               length(sp_names) , ") ##############"))
  
  name_dot <- str_replace(sp_name, " ", ".")
  # skips the EM that were already created
  em_path <- paste0(name_dot, "/", name_dot, ".loop.ensemble.models.out")
  if(file.exists(em_path)){
    print("EM already computed.")
    counter <- counter + 1
    next
  }
  
  # creating the biomod2 obj
  myRespName <- sp_name
  myResp <- subset(occurrence, occurrence[[myRespName]] == 1)[,myRespName]
  myRespXY <- subset(occurrence, occurrence[[myRespName]] == 1)[,c("x", "y")]
  myExpl <- read.table("climate_grid_current.dat", header=T)
  #TODO. add test data
  biomod_data <- BIOMOD_FormatingData(resp.var = myResp,
                                      expl.var = rast(myExpl,
                                                      crs="+proj=latlon +ellps=WGS84"),
                                      resp.xy = myRespXY,
                                      resp.name = myRespName,
                                      PA.nb.rep = 2,
                                      PA.nb.absences = sum(myResp),
                                      PA.strategy = 'random',
                                      filter.raster = T)
  
  # selecting the highly correlated variables
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
  print(pred_sel)
  #only use the selected variables for modeling
  myExpl_sel <- myExpl[,c("x", "y", pred_sel)]
  
  # update the biomod2 data object with new variables
  biomod_data_sel <- BIOMOD_FormatingData(resp.var = myResp,
                                          expl.var = rast(myExpl_sel,
                                                          crs="+proj=latlon +ellps=WGS84"),
                                          resp.xy = myRespXY,
                                          resp.name = myRespName,
                                          PA.nb.rep = 2,
                                          PA.nb.absences = sum(myResp),
                                          PA.strategy = 'random',
                                          filter.raster = T)
  # create models
  BiomodModelOut <- BIOMOD_Modeling(bm.format = biomod_data_sel,
                                    modeling.id="loop",
                                    models =c("GLM", "RF",
                                              "MAXNET"),
                                    CV.strategy = 'random',
                                    CV.nb.rep = 2,
                                    CV.perc = 0.8,
                                    var.import = 3,
                                    metric.eval = c('TSS','ROC'))
  
  # combine to ensemble
  myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = BiomodModelOut,
                                        models.chosen = 'all',
                                        em.by = 'all',
                                        em.algo = c('EMca', 'EMwmean'),
                                        metric.select = 'TSS',
                                        metric.select.thresh = c(0.7),
                                        var.import = 3,
                                        metric.eval = c("TSS", "ROC"))
  
  print(get_evaluations(myBiomodEM))
  
  counter <- counter + 1
}
