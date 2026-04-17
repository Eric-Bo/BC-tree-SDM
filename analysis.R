# this script does all the analysis for the different pinales species


# set working directory to directory of script
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

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


# study region
extent_environ <- ext(-139.06,-114.03,48.30,60.00)

# bioclim files
download_current <- rast(paste0("climate/wc2.1_2.5m/wc2.1_2.5m_bio_", 1:19, ".tif"))

current_climate_bc <- crop(download_current, extent_environ)

# load occurrence records
occurrence <- readRDS(file = "occurrence_iNat_100obs_y2000.Rds")

# species names
sp_names <- setdiff(names(occurrence), c("x", "y"))

# Get Canada boundaries
canada <- ne_states(country = "Canada", returnclass = "sf")

# Filter for British Columbia
bc_boundary <- canada %>% filter(name == "British Columbia")

# create df that will contain all projections
myExpl <- read.table("climate_grid_current.dat", header=T)
projections <- data.frame()

# collect the evaluations of the models in a df
eval_df <- data.frame()

# collect the expl variables
expl_vars <- list()

counter <- 1

# looping through all the species
for(sp_name in sp_names){
  
  # check if the ensemble model exists
  name_dot <- str_replace(sp_name, " ", ".")
  # skips the EM that were already created
  em_path <- paste0(name_dot, "/", name_dot, ".loop.models.out")
  if(!file.exists(em_path)){
    counter <- counter + 1
    next
  }
  
  # load ensemble models
  my_species <- str_replace(sp_name, " ", ".")
  # use a seperate env bc load just dumps the obj into namespace
  model_env <- new.env()
  load(paste0(my_species, "/", my_species, ".loop.models.out"), envir = model_env)
  loaded_mod <- model_env[[ls(model_env)[1]]]
  
  
  break
  
  #--------evaluations-----
  
  
  
  # evals <- get_evaluations(loaded_mod)
  # evals <- evals[,c("algo", "metric.eval", "calibration")]
  # evals$species <- sp_name
  # if(counter == 1){
  #   eval_df <- evals
  #   counter <- counter + 1
  # } else {
  #   eval_df <- rbind(eval_df, evals)
  # }
  # next
  
  # get the eval metrics and corresponding plots
  # pdf(paste0("results/evals/", name_dot, ".pdf"))
  # bm_PlotEvalMean(loaded_mod)
  # dev.off()
  
  #------projecting----------
  
  # # get the explaining variables used in this model
  # expl_var <- names(table(get_variables_importance(loaded_mod)$expl.var))
  # 
  # # myExpl_sel <- myExpl[,c("x", "y", expl_var)]
  # # 
  # # myBiomodEMProj_curr <- BIOMOD_EnsembleForecasting(bm.em = loaded_mod,
  # #                                                   proj.name = 'CurrentEM',
  # #                                                   new.env = rast(myExpl_sel,
  # #                                                                  crs="+proj=latlon +ellps=WGS84"),
  # #                                                   models.chosen = 'all',
  # #                                                   metric.binary = 'all',
  # #                                                   metric.filter = 'all')
  # # 
  # myExpl_fut <- read.table("climate_grid_future585.dat", header=T)
  # myExpl_sel_fut <- myExpl_fut[,c("x", "y", expl_var)]
  # myBiomodEMProj_fut <- BIOMOD_EnsembleForecasting(bm.em = loaded_mod,
  #                                                  proj.name = 'Future585EM',
  #                                                  new.env = rast(myExpl_sel_fut,
  #                                                                 crs="+proj=latlon +ellps=WGS84"),
  #                                                  models.chosen = 'all',
  #                                                  metric.binary = 'all',
  #                                                  metric.filter = 'all')
  # 
  #-----plotting---------
  
  # # recovering projections from file
  # # use a separate env bc load just dumps the obj into namespace
  # cur_proj_env <- new.env()
  # load(paste0(my_species, "/proj_CurrentEM/", my_species, ".CurrentEM.ensemble.projection.out"), envir = cur_proj_env)
  # myBiomodEMProj_cur <- cur_proj_env[[ls(cur_proj_env)[1]]]
  # 
  # # fut_proj_env <- new.env()
  # # load(paste0(my_species, "/proj_FutureEM/", my_species, ".FutureEM.ensemble.projection.out"), envir = fut_proj_env)
  # # myBiomodEMProj_fut <- fut_proj_env[[ls(fut_proj_env)[1]]]
  # 
  # # discrete:
  # # 1 = committee average
  # # 2 = weighted average
  # algo <- 1
  # ensemble_current <- get_predictions(myBiomodEMProj_cur, full.name=myBiomodEMProj_cur@models.projected[algo])
  # ensemble_future  <- get_predictions(myBiomodEMProj_fut, full.name=myBiomodEMProj_fut@models.projected[algo])
  # 
  # 
  # # Convert to data frames for ggplot2
  # df_current <- as.data.frame(ensemble_current, xy = TRUE)
  # df_future  <- as.data.frame(ensemble_future, xy = TRUE)
  # 
  # # Rename columns for clarity
  # colnames(df_current) <- c("Longitude", "Latitude", "Presence")
  # colnames(df_future)  <- c("Longitude", "Latitude", "Presence")
  # 
  # # Add scenario labels
  # df_current$Scenario <- "Current"
  # df_future$Scenario  <- "Future (SSP5-8.5, 2081-2100)"
  # 
  # # Combine both datasets
  # df_all <- bind_rows(df_current, df_future)
  # 
  # # convert to discrete prediction
  # evals <- get_evaluations(loaded_mod)
  # cutoff <- evals[evals$full.name==myBiomodEMProj_cur@models.projected[algo] & 
  #                   evals$metric.eval=="TSS", ]$cutoff
  # 
  # df_all_bin <- df_all %>% mutate(Presence = ifelse(Presence>=cutoff,"Present","Absent"))
  # 
  # # combine projections
  # if(counter == 1){
  #   projections <- df_all_bin
  #   colnames(projections)[3] <- sp_name
  # } else {
  #   projections[sp_name] <- df_all_bin$Presence
  # }
  # counter <- counter + 1
  
  # -----------discrete plots for single species------------
  
  # # get the occ x and y for the species
  # occ_sp <- subset(occurrence, occurrence[[sp_name]] == 1)
  # occ_xy <- occ_sp[,c("x", "y")]
  # 
  # 
  # # plot for discrete distribution comparisons for current and future conditions
  # p <- ggplot() +
  #     geom_tile(data = df_all_bin, aes(x = Longitude, y = Latitude, fill = Presence), alpha = 0.8) +
  #     scale_fill_manual(values=c("Present"="darkgreen", "Absent"="grey"))+
  #     geom_sf(data = bc_boundary, fill = NA, color = "black", size = 2) +  # Add BC boundary
  #     facet_wrap(~Scenario) +
  #     coord_sf() +  # Ensure spatial alignment
  #     geom_point(data = occ_xy %>% mutate(Scenario = "Current"), aes(x=x, y=y, color = "Observation"), size=0.1, alpha = 0.3) +
  #     scale_color_manual(values = c("Observation" = "red")) +  
  #   guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +  
  #   theme_minimal() +
  #     labs(title = paste0("Species Distribution Model for ", sp_name ," in B.C."),
  #          x = "Longitude", y = "Latitude") +
  #     theme(legend.position = "bottom",
  #           text = element_text(size = 14),
  #           panel.spacing = unit(1, "lines"))
  # ggsave(paste0("results/plots/discrete_ca_", name_dot, ".pdf"),
  #        width = 23,
  #        height = 12,
  #        units = "cm")
  
  
  
  
}


#-----total combined species count--------

#projections <- readRDS("projections.rds")
# to count the number of occurences in each cell a numeric df is created
projections_num <- projections %>%
  mutate(across(everything(), ~ ifelse(. == "Present", 1,
                                       ifelse(. == "Absent", 0, .)))) %>%
  mutate(across(where(~ all(suppressWarnings(!is.na(as.numeric(.))))), as.numeric))

projections_num$total <- rowSums(projections_num[,c(-1,-2,-4)])
projections_totals <- projections_num[,c("Longitude", "Latitude", "Scenario", "total")]

projections_totals <- projections_totals %>% 
  pivot_wider( names_from = Scenario, values_from = total)

projections_totals$change <- projections_totals$`Future (SSP5-8.5, 2081-2100)` - projections_totals$Current

saveRDS(projections_totals, "projections_585.rds")

p <- ggplot() +
    geom_tile(data = projections_totals, aes(x = Longitude, y = Latitude, fill = change)) +
    scale_fill_gradient2(low = "darkred", mid = "white", high = "darkgreen", midpoint=0) +
    geom_sf(data = bc_boundary, fill = NA, color = "black", size = 2) +
    coord_sf() +  # Ensure spatial alignment
    #geom_point(data = myRespXY, aes(x=x, y=y, color="red"), size=0.5, alpha = 0.7) +
    theme_minimal() +
  ggtitle("Species richness change for 17 Pinales species in B.C.\nfor SSP3-7.0 in 2081-2100")+
    labs(x = "Longitude", y = "Latitude") +
    theme(legend.position = "bottom",
          text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))

ggsave("results/combined_total_SSP3.pdf",
       width = 17,
       height = 12,
       units = "cm", dpi=300)

# combine all scenarios into one plot
proj1 <- readRDS("projections_126.rds")
proj2 <- readRDS("projections_245.rds")
proj3 <- readRDS("projections_370.rds")
proj4 <- readRDS("projections_585.rds")

proj1 <- proj1[,c("Longitude", "Latitude", "change")]
proj2 <- proj2[,c("Longitude", "Latitude", "change")]
proj3 <- proj3[,c("Longitude", "Latitude", "change")]
proj4 <- proj4[,c("Longitude", "Latitude", "change")]

proj1$SSP <- "SSP1-2.6"
proj2$SSP <- "SSP2-4.5"
proj3$SSP <- "SSP3-7.0"
proj4$SSP <- "SSP5-8.5"


proj_all <- rbind(proj1, proj2, proj3, proj4)

p <- ggplot(proj_all) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = change)) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "darkgreen", midpoint=0) +
  geom_sf(data = bc_boundary, fill = NA, color = "black", size = 2) +
  coord_sf() +  # Ensure spatial alignment
  theme_minimal() +
  ggtitle("Species richness change for 17 Pinales species in B.C. (2081-2100)") +
  labs(x = "Longitude", y = "Latitude", fill = "Change in Richness") +
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        panel.spacing.x = unit(2, "cm"),
        strip.text = element_text(size = 14)) +
  facet_wrap(~SSP, ncol = 2)

# Compute mean change for each SSP
summary_stats <- proj_all %>%
  group_by(SSP) %>%
  summarize(
    mean_change = mean(change, na.rm = TRUE),  
    max_lat = max(Latitude)  
  )

# Adjust y-position to be slightly above the plots
summary_stats$label_y <- summary_stats$max_lat + 0.5  # Adjust this value as needed

# Create the text label
summary_stats$label <- paste0("Mean change: ", round(summary_stats$mean_change, 2))

# Add text above each subplot
p <- p + 
  geom_text(data = summary_stats, aes(x = mean(range(proj_all$Longitude)), y = label_y, label = label),
            size = 5, hjust = 0.5, vjust = 0, color = "black")


ggsave("results/species_richness_change.pdf",
       width = 25,
       height = 24,
       units = "cm", dpi=300)

# boxplot doesnt look good tho
q <- ggplot(proj_all, aes(x = SSP, y = change, fill = SSP)) +
  geom_boxplot(outlier.shape = NA) +  # Avoid too many outliers cluttering the plot
  geom_jitter(alpha = 0.3, width = 0.2) +  # Adds data points for clarity
  theme_minimal() +
  scale_fill_manual(values = c("SSP3-7.0" = "red", "SSP5-8.5" = "blue", "SSP1-2.6" = "green", "SSP2-4.5" = "orange")) +
  ggtitle("Distribution of Species Richness Change by SSP") +
  labs(x = "SSP Scenario", y = "Species Richness Change") +
  theme(legend.position = "none", text = element_text(size = 14), plot.title = element_text(hjust = 0.5))

ggsave("results/species_richness_boxplot.pdf",
       width = 25,
       height = 14,
       units = "cm", dpi=300)


#------plotting evals----
eval_wide <- eval_df %>% 
  pivot_wider(names_from = metric.eval, values_from = calibration)

ggplot(eval_wide, aes(x = TSS, y = ROC, color = species, shape = algo)) +
  geom_point(size = 4, stroke = 1.2) +  # Larger points with an outline
  theme_minimal() +  
  labs(
    x = "True Skill Statistic (TSS)",
    y = "Area under the Curve (AUC)",
    title = "Evaluation of Model Performance: TSS vs. AUC",
    shape = "ensemble selection algorithm"
  ) +
  scale_shape_manual(values = c(16, 17), labels=c("Committee averaging", "weighted mean")) +  # Different shapes for EMca and EMwmean
  scale_color_viridis_d() +  # Color palette for 17 species
  theme(
    plot.title = element_text(face = "bold", size = 16),
    text = element_text(size = 14),
    legend.position = "right",
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(ncol = 2))

ggsave(paste0("results/evaluation_metrics.pdf"),
       width = 23,
       height = 12,
       units = "cm")


table(occurrence)
species_counts <- colSums(occurrence[, -c(1,2)])
print(species_counts)


#--------plotting bioclim changes---------

clim_cur <- read.table("climate_grid_current.dat", header=T)
clim_fut <- read.table("climate_grid_future585.dat", header=T)

# Define a custom function to calculate the percentage change
calculate_change <- function(current, future, var) {
  # For b1, calculate the absolute difference
  if (var == "b1") {
    return(future - current)  # Absolute difference for b1
  } else {
    # For other variables, calculate the percentage change
    current[current == 0] <- NA  # Avoid division by zero
    perc_change <- (future - current) / current * 100
    perc_change[is.infinite(perc_change)] <- NA  # Handle infinite values
    return(perc_change)
  }
}

# Identify the four most important variables
top_vars <- c("b15", "b3", "b18", "b1")  # Adjust as needed

# Compute change (absolute difference or percentage) between future and current climate for selected variables
climate_change <- clim_cur %>%
  select(x, y, all_of(top_vars)) %>%
  rename_with(~ paste0(.x, "_cur"), -c(x, y)) %>%
  left_join(
    clim_fut %>%
      select(x, y, all_of(top_vars)) %>%
      rename_with(~ paste0(.x, "_fut"), -c(x, y)),
    by = c("x", "y")
  )

# Calculate the change manually for each variable using the custom function
for (var in top_vars) {
  curr_col <- paste0(var, "_cur")
  fut_col <- paste0(var, "_fut")
  
  climate_change[[paste0("change_", var)]] <- mapply(calculate_change,
                                                     climate_change[[curr_col]],
                                                     climate_change[[fut_col]],
                                                     MoreArgs = list(var))
}

# Check for NA or Inf values
summary(climate_change)

# Clip extreme values if needed
climate_change <- climate_change %>%
  mutate(across(starts_with("change"),
                ~ pmin(pmax(., -1000), 1000),
                .names = "clipped_{.col}"))

# Pivot to long format for plotting
climate_change_long <- climate_change %>%
  select(x, y, starts_with("clipped_change")) %>%
  pivot_longer(cols = starts_with("clipped_change"), 
               names_to = "Variable", 
               values_to = "change") %>%
  mutate(Variable = gsub("clipped_change_", "BIO", Variable))  # Format variable names

# Generate and save individual plots for each variable
# Generate and save individual plots for each variable
# Generate and save individual plots for each variable
for (var in top_vars) {
  # Filter data for the current variable
  plot_data <- subset(climate_change_long, Variable == paste0("BIO", var))
  
  if (var == "b1") {
    # Custom color scale for b1 (light red to dark red for increases)
    p <- ggplot(plot_data) +
      geom_tile(aes(x = x, y = y, fill = change)) +
      scale_fill_gradient(low = "white", high = "darkred", 
                          limits = c(min(plot_data$change, na.rm = TRUE), max(plot_data$change, na.rm = TRUE))) +
      geom_sf(data = bc_boundary, fill = NA, color = "black", size = 2) +
      coord_sf() +
      theme_minimal() +
      ggtitle("Absolute Change in b1") +
      labs(x = "Longitude", y = "Latitude", fill = "Change") +
      theme(legend.position = "bottom",
            legend.key.size = unit(1.5, "cm"),  # Adjust the length of legend color bar
            text = element_text(size = 16),
            plot.title = element_text(hjust = 0.5),
            panel.spacing.x = unit(2, "cm"),
            strip.text = element_text(size = 14))
  } else {
    # For other variables, use a broader scale (percentage change)
    p <- ggplot(plot_data) +
      geom_tile(aes(x = x, y = y, fill = change)) +
      scale_fill_gradient2(low = "darkgreen", mid = "white", high = "darkred", midpoint = 0, 
                           limits = c(min(plot_data$change, na.rm = TRUE), max(plot_data$change, na.rm = TRUE))) +
      geom_sf(data = bc_boundary, fill = NA, color = "black", size = 2) +
      coord_sf() +
      theme_minimal() +
      ggtitle(paste("Percentage Change in", var)) +
      labs(x = "Longitude", y = "Latitude", fill = "% Change") +
      theme(legend.position = "bottom",
            legend.key.size = unit(1.5, "cm"),  # Adjust the length of legend color bar
            text = element_text(size = 16),
            plot.title = element_text(hjust = 0.5),
            panel.spacing.x = unit(2, "cm"),
            strip.text = element_text(size = 14))
  }
  
  # Save the plot as a file (adjust file format and path as needed)
  ggsave(paste0("results/change_", var, ".pdf"), plot = p, width = 8, height = 6)
}









