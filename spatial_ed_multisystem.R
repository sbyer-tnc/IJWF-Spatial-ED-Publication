# Multi-system (aka Real Estate) Ecological Departure
# Same as single-system spatial ED BUT multiplied by proportion of focal system present in the window

# Requires another moving window, but not a total re-calculation of spatial ED
# Use outputs from single-system ED script
# spatial_ed_singlesystem.R

#-----------------------------------------------------------------
#-----------------------------------------------------------------
# Libraries

library(rgdal)
library(sp)
library(raster)
library(dplyr)
library(rgeos)
library(spdep)
library(colorRamps)

#-----------------------------------------------------------------
#-----------------------------------------------------------------
# Local window version

# User defines window radius (should be same as in single-system)
# Define raster cell resolution
# Outputs matrix of systems found within window
simpleWindow <- function(raster, cell_id, cell_res, pixel_radius){
  rbuf <- setValues(raster, NA) # Set all values to NA
  rbuf[cell_id] <- 1 # Set our pixel value to 1
  buff <- buffer(rbuf, width = cell_res*pixel_radius) # From that pixel, buffer by our pixel radius (multiply by 60 for distance)
  win <- rasterToPolygons(buff, dissolve=TRUE) # Convert the bufer to a polygon
  winx <- crop(mask(raster, buff), win)
  winm <- as.matrix(winx)
  return(winm)
}

# Real estate function
######################
# REQUIRES SYSTEM-ONLY RASTER!!!
######################
# Removes edge effects
realEstate <- function(sys_code){
  r_path <- "TBD\\spatialed_mean_"
  r_ext <- ".tif"
  r_full_path <- paste(paste(r_path, as.character(sys_code), sep=""), r_ext, sep="")
  print(r_full_path)
  ed <- raster(r_full_path)
  plot(ed)
  df_re_values <- data.frame("cell_id" = as.numeric(), "og_ed" = as.numeric(), 
                             "prop_sys" = as.numeric(),  "re_ed" = as.numeric())
  docells <- which(!is.na(ed[]))
  print(length(docells))
  for(i in docells){
    print(i)
    sys_raster <- raster("TBD_ReferenceOnly_SYS.tif")
    w <- simpleWindow(raster = sys_raster, cell_id = i, pixel_radius = 15)
    # Need to account for number of white-space (NAs) pixels surrounding the window
    # Using the SYS raster removes need for subtracting 253 whitespace pixels
    prop_sys <- length(which(w == sys_code))/sum(!is.na(w))
    re_ed <- ed[i]*prop_sys
    re_tbl <- data.frame("cell_id" = i, "og_ed" = ed[i], 
                         "prop_sys" = prop_sys,  "re_ed" = re_ed)
    df_re_values <- rbind(df_re_values, re_tbl)
    # Assign new value to pixel and move on
    ed[i] <- re_ed
  }
  plot(ed)
  names(ed) <- paste("realestate_", as.character(sys_code), sep="")
  outlist <- list("tbl" = df_re_values, "raster" = ed)
  return(outlist)
}


#-----------------------------------------------------------------
#-----------------------------------------------------------------
# Apply to systems

# 10190 - PJ woodland
re_sys <- realEstate(sys_code = 10190)
writeRaster(re_sys$raster, "TBD\\multi_spatial_ed_10190.tif",
            overwrite=TRUE)
write.csv(re_sys$tbl, "TBD\\multi_spatial_ed_10190.csv")


# More systems...
# Copy above snippet, change 5-digit SYS code to get table for your desired system

#-----------------------------------------------------------------
#-----------------------------------------------------------------
# Use single-system spatial ed weights and values to get multi-system variance

propRaster <- function(sys_code){
  # Real estate raster
  r_path <- "path-to-raster\\realestate_ed_"
  r_ext <- ".tif"
  r_full_path <- paste(paste(r_path, as.character(sys_code), sep=""), r_ext, sep="")
  print(r_full_path)
  fill_raster <- raster(r_full_path)
  # Table with proportions at each pixel
  t_path <- "path-to-table\\realestate_ed_"
  t_ext <- ".csv"
  t_full_path <- paste(paste(t_path, as.character(sys_code), sep=""), t_ext, sep="")
  print(t_full_path)
  ed_tbl <- read.csv(t_full_path)
  # Real estate ED table with proportions already calculated at pixels
  # Assign cell proportion values to raster cells by ID
  dum_tbl <- data.frame(Dum_ID = as.numeric(seq(1:ncell(fill_raster))))
  new_tbl <- merge(dum_tbl, ed_tbl, by.x = 'Dum_ID', by.y = 'cell_id', all.x=TRUE)
  r <- fill_raster
  r[] <- new_tbl$prop_sys
  plot(r)
  return(r)
}

run_ED <- function(in_df, names, ref_raster, out_name){
  sub_df <- in_df[in_df$NRV_ID %in% names,]
  group_prop <- sub_df %>% group_by(SYSXCLA, Pixel_ID) %>% summarise(MinProp = min(`Weighted Proportion`, na.rm=TRUE)) %>%
    group_by(Pixel_ID) %>% summarise(SumProp = sum(MinProp, na.rm=TRUE))
  group_prop$ED_metric <- 1 - group_prop$SumProp
  # Create raster
  ed_raster <- ref_raster
  names(ed_raster) <- out_name
  # Account for NAs then assign ED values to pixels
  dum_tbl <- data.frame(Dum_ID = as.numeric(seq(1:ncell(ed_raster))))
  new_tbl <- merge(dum_tbl, group_prop, by.x = 'Dum_ID', by.y = 'Pixel_ID', all.x=TRUE)
  ed_raster[] <- new_tbl$ED_metric
  return(ed_raster)
}

# Theoretically, could use this to write mean real estate ED too...
realEstateVar <- function(sys_code, prop_sys_raster){
  # Raster to populate with variance values
  r_path <- "K:\\GIS3\\Projects\\Bretzlaff\\Geodata\\Raster\\tjr_ed\\spatialed_mean_"
  r_ext <- ".tif"
  r_full_path <- paste(paste(r_path, as.character(sys_code), sep=""), r_ext, sep="")
  print(r_full_path)
  fill_raster <- raster(r_full_path)
  # Get window proportion value - need to multiple all NRV rasters by this proportion value
  # df_re_values <- data.frame("cell_id" = as.numeric(), "prop_sys" = as.numeric())
  # docells <- which(!is.na(fill_raster[]))
  # print(length(docells))
  # for(i in docells){
  #   print(i)
  #   w <- simpleWindow(raster = fill_raster, cell_id = i, pixel_radius = 15)
  #   prop_sys <- sum(!is.na(w))/length(w)
  #   fill_raster[i] <- prop_sys
  # }
  # Table with weights at each pixel calculated for all NRVs
  t_path <- "path-to-table-with-weights\\df_weights_nrvs_"
  t_ext <- ".csv"
  t_full_path <- paste(paste(t_path, as.character(sys_code), sep=""), t_ext, sep="")
  print(t_full_path)
  ed_tbl <- read.csv(t_full_path)
  ed_tbl <- ed_tbl %>% dplyr::select(-c(X))
  colnames(ed_tbl) <- c("NRV_ID", "Pixel_ID", "SYSXCLA", "Number of Pixels", "Raw Proportion", "Weighted Proportion", "Weighted Count")
  it_names <- unique(ed_tbl$NRV_ID)[-1] # Remove initial conditions results from NRV list
  nrvi_stack <- stack()
  for(i in it_names){
    print(i)
    nrvi <- run_ED(in_df = ed_tbl, names = c("TJR_SYSXCLA_60m_50415", i),
                   ref_raster = fill_raster, out_name = i)
    nrvi_stack<- stack(nrvi_stack, nrvi)
  }
  return(nrvi_stack*prop_sys_raster)
}

#-----------------------------------------------------------------
#-----------------------------------------------------------------
# Apply to systems

# PJ woodland
rx <- propRaster(sys_code = 10190)
plot(rx)
r_stack <- realEstateVar(sys_code = 10190, prop_sys_raster = rx) 
r_var <- calc(r_stack, fun=var)
plot(r_var)
writeRaster(r_var, "TBD\\realestate_variance_10190.tif")

# End multi-system ED