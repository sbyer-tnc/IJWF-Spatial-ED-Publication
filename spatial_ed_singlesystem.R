#-----------------------------------------------------------------
# Create Spatial ED 
# Compare Year-0 TJR with multiple NRV runs
# ED value for each pixel/NRV run combo
# Get the mean ED value for each pixel

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
# Define functions

# Get NRV - read in table of NRV proportions and select the system you want
# outputs a table in correct format for calculating NRV
get_nrv <- function(nrv_tbl, sys_code){
  tbl <- read.csv(nrv_tbl)
  sub_tbl <- tbl[tbl$SYS_CODE == sys_code,]
  # Collapse uncharacteristic classes to one code
  sub_tbl[sub_tbl$CLA_CODE > 99,]$SYSXCLA <- as.numeric(as.character(paste(sys_code, "999", sep = "")))
  sub_tbl <- sub_tbl %>% group_by(SYSXCLA) %>% summarise(ExpectProp = sum(NRV))
  sub_tbl <- as.data.frame(sub_tbl)
  # Force uncharacteristic code NRV to be slightly > 0
  sub_tbl[sub_tbl$SYSXCLA == as.numeric(as.character(paste(sys_code, "999", sep = ""))),]$ExpectProp <- 0.00001
  sub_tbl$SYSXCLA <- as.character(sub_tbl$SYSXCLA)
  return(sub_tbl)
} 


# Input system code and raster, function will clear all pixels that are not of that system
select_sys <- function(sys_code, in_raster){
  sys_code_char <- as.character(sys_code)
  sys_code_low <- as.numeric(paste(sys_code_char, "000", sep=""))
  sys_code_high <- as.numeric(paste(sys_code_char, "999", sep=""))
  out_raster <- in_raster
  out_raster[out_raster > sys_code_high] <- NA
  out_raster[out_raster < sys_code_low] <- NA
  return(out_raster)
}

# Reclassify uncharacteristic classes in single-system raster
unchar_reclass <- function(sys_code, in_raster){
  sys_code_char <- as.character(sys_code)
  sys_code_low <- as.numeric(paste(sys_code_char, "000", sep=""))
  sys_code_high <- as.numeric(paste(sys_code_char, "999", sep=""))
  raster_reclass <- reclassify(in_raster, c(sys_code_low, sys_code_high, sys_code_high),
                               include.lowest=TRUE)
  return(raster_reclass)
}


# Create a circular window around the pixel and return a table of SYSXCLA values in that window and a matrix of those values in the window
# Also returns inputs for creating a euclidean distance raster/matrix
circleWindow <- function(raster, cell_id, cell_res, pixel_radius){
  rbuf <- setValues(raster, NA) # Set all values to NA
  rbuf[cell_id] <- 1 # Set our pixel value to 1
  buff <- buffer(rbuf, width = cell_res*pixel_radius) # From that pixel, buffer by our pixel radius (multiply by 60 for distance)
  win <- rasterToPolygons(buff, dissolve=TRUE) # Convert the bufer to a polygon
  win_sysxcla <- data.frame(table(extract(raster, win), useNA="ifany")) # Use buffer polygon to get sysxcla values in window
  colnames(win_sysxcla) <- c("SYSXCLA", "Freq")
  win_sysxcla$SYSXCLA <- as.character(win_sysxcla$SYSXCLA)
  win_sysxcla$Proportion <- win_sysxcla$Freq/sum(win_sysxcla$Freq) # Create a frequency table of the raw proportions in the window
  winx <- crop(mask(raster, buff), win)
  winm <- as.matrix(winx)
  out_list <- list("tbl" = win_sysxcla, "Window_poly" = win, "Pixel" = rbuf, "Window" = buff, "SYSXCLA_matrix" = winm)
  return(out_list)
}


# Calculate distance from the pixel within the window
getDist <- function(pixel_layer, window, window_polygon){
  wd <- mask(distance(pixel_layer), window) # distance() defaults to Euclidean distance
  # If we want a distance decay, need a decay function (Wd = e^(decayrate*d)) to canned distance function
  # OPTIONAL - add decay to distance raster
  wdx <- 1.0171*exp(-0.004*wd) # Taken from CTZ_FireStarts.py for decaying distance from Roads
  sumwdx <- sum(wdx[], na.rm=TRUE) # Use sum of all distance values to normalize weights
  wraster <- setValues(wdx, wdx[]/sumwdx) # Normalize weights in raster
  wraster <- crop(mask(wraster, window), window_polygon) # Crop/mask raster to window
  # Optional, but do it anyway - convert distance values to 0 to 1
  wraster <- wraster*(1/max(wraster[], na.rm=TRUE))
  wmatrix <- as.matrix(wraster) # Convert raster to matrix
  #outlist <- list("decay_raster" = wraster, "decay_matrix" = wmatrix) # Return raster and matrix
  outlist <- list("decay_matrix" = wmatrix) # Return only matrix
  return(outlist)
}



# Function to create an ecological distance raster
ecoRaster <- function(sys_raster, sys_code, lut){
  #sys_raster[pixel] <- 1
  # Reclassify SYS raster
  lut_select <- lut[lut$SYS_CODE1 == sys_code,] # Subset of full lookup table
  lut_rcl <- cbind(lut_select$SYS_CODE2, lut_select$Value) # 2-column matrix, from-to
  r_reclass <- reclassify(sys_raster, lut_rcl)
  #eco_raster <- abs(r_reclass) # Take the absolute value of the ecological distance raster's values
  eraster <- (r_reclass + 1)/2
  return(eraster)
}

ecoDist <- function(eco_raster, window, window_polygon){
  r_window <- crop(mask(eco_raster, window), window_polygon)
  ematrix <- as.matrix(r_window) # Convert raster to matrix
  #outlist <- list("eco_raster" = r_window, "eco_matrix" = ematrix) # Return raster and matrix
  outraster <- list("eco_matrix" = ematrix) # Return only matrix
  return(outraster)
}


# Function to create BINARY ecological distance raster
binary_ecoRaster <- function(sys_raster, sys_code, lut){
  #sys_raster[pixel] <- 1
  # Reclassify SYS raster
  lut_select <- lut[lut$SYS_CODE1 == sys_code,] # Subset of full lookup table
  lut_rcl <- cbind(lut_select$SYS_CODE2, lut_select$Value) # 2-column matrix, from-to
  r_reclass <- reclassify(sys_raster, lut_rcl)
  #eco_raster <- abs(r_reclass) # Take the absolute value of the ecological distance raster's values
  eraster <- (r_reclass + 1)/2
  eraster[eraster != 1] <- 0
  return(eraster)
  
}

binary_ecoDist <- function(eco_raster, window, window_polygon){
  r_window <- crop(mask(eco_raster, window), window_polygon)
  ematrix <- as.matrix(r_window) # Convert raster to matrix
  #outlist <- list("eco_raster" = r_window, "eco_matrix" = ematrix) # Return raster and matrix
  outraster <- list("eco_matrix" = ematrix) # Return only matrix
  return(outraster)
}


# Combine ecological distance and euclidean distance to make the final weights matrix
# Distance values of 1 = more ecologically similar/closer to pixel
# Distance values of 0 = less ecologically similar/farther from pixel
mapWeights <- function(decay_matrix, eco_matrix){
  weight_matrix <- decay_matrix*eco_matrix
  rw <- raster(weight_matrix)
  w_sum <- sum(rw[], na.rm = TRUE) # Not expected to equal 1
  outlist <- list('wgt_raster' = rw, 'wgt_matrix' = weight_matrix, 'wgt_sum' = w_sum) # Return all things
  return(outlist)
}



# REMOVES non-system pixels from the weight consideration
binary_df_proportions <- function(in_raster, cell_id, system, window_radius, nrv_table){
  # Raster will be defined outside function
  # Make window
  w <- circleWindow(raster = in_raster, cell_id = cell_id, cell_res = 60, pixel_radius = window_radius)
  # Create decaying distance raster in window
  decay <- getDist(pixel_layer = w$Pixel, window = w$Window, window_polygon = w$Window_poly)
  # Create ecological distance raster in window
  eco <- binary_ecoDist(eco_raster = eco_rast, window = w$Window, window_polygon = w$Window_poly)
  #plot(eco$eco_raster)
  # Combine decay distance and ecological distance to make the window weights matrix
  weights <- mapWeights(decay_matrix = decay$decay_matrix, eco_matrix = eco$eco_matrix)
  all_weights <- sum(weights$wgt_matrix, na.rm=TRUE) # Sum of all weights used to proportion each SYSXCLA
  # Adjust proportions so we don't look at non-system pixels in the window
  adj_tbl <- w$tbl[!is.na(w$tbl$SYSXCLA),]
  df <- data.frame()
  for(j in seq(1, nrow(adj_tbl))){
    sysxcla_value <- as.numeric(as.vector(adj_tbl$SYSXCLA[j])) # Get first sysxcla code
    sysxcla_pixels <- which(w$SYSXCLA_matrix == sysxcla_value) # Get IDs of pixels of sysxcla_value
    sysxcla_proportion <- length(sysxcla_pixels)/sum(adj_tbl$Freq)# Get actual proportion of pixels in wth window
    wt <- weights$wgt_matrix[sysxcla_pixels] # Get weights of pixels with sysxcla_Value
    sysxcla_weights <- sum(wt) # Get sum of weights of sysxcla_pixels
    sysxcla_calc <- sysxcla_weights/all_weights # Sum of weights of our SYSXCLA/Sum of all weights in the window = proportion
    sysxcla_out <- c(sysxcla_value, length(sysxcla_pixels), sysxcla_proportion, sysxcla_calc)
    #print(sysxcla_out)
    df <- rbind(df, sysxcla_out)
  }
  colnames(df) <- c("SYSXCLA", "Number of Pixels", "Raw Proportion", "Weighted Proportion")
  # Get weighted counts
  df$`Weighted Count` <- sum(df$`Number of Pixels`)*df$`Weighted Proportion`
  current <- merge(nrv_table, df, by.x = "SYSXCLA", by.y = "SYSXCLA", all=TRUE) # This is kind of moot with the NRV runs we now have...
  current[is.na(current)] <- 0 # Remove any possible NA rows
  new_tbl <- current[!(current$ExpectProp == 0 & current$`Number of Pixels`== 0),] # Remove classes where NRV has 0% and has no pixels
  new_tbl$Pixel_ID <- cell_id
  new_tbl$NRV_ID <- names(in_raster)
  new_tbl <- new_tbl %>% select(NRV_ID, Pixel_ID, SYSXCLA, `Number of Pixels`, `Raw Proportion`,
                                `Weighted Proportion`, `Weighted Count`)
  return(new_tbl)
} 



# Function to make raster of Ecological Departure for any 2 runs
# First need to get df of proprtions in window from above
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


# Create mean ED and variance of ED across all combinations of initial condition X future sims
prep_ED <- function(sys_code, results_file, out_folder){
  fill_raster <- select_sys(sys_code = sys_code, in_raster = tjr) 
  df_weights <- read.csv(results_file)
  df_weights <- df_weights %>% dplyr::select(-c(X))
  colnames(df_weights) <- c("NRV_ID", "Pixel_ID", "SYSXCLA", "Number of Pixels", "Raw Proportion", "Weighted Proportion", "Weighted Count")
  it_names <- unique(df_weights$NRV_ID)[-1] # Remove initial conditions results from NRV list
  nrvi_stack <- stack()
  for(i in it_names){
    print(i)
    nrvi <- run_ED(in_df = df_weights, names = c("TJR_SYSXCLA_60m_50415", i),
                   ref_raster = fill_raster, out_name = i)
    nrvi_stack<- stack(nrvi_stack, nrvi)
  }
  nrvi_mean <- mean(nrvi_stack, na.rm=TRUE)
  names(nrvi_mean) <- paste("spatialed_mean_", sys_code, sep="")
  filename_mean <- paste(paste(out_folder, names(nrvi_mean), sep="\\"), ".tif", sep="")
  writeRaster(nrvi_mean, filename_mean, overwrite=TRUE)
  nrvi_var <- calc(nrvi_stack, fun=var)
  names(nrvi_var) <- paste("spatialed_variance_", sys_code, sep="")
  filename_var <- paste(paste(out_folder, names(nrvi_var), sep="\\"), ".tif", sep="")
  writeRaster(nrvi_var, filename_var, overwrite=TRUE)
  nrvi_out <- stack(nrvi_mean, nrvi_var)
  return(nrvi_out)
}

#-----------------------------------------------------------------
#-----------------------------------------------------------------
# Run on individual systems

# Table with expected proportions for each class
# Either provided via literature/expert review OR create using 700-year NRV runs
snake_nrv <- read.csv('path-to-table')

#-----------------------------------------------------------------
# Ex. PJ woodland

# Initial conditions
syscode <- 10190
snake <- raster("path-to-veg-raster", RAT=TRUE)
snake.select <- select_sys(syscode, snake)


# Expected NRV proportions
nrv.select <- get_nrv(nrv_tbl = snake_nrv, sys_code = syscode)
print(nrv.select)


# Reclassify uncharacteristic class to code 999
snake.reclass <- unchar_reclass(sys_code = syscode, in_raster = snake.select)
plot(snake.reclass, col = colorRamps::blue2red(5))
sum(!is.na(snake.reclass[]))  # Number of pixels to process


# Future simulations 
setwd("path-to-future-simulation-veg-rasters")
rs <- stack(list.files(pattern=".tif"))
nrvx <- select_sys(sys_code = syscode, in_raster = rs)
plot(nrvx[[1:2]], col = colorRamps::blue2red(5)) # Maps should look different
# No need to reclassify - future sims only contain reference classes


# Create ecological distance raster from table of system "likenesses"
# Table provided by ecologists
# Ranges -1 to 1, where -1 is most unlike system, 1 is the same system;
# e.g. PJ woodland similar to curl-leaf mtn mahogany (value = 0.75);
# e.g. PJ woodland most different from mine-active, desert wash, etc. (values = -1)
#################################################
# WILL NEED TO UPDATE FOR SNAKE RANGE NEW SYSTEMS
lut <- read.csv("path\\sys_ecodist_lut.csv")
sys_raster <- raster("TBD_ReferenceOnly_SYS.tif")
eco_rast <- binary_ecoRaster(sys_raster = sys_raster, sys_code = syscode, lut=lut)
plot(eco_rast, col = colorRamps::blue2green(5))


# Stack initial conditions and NRV rasters
x_stack <- stack(snake.reclass, nrvx)
nlayers(x_stack)
names(x_stack)
# Rename layers if helpful


# Create empty data frame to store all pixel data for each NRV
df_weights <- data.frame("NRV_ID" = character(), "Pixel_ID" = character(), "SYSXCLA" = vector(), 
                         "Weighted Count" = numeric(), "Weighted Proportion" = numeric(), 
                         "Raw Proportion" = numeric())

# For each NRV run, calculate ecological departure at each pixel
docells <- which(!is.na(x_stack[[1]][]))
print(length(docells))

# ADJUST WINDOW IF NEEDED
for(s in seq(1, nlayers(x_stack))){
  r <- x_stack[[s]]
  print(names(r))
  for(i in docells){
    pixel <- r[i]
    print(paste(i, pixel, sep=" with code: "))
    df_test <- binary_df_proportions(in_raster = r, cell_id = i, system = syscode, 
                                     window_radius = 15, nrv_table = nrv.select)
    df_weights <- rbind(df_weights, df_test)
  }
}

unique(df_weights$NRV_ID)
write.csv(df_weights, "TBD\\df_weights_nrvs_10190.csv")


#---------------------------------------------------------------------
# apply above to more systems...


#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Use year 0 - mapped TJR raster
snake <- raster("TBD.tif", RAT=TRUE)

# pinon juniper woodland
pjw <- prep_ED(sys_code = 10190, 
               results_file = "TBD\\df_weights_nrvs_10190.csv",
               out_folder = "TBD\\tbd_ed")
plot(pjw)


# END single-system spatial ecological departure metric
# Use outputs in Multi-system ED script if using (spatial_ed_multisystem.R)