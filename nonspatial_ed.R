# Non-spatial ecological departure of systems

# How out-of-whack is each system across the landscape?
# Compare OBSERVED distribution of classes to EXPECTED districution of classes

library(rgdal)
library(sp)
library(raster)
library(dplyr)
library(rgeos)
library(spdep)
library(colorRamps)
library(readxl)
library(dplyr)
library(data.table)

#-----------------------------------------------------------------
# Functions

select_sys <- function(sys_code, in_raster){
  sys_code_char <- as.character(sys_code)
  sys_code_low <- as.numeric(paste(sys_code_char, "000", sep=""))
  sys_code_high <- as.numeric(paste(sys_code_char, "999", sep=""))
  out_raster <- in_raster
  out_raster[out_raster > sys_code_high] <- NA
  out_raster[out_raster < sys_code_low] <- NA
  return(out_raster)
}

# Get NRV - read in table of NRV proportions and select the system you want
# outputs a table in correct format for calculating ecological departure
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

#-----------------------------------------------------------------
# Get non-spatial ecological departure of user-defined map area
# Developed for Snake Range area, hence "snake"-named variables.

# Raster of SYS and CLA codes
snake <- raster('path-to-raster')
plot(snake)

# Table with expected proportions for each class
# Either provided via literature/expert review OR create using 700-year NRV runs
# Would need to do latter with a different script
snake_nrv <- read.csv('path-to-table')

# May not be interested in every system
# Loop thru systems of interest by 5-digit SYS code
sys_list <- c(10110, 10190, 10620, 10790, 10791, 10801, 
              10802, 10803, 10804, 10810, 10812, 11060, 
              11230, 11260, 11261, 11350, 11450, 11451, 
              11530, 11540, 11543, 11544)

# Empty dataframes to populate
prop_all <- data.frame(SYSXCLA = vector(), Count = vector(),
                       Prop = vector(), ExpectProp = vector())
ns_df_all <- data.frame(SYS_CODE = vector(), NS_ED = vector())

# For each system, get proportions of each class - including uncharacteristic classes
# U-classes take away from reference classes
for(i in sys_list){
  print(i)
  # Get current proportions
  my_raster <- select_sys(sys_code = i, in_raster = snake)
  sys_lo <- as.numeric(paste(i, "100", sep = ""))
  sys_hi <- as.numeric(paste(i, "999", sep = ""))
  my_raster <- reclassify(my_raster, c(sys_lo, sys_hi, sys_hi), include.lowest=TRUE)
  code_df <- data.frame(table(my_raster[]))
  colnames(code_df) <- c("SYSXCLA", "Count")
  code_df$Prop <- code_df$Count/(sum(code_df$Count))

  # Get expected proportions
  nrv <- get_nrv(nrv_tbl = snake_nrv, sys_code = i)
  #nrv <- nrv[-c(which(nrv$ExpectProp==0)),]
  #print(nrv)
  # Compare proportions - get non-spatial ED
  mdf <- full_join(code_df, nrv, by="SYSXCLA")
  rownames(mdf) <- mdf$SYSXCLA
  prop_all <- rbind(prop_all, mdf)
  mdf <- mdf[,-c(1, 2)]
  if(length(which(is.na(mdf$Prop)) > 0)){
    print("replace NA in current proportion...")
    mdf[is.na(mdf$Prop),]$Prop <- 0
  } else{
    print("no changes to current proportion")
  }
  mdf$ExpectProp <- mdf$ExpectProp*0.01
  print(mdf)
  dat <- transform(mdf, min = pmin(Prop, ExpectProp))
  print(dat)
  ns_df <- data.frame(SYS_CODE = i, NS_ED = 1 - sum(dat$min))
  print(ns_df)
  ns_df_all <- rbind(ns_df_all, ns_df)
  # Reclassify raster
  sys_raster <- reclassify(sys_raster, c(i, i, ns_df$NS_ED), include.lowest=TRUE) 
}
out_raster <- reclassify(sys_raster, c(1, 99999999, NA))
plot(out_raster)
writeRaster(out_raster, "TBD.tif",
            overwrite=TRUE)


print(prop_all)
colnames(prop_all)[3] <- "CurrentProp"
prop_all[is.na(prop_all$Count),]$Count <- 0
prop_all[is.na(prop_all$CurrentProp),]$CurrentProp <- 0
prop_all$ExpectProp <- prop_all$ExpectProp*0.01
write.csv(prop_all, "TBD.csv")

write.csv(ns_df_all, "TBD.csv")
# Save to XLSX and color-code to get traditional high/med/low departure table


# END