
#FIRST ACTIVATE r-with-gdal ENVIRONMENT

#For optmatch, paste into terminal:
#First time only, run this first:
#wget -O optmatch_0.10.6.tar.gz https://packagemanager.posit.co/cran/__linux__/bullseye/latest/src/contrib/optmatch_0.10.6.tar.gz?r_version=4.3&arch=x86_64

#R CMD INSTALL optmatch_0.10.6.tar.gz


#IN TERMINAL:
#mamba install -c conda-forge r-rlemon r-svd r-sparsem r-survival r-RItools r-geojsonio
# install.packages('RItools')
# install.packages("spatialEco", version='1.3.2')

# library("devtools")
# library('RItools')
# library("sfheaders")
# library("spatialEco")

#Use mamba install -c conda-forge r-____ to install packages
#mamba install -c conda-forge r-terra r-optmatch r-sp r-rgdal r-sf r-rgeos r-dplyr r-plyr r-ggplot2 r-raster r-mapview r-stringr r-maptools r-gridExtra r-lattice r-MASS r-foreach r-doParallel r-rlang r-tidyr r-magrittr r-viridis r-ggmap r-Hmisc r-hrbrthemes r-spatialEco r-bit64 r-randomForest r-modelr r-geojsonio r-ritools r-spatialeco

#!/usr/bin/env Rscript

# This global processing script is derived from the global processing notebook 
#the input can be the iso3 code (3-character) for one or multiple countries 

options(warn=-1)
options(dplyr.summarise.inform = FALSE)

packages <- c("sp","rgdal","sf","rgeos","dplyr","plyr","ggplot2","raster","mapview","stringr","terra",
              "maptools","gridExtra","lattice","MASS","foreach","doParallel","RItools","optmatch",
              "rlang","tidyr","magrittr","viridis","ggmap","Hmisc","hrbrthemes","spatialEco","bit64","randomForest", "modelr","geojsonio")
#as.factor(paste0('r-', packages," "))



package.check <- lapply(packages, FUN = function(x) {
    suppressPackageStartupMessages(library(x, character.only = TRUE))
})


#To test, we define the variables manually. For final version, run the commented out section below
# gediwk <- 24
# iso3 <-"SEN"
# mproc <- 1
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
#Add more arguments in this list, can set og pathfiles as defaults  
  iso3 <- args[1]  #country to process
  gediwk <- args[2]   #the # of weeks GEDI data to use
  mproc <- as.integer(args[3])#the number of cores to use for macthing 
}

cat("Step 0: Loading global variables for", iso3,"with wk", gediwk, "data \n")

#f.path<-"https://maap-ops-workspace.s3.amazonaws.com/shared/abarenblitt/GEDI_PA/Matching_Layers/SEN/"

f.path <- "~/GEDI_PA/Matching_Layers/SEN/"
matching_tifs <- c("wwf_ecoreg","wwf_biomes","d2roads", "dcities","dem","pop_cnt_2000","pop_den_2000","slope", "tt2cities_2000", "wc_prec_1990-1999",
                   "lc2000","wc_tmax_1990-1999","wc_tavg_1990-1999","wc_tmin_1990-1999" )

list.files(path=f.path, pattern=NULL, all.files=FALSE, 
    full.names=FALSE)

#allPAs will need to come from Celio once we have full IUCN list, follow format of Amber's example rds file
#MOVE ALL THE BELOW TO SHARED BUCKET AND CHANGE TO S3
#*OR* Make these arguments and call these in the terminal and pass in these paths as different arguments *BETTER*
#Need to test that these work in S3 paths
ecoreg_key <- read.csv(paste("~/GEDI_PA/wwf_ecoregions_key.csv",sep=""))
allPAs <- readRDS(paste("~/GEDI_PA/SEN_PA_poly.rds",sep=""))
MCD12Q1 <- raster(paste("~/GEDI_PA/","GEDI_ANCI_PFT_r1000m_EASE2.0_UMD_v1_projection_defined_6933.tif",sep=""))
projection(MCD12Q1) <- sp::CRS(paste("+init=epsg:",6933,sep=""))
world_region <- raster(paste("~/GEDI_PA/","GEDI_ANCI_CONTINENT_r1000m_EASE2.0_UMD_v1_revised_projection_defined_6933.tif",sep=""))
projection(world_region) <- sp::CRS(paste("+init=epsg:",6933,sep=""))
adm <- readOGR(paste("~/GEDI_PA/Matching_Layers/SEN/SEN_admin.geojson"),verbose=F)
adm_prj <- spTransform(adm, "+init=epsg:6933") 
load("~/GEDI_PA/rf_noclimate.rdata")
source("~/GEDI_PA/matching_func-Copy1.r")

gedi_folder <- paste("~/GEDI_PA/Matching_Layers/SEN/SEN_Tiles/")
gedi<-list.files(path=gedi_folder, pattern="L2A", all.files=FALSE, 
    full.names=FALSE)
length(gedi)
# gedi_data <- readOGR(list.files(gedi_folder,full.names=TRUE))
#gedi_data <- readOGR(list.files(gedi_folder,full.names=TRUE)[i], "SEN_admin_L2A")

#Narrow down do L2A 
gedi_folder <- paste("~/GEDI_PA/Matching_Layers/SEN/SEN_Tiles/")
length(dir(gedi_folder, pattern="L2A"))

# STEP1. Create 1km sampling grid with points only where GEDI data is available; first check if grid file exist to avoid reprocessing 
if(!file.exists(paste("~/GEDI_PA/Matching_Layers/SEN/*.rds"))){
  cat("Step 1: Creating 1km sampling grid filter GEDI data for", iso3,"\n")
  GRID.lats <- raster("~/GEDI_PA/EASE2_M01km_lats.tiff")
  GRID.lons <- raster("~/GEDI_PA/EASE2_M01km_lons.tiff")
  GRID.lats.adm   <- crop(GRID.lats, adm_prj)
  GRID.lats.adm.m <- raster::mask(GRID.lats.adm, adm_prj)
  GRID.lons.adm   <- crop(GRID.lons, adm_prj)
  GRID.lons.adm.m <- raster::mask(GRID.lons.adm, adm_prj)
  rm(GRID.lats, GRID.lons, GRID.lats.adm, GRID.lons.adm)
  
  #1.3) extract coordinates of raster cells with valid GEDI data in them
  gedi_folder <- paste("~/GEDI_PA/Matching_Layers/SEN/SEN_Tiles/")

  GRID.coords <- data.frame()
  for(i in 1:length(dir(gedi_folder,pattern="L2A"))){
    # print(list.files(gedi_folder)[i])
    gedi_data <- read_sf(list.files(gedi_folder,pattern="L2A",full.names=TRUE)[i]) %>%
      dplyr::select(lon_lowestmode,lat_lowestmode)
    gedi_data<- gedi_data %>% st_drop_geometry()
    gedi_pts  <- SpatialPoints(coords=gedi_data[,c("lon_lowestmode","lat_lowestmode")],
                               proj4string=CRS("+init=epsg:4326"))
    gedi_pts_prj <- spTransform(gedi_pts, "+init=epsg:6933")
    
    gcount_ras <- rasterize(coordinates(gedi_pts_prj),GRID.lons.adm.m , fun="count",background=NA)
    names(gcount_ras) <- "gshot_counts"
    pxid <- raster::extract(gcount_ras,  gedi_pts_prj)
    gedi_pts_prj %>% 
      SpatialPointsDataFrame(., data=data.frame(pxid)) ->gedi_pts_prj_sp
    gedi_pts_prj_sp$pxid[is.na(gedi_pts_prj_sp$pxid)] <- 0
    gedi_pts_prj_sp[gedi_pts_prj_sp$pxid>5,]->gedi_pts_prj_filtered  #change the numeric threshold to filter with a different min # of GEDI shots in each 1km cell
    
    GRID.lons.overlap <- GRID.lons.adm.m[gedi_pts_prj_filtered]
    GRID.lats.overlap <- GRID.lats.adm.m[gedi_pts_prj_filtered]
    
    x.overlap <- GRID.lons.overlap[!is.na(GRID.lons.overlap)]
    y.overlap <- GRID.lats.overlap[!is.na(GRID.lats.overlap)]
    
    xy.overlap <- cbind(x.overlap,y.overlap)
    xy.overlap.clean <- unique(xy.overlap)
    
    GRID.coords <- rbind(GRID.coords, xy.overlap.clean) 
  }
  GRID.for.matching <- SpatialPoints(coords = GRID.coords, proj4string=CRS("+init=epsg:4326"))
  saveRDS(GRID.for.matching, file = paste("~/GEDI_PA/Matching_Layers/SEN/",gediwk,".RDS", sep=""))
} else if (file.exists(paste("~/GEDI_PA/Matching_Layers/SEN/",gediwk,".RDS", sep=""))) {
  cat(paste("STEP 1: Grid file exists, no need to process grids for ",iso3, "\n"))
}

df <- readRDS("~/GEDI_PA/Matching_Layers/SEN/24.RDS")
length(df)

# STEP2. Clip sampling grid to nonPA areas within country & sample raster layers on nonPA grid
cat("Step 2.0: Reading 1k GRID from RDS for " ,iso3, "\n")
GRID.for.matching <- df

#GRID.for.matching2 <- as(GRID.for.matching,'Spatial')
print(GRID.for.matching)
#GRID.pts.nonPA <- GRID.for.matching2 %>% spTransform(.,"+init=epsg:4326")

#head(GRID.pts.nonPA)


if(file.exists(paste("~/GEDI_PA/Matching_Layers/SEN/24.RDS"))){
  cat("Step 2.1: Preparing control dataset for", iso3, "\n")
  GRID.pts.nonPA <-GRID.for.matching %>% spTransform(., "+init=epsg:4326")
    for(i in 1:length(allPAs)){
  #print(allPAs)
  # print(i)
      PA          <- allPAs[i,]
      PA_prj      <- spTransform(PA, "+init=epsg:6933")
      PA_prj_buff <- gBuffer(PA_prj, width = 10000) #10km buffer
      PA2         <- spTransform(PA_prj_buff, "+init=epsg:4326")
  # GRID.pts.nonPA <- GRID.pts.nonPA %>% SpatialPoints(.)
  # print(class(PA2))
  # print(class(GRID.pts.nonPA))
  # print(length(GRID.pts.nonPA))
  
      overlap     <- GRID.pts.nonPA[PA2]
        if(length(overlap)>0){
            GRID.pts.nonPA0 <- st_difference(sf::st_as_sf(GRID.pts.nonPA), sf::st_as_sf(PA2)) ##remove pts inside poly
            GRID.pts.nonPA <- as(GRID.pts.nonPA0$geometry,'Spatial') %>% spTransform(., "+init=epsg:4326")
      }
    }

  nonPA_xy  <- coordinates(GRID.pts.nonPA)
  colnames(nonPA_xy)  <- c("x","y")
  nonPA_spdf  <- tryCatch(SpatialPointsDataFrame(nonPA_xy, data=data.frame(nonPA_xy),
                                        proj4string=CRS("+init=epsg:4326")),
                          error=function(cond){
                            cat("Country too small, after buffer no grid left, so quit processing country", iso3, dim(nonPA_xy),"\n")
                            writeLines("Country too small, after buffer no grid left", paste("~/GEDI_PA/Matching_Layers/SEN/","WDPA_log/",iso3,"_log_control.txt", sep=""))
                            return(quit(save="no"))})

    
  for (j in 1:length(matching_tifs)){
     # print(matching_tifs[j])
    ras <- raster(paste("~/GEDI_PA/Matching_Layers/SEN/",matching_tifs[j],".tif", sep=""))
    ras_ex <- raster::extract(ras, nonPA_spdf@coords, method="simple", factors=FALSE)
    nm <- names(ras)
    nonPA_spdf <- cbind(nonPA_spdf, ras_ex)
    # print(nonPA_spdf)
    names(nonPA_spdf)[j+2] <- matching_tifs[j]
    
  }
  d_control <- nonPA_spdf
  d_control$status <- as.logical("FALSE")
  names(d_control) <- make.names(names(d_control), allow_ = FALSE)
  d_control <- data.frame(d_control) %>%
    dplyr::rename(
      land_cover = lc2000,
      slope = slope,
      elevation = dem,
      popden = pop.den.2000,
      popcnt=pop.cnt.2000,
      min_temp=wc.tmin.1990.1999,
      max_temp=wc.tmax.1990.1999,
      mean_temp = wc.tavg.1990.1999,
      prec = wc.prec.1990.1999,
      tt2city= tt2cities.2000,
      wwfbiom = wwf.biomes,
      wwfecoreg = wwf.ecoreg,
      d2city = dcities,
      d2road = d2roads,
      lon = x,
      lat = y)
 d_control$land_cover <- factor(d_control$land_cover, levels=sequence(7),
                                 labels = c("l1_forest",
                                            "l2_grassland",
                                            "l3_agriculture",
                                            "l4_wetlands",
                                            "l5_artificial",
                                            "l6_other land/bare",
                                            "l7_water"))
  d_control$wwfbiom <- factor(d_control$wwfbiom,
                           levels = as.vector(unique(ecoreg_key[,"BIOME"])),
                           labels = as.vector(unique(ecoreg_key[,"BIOME_NAME"])))
  d_control$wwfecoreg <- factor(d_control$wwfecoreg,
                             levels = as.vector(ecoreg_key[,"ECO_ID"]),
                             labels = as.vector(ecoreg_key[,"ECO_NAME"]))
  
  
  d_control$UID <-  seq.int(nrow(d_control))
    
    
    
  
  saveRDS(d_control, file=paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/SEN_prepped_control_wk.RDS",sep="")) 
  
} else if (file.exists(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/SEN_prepped_control_wk.RDS",sep=""))){
  cat("Step 2.1: preppred control dataset is already exist for", iso3, "no need for reprocessing\n")
}

cat("Step 3.0: Reading 1k GRID from RDS for " ,iso3)

#STEP3. Loop through all PAs in iso3 country:
# - clip sampling grid to each PA
# - sample raster layers on each PA grid
# - save each PA sample into prepped_pa_##.RDS file

cat("Step 3.0: Reading 1k GRID from RDS for " ,iso3)
GRID.for.matching <- readRDS("~/GEDI_PA/Matching_Layers/SEN/24.RDS")

if(length(dir(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_testPAs",sep=""),pattern = paste(gediwk,".RDS",sep="")))==0){
  cat("Step 3.1: Processing prepped PA treatment dataset for ", iso3, "\n")
  for(i in 1:length(allPAs)){
    cat(iso3, i, "out of ", length(allPAs), "\n")
    testPA <- allPAs[i,]
    testPA <- spTransform(testPA, "+init=epsg:4326")
    GRID.pts.testPA <- GRID.for.matching[testPA]
    
    if(length(GRID.pts.testPA)>0){
      testPA_xy <- coordinates(GRID.pts.testPA)
      colnames(testPA_xy) <- c("x","y")
      testPA_spdf <- SpatialPointsDataFrame(testPA_xy, data=data.frame(testPA_xy),
                                            proj4string=CRS("+init=epsg:4326"))
      for (j in 1:length(matching_tifs)){
        ras <- raster(paste("~/GEDI_PA/Matching_Layers/SEN/",matching_tifs[j],".tif", sep=""))
        ras <- crop(ras, testPA)
        ras_ex <- raster::extract(ras, testPA_spdf@coords, method="simple", factors=F)
        nm <- names(ras)
        testPA_spdf <- cbind(testPA_spdf, ras_ex)
        names(testPA_spdf)[j+2] <- matching_tifs[j]
        
      }
      d_pa <- testPA_spdf
      d_pa$status <- as.logical("TRUE")
      d_pa$DESIG_ENG <- testPA$DESIG_ENG
      d_pa$REP_AREA <- testPA$REP_AREA
      d_pa$PA_STATUS <- testPA$STATUS
      d_pa$PA_STATUSYR <- testPA$STATUS_YR
      d_pa$GOV_TYPE <- testPA$GOV_TYPE
      d_pa$OWN_TYPE <- testPA$OWN_TYPE
      d_pa$MANG_AUTH <- testPA$MANG_AUTH
      names(d_pa) <- make.names(names(d_pa), allow_ = FALSE)
      d_pa <- data.frame(d_pa) %>%
        dplyr::rename(
          land_cover = lc2000,
          slope = slope,
          elevation = dem,
          popden = pop.den.2000,
          popcnt=pop.cnt.2000,
          min_temp=wc.tmin.1990.1999,
          max_temp=wc.tmax.1990.1999,
          mean_temp = wc.tavg.1990.1999,
          prec = wc.prec.1990.1999,
          tt2city= tt2cities.2000,
          wwfbiom = wwf.biomes,
          wwfecoreg = wwf.ecoreg,
          d2city = dcities,
          d2road = d2roads,
          lon = x,
          lat = y)
      d_pa$land_cover <- factor(d_pa$land_cover, levels=sequence(7),
                                labels = c("l1_forest",
                                           "l2_grassland",
                                           "l3_agriculture",
                                           "l4_wetlands",
                                           "l5_artificial",
                                           "l6_other land/bare",
                                           "l7_water"))
      d_pa$wwfbiom <- factor(d_pa$wwfbiom,
                          levels = as.vector(unique(ecoreg_key[,"BIOME"])),
                          labels = as.vector(unique(ecoreg_key[,"BIOME_NAME"])))
      d_pa$wwfecoreg <- factor(d_pa$wwfecoreg,
                            levels = as.vector(ecoreg_key[,"ECO_ID"]),
                            labels = as.vector(ecoreg_key[,"ECO_NAME"]))
      
      d_pa$UID <- seq.int(nrow(d_pa))
      saveRDS(d_pa, file = paste(f.path,"WDPA_matching_results/",iso3,"_testPAs","/","prepped_pa_",
                                 testPA$WDPAID,"_wk",gediwk,".RDS", sep="")) 
    }
  }
} else if (length(dir(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_testPAs","/",sep=""),pattern = paste(gediwk,".RDS",sep="")))>0){
  cat("Step 3.1: prepped PA treatment dataset is already exist for ", iso3, "no need for reprocessing\n")
}


test<-readRDS("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/SEN_testPAs/prepped_pa_352618_wk24.RDS")[1,]
colnames(test)

#STEP4. Set up spatial points data frames (control + each PA) for point matching
# if (file.exists(paste(f.path,"WDPA_matching_results/",iso3,"_wk",gediwk,"/",iso3,"_matching_output_wk",gediwk,".RDS", sep=""))){
cat("Step 4: Performing matching for", iso3,"\n")
d_control_local <- readRDS(file=paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_prepped_control_wk.RDS", sep=""))
d_control_local <-d_control_local[complete.cases(d_control_local), ]  #filter away non-complete cases w/ NA in control set

if(!dir.exists(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_wk",gediwk,"/",sep=""))){
  # cat("Matching result dir does not EXISTS\n")
  dir.create(file.path(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_wk",gediwk,"/",sep="")))
  d_PAs <- list.files(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_testPAs/", sep=""), pattern=paste("wk",gediwk,sep=""), full.names=FALSE)
} else if (dir.exists(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_wk",gediwk,"/",sep=""))){   #if matching result folder exists, check for any PAs w/o matched results
  pattern1 = c(paste("wk",gediwk,sep=""),"RDS")
  matched_PAid <- list.files(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_wk",gediwk,"/",sep=""), full.names = FALSE, pattern=paste0(pattern1, collapse="|"))%>%
    readr::parse_number() %>% unique()
  d_PAs<- list.files(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_testPAs/", sep=""), pattern=paste("wk",gediwk,sep=""), full.names=FALSE)
  d_PA_id <- d_PAs %>% readr::parse_number()
  runPA_id1 <- d_PA_id[!(d_PA_id %in% matched_PAid)]
  
  matched_all <- list.files(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = FALSE)
  registerDoParallel(3)
  matched_PAs <- foreach(this_rds=matched_all, .combine = c, .packages=c('sp','magrittr', 'dplyr','tidyr','raster')) %dopar% {   #non-NA matched results
    matched_PAs=c()
    # print(this_rds)
    if(nchar(iso3)>3){
      id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[4]  
    } else {
      id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[3]
    }
    matched <- readRDS(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep=""))
    if(!is.null(matched)){
      if(nrow(matched)!=0){
        matched_PAs=c(matched_PAs,this_rds) 
      }
    }else {
      # print(this_rds)
      matched_PAs=matched_PAs
    }
    return(matched_PAs)
  }
  stopImplicitCluster()
  
  if(!is.null(matched_PAs)){
    fullmatch_ids <- matched_PAs %>% readr::parse_number()
    runPA_id2 <- d_PA_id[!(d_PA_id %in% fullmatch_ids)]
    runPA_id <- c(runPA_id1,runPA_id2)
    
  } else{
    fullmatch_ids <- d_PAs %>% readr::parse_number()
    runPA_id2 <- fullmatch_ids#d_PA_id[!(d_PA_id %in% fullmatch_ids)]
    runPA_id <- c(runPA_id1,runPA_id2)
    
  }
  
  if (length(runPA_id)>0){
    # Pattern2 <-  paste(runPA_id, collapse="|")
    t <- d_PA_id %in% runPA_id
    runPA <-  d_PAs[t]
    d_PAs <- runPA
  } else {
    d_PAs <- NULL
  }
  write.csv(d_PAs, paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/", iso3, "_wk_", gediwk, "_null_matches_rerun.csv",sep=""))
  cat("Step 4: need to rerun ", length(d_PAs),"PAs\n")
}

registerDoParallel(mproc)
# cat("Parallel processing",getDoParWorkers(),"PAs \n")
startTime <- Sys.time()
foreach(this_pa=d_PAs,.combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr','optmatch','doParallel')) %dopar% {
  pa <- this_pa
  id_pa <-pa %>%str_split("_") %>% unlist %>% .[3]
  # cat(id_pa, "in",iso3,"\n")
  cat("No.", match(pa,d_PAs),"of total",length(d_PAs),"PAs in ", iso3, "\n" )
  d_pa <- readRDS(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_testPAs/",pa, sep=""))
  d_filtered_prop <- tryCatch(propensity_filter(d_pa, d_control_local), error=function(e) return(NA))  #return a df of control and treatment after complete cases and propensity filters are applied
  # cat("Propensity score filtered DF dimension is",dim(d_filtered_prop),"\n")
  d_wocat_all <- tryCatch(filter(d_filtered_prop, status),error=function(e) return(NA))
  d_control_all <- tryCatch(filter(d_filtered_prop, !status),error=function(e) return(NA))
  
  n_control <- dim(d_control_all)[1]
  # ids_all <- d_control_all$UID   #seq(1,n_control)
  ids_all0 <- tryCatch(d_control_all$UID, error=function(e) return(NA))
  ids_all <- d_control_all$UID
  set.seed(125)
  # cat("Using number of cores:",getDoParWorkers(),"\n")
  N <- ceiling(nrow(d_wocat_all)/300)
  l <- tryCatch(split(d_wocat_all, sample(1:N, nrow(d_wocat_all), replace=TRUE)),error=function(e) return(NULL))
  # l <- tryCatch(split(d_wocat_all, (as.numeric(rownames(d_wocat_all))-1) %/% 300),error=function(e) return(0))
  
  if (length(l)<900 && length(l)>0 ){
    pa_match <- data.frame()
    for (pa_c in 1:length(l)){
      ids_all <- d_control_all$UID
      cat("chunk",pa_c,"out of ",length(l), "chunks of PA", id_pa,"\n")
      
      d_wocat_chunk <- l[[pa_c]]
      # #sample the control dataset to the size of the sample dataset, keep unsampled ids to iterate until full number of matches found
      n_treatment <- dim(d_wocat_chunk)[1]
      
      t <- ifelse(floor(n_control/n_treatment)<=7, ifelse(floor(n_control/n_treatment)<1, 1,floor(n_control/n_treatment)),7)   #floor(n_control/n_treatment))
      n_sample <- round(n_treatment*t)    #now the n_control is 1.4 times the number of n_treatment, 7 will set the if ststament below to flase
      m_all2_out <- data.frame()
      Bscore <- data.frame()
      n_matches <- 0
      tryCatch(
        while(n_matches < n_treatment){
          n_ids <- length(ids_all)
          # cat("n ids",n_ids,"\n")
          if(n_ids > n_sample){
            set.seed(125)
            sample_ids_bar <- sample(ids_all, n_sample)
            sample_ids <- sample(ids_all0, n_sample)
            d_control_sample <- d_control_all[d_control_all$UID %in% sample_ids,]
            ids_all <-setdiff(ids_all, sample_ids_bar)    #ids_all[-sample_ids]
            # cat("protected uid", head(d_wocat_chunk$UID),"\n")
            # All approaches
            new_d <- tryCatch(rbind(d_wocat_chunk,d_control_sample),error=function(e) return(NULL))
            # new_d <- tryCatch(rbind(d_wocat_chunk,d_control_all),error=function(e) return(NULL))
            #create a smaller distance matrix
            m_all <- tryCatch(match_wocat(new_d, pid=id_pa),error=function(e) return(NULL))
            # m_all <- match_wocat(new_d)
            m_all2 <- tryCatch(m_all[1,],error=function(e) return(NULL))
            # m_all2 <- m_all[1,]
            n_matches_temp <- tryCatch(nrow(m_all2$df),error=function(e) return(NULL))
            # n_matches_temp <- nrow(m_all2$df)
            if(!is.null(n_matches_temp)){
              # n_matches <- n_matches + nrow(m_all2$df)
              m_all2$df$pa_id <- rep(id_pa,n_matches_temp)
              m_all2_out <- rbind(m_all2_out, m_all2$df)
              matched_protected <- m_all2$df %>% dplyr::filter(status==TRUE)
              matched_control <- m_all2$df %>% dplyr::filter(status==FALSE)
              cat("matched_protected", nrow(matched_protected),"\n")
              n_matches <- n_matches + nrow(matched_protected)
              d_wocat_chunk <- d_wocat_chunk[-(match(matched_protected$UID,d_wocat_chunk$UID)),]
              # d_control_all <- d_control_all[-(match(matched_control$UID,d_control$UID)),]
            } 
            # ids_all <-setdiff(ids_all, sample_ids)
            ids_all0 <-setdiff(ids_all0, matched_control$UID)
            # else {
            #   n_treatment <- 0  #if not macthes are found in this sampling
            # }
          } else {n_treatment <- n_matches}
        }, error=function(e) return(NULL))
      # ids_all0 <-setdiff(ids_all0, matched_control$UID)
      match_score <- m_all2_out
      cat(table(match_score$status),"\n")
      pa_match <- rbind(pa_match,match_score)
    }
  } else if (length(l)>=900){
    registerDoParallel(4)
    pa_match <- foreach(pa_c=1:length(l), .combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr','optmatch','doParallel'))%dopar%{
      # cat("Matching treatment chunk", pa_c, "out of", length(l), "for PA", id_pa,"\n")
      cat("chunk",pa_c,"out of ",length(l), "chunks of PA", id_pa,"\n")
      # cat("head control",head(ids_all0),"\n")
      d_wocat_chunk <- l[[pa_c]]
      # #sample the control dataset to the size of the sample dataset, keep unsampled ids to iterate until full number of matches found
      n_treatment <- dim(d_wocat_chunk)[1]
      # cat( "n control", length(ids_all0),"\n")
      t <- ifelse(floor(n_control/n_treatment)<=7, ifelse(floor(n_control/n_treatment)<1, 1,floor(n_control/n_treatment)),7)   #floor(n_control/n_treatment))
      n_sample <- round(n_treatment*t)    #now the n_control is 1.4 times the number of n_treatment, 7 will set the if ststament below to flase
      m_all2_out <- data.frame()
      Bscore <- data.frame()
      n_matches <- 0
      
      tryCatch(
        while(n_matches < n_treatment){
          
          n_ids <- length(ids_all0)
          # cat("n ids",n_ids,"\n")
          if(n_ids > n_sample){
            set.seed(125)
            sample_ids_bar <- sample(ids_all, n_sample)
            sample_ids <- sample(ids_all0, n_sample)
            d_control_sample <- d_control_all[d_control_all$UID %in% sample_ids,]
            ids_all <-setdiff(ids_all, sample_ids)    #ids_all[-sample_ids]
            # cat("protected uid", head(d_wocat_chunk$UID),"\n")
            # All approaches
            new_d <- tryCatch(rbind(d_wocat_chunk,d_control_sample),error=function(e) return(NULL))
            # new_d <- tryCatch(rbind(d_wocat_chunk,d_control_all),error=function(e) return(NULL))
            
            #create a smaller distance matrix
            m_all <- tryCatch(match_wocat(new_d, pid=id_pa),error=function(e) return(NULL))
            # m_all <- match_wocat(new_d)
            m_all2 <- tryCatch(m_all[1,],error=function(e) return(NULL))
            # m_all2 <- m_all[1,]
            n_matches_temp <- tryCatch(nrow(m_all2$df),error=function(e) return(NULL))
            # n_matches_temp <- nrow(m_all2$df)
            if(!is.null(n_matches_temp)){
              # n_matches <- n_matches + nrow(m_all2$df)
              m_all2$df$pa_id <- rep(id_pa,n_matches_temp)
              m_all2_out <- rbind(m_all2_out, m_all2$df)
              matched_protected <- m_all2$df %>% dplyr::filter(status==TRUE)
              matched_control <- m_all2$df %>% dplyr::filter(status==FALSE)
              cat("matched_protected", nrow(matched_protected),"\n")
              n_matches <- n_matches + nrow(matched_protected)
              d_wocat_chunk <- d_wocat_chunk[-(match(matched_protected$UID,d_wocat_chunk$UID)),]
              # d_control_all <- d_control_all[-(match(matched_control$UID,d_control$UID)),]
              # 
            } 
            ids_all0 <-setdiff(ids_all0, matched_control$UID)
            # cat( "n control", length(ids_all0),"\n")
            
            # else {
            #   n_treatment <- 0  #if not macthes are found in this sampling
            # }
          } else {n_treatment <- n_matches}
        }, error=function(e) return(NULL))
      # ids_all0 <-setdiff(ids_all0, matched_control$UID)
      match_score <- m_all2_out
      # cat(table(match_score$status),"\n")
      return(match_score)
    }
    stopImplicitCluster()
  } else{
    pa_match <- NULL
  }
    saveRDS(pa_match, file=paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_wk",gediwk,"/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep=""))
  # cat("Results exported for PA", id_pa,"\n")
  rm(pa_match)
  return(NULL)
}

tElapsed <- Sys.time()-startTime
# cat(tElapsed, "for matching all PAs in", iso3,"\n")
stopImplicitCluster()
cat("Done matching for",iso3,". Finishing...\n")
 

 # writeLines(paste("Full data balanced and exported GEDI extracts using GEDI data until week", gediwk, sep=""), paste(f.path,"WDPA_log/",iso3,"_log_success_wk",gediwk,".txt", sep=""))

# test2<- read.table("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_log/SEN/SEN_pa_2332_matching_used_covar_log_wk24.txt")
# test2

test<- readRDS("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/SEN_wk24/SEN_pa_2332_matching_results_wk24.RDS")
print(length(test$propensity_scoreA))
#NEED TO PULL PROPENSITY SCORE BEFORE MATCHING!!! 

# test2<- readRDS("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_points/SEN_testPAs/prepped_pa_2332_wk24.RDS")
# print(head(test2))

library(sf)
library(ggplot2)

my_sf <- st_as_sf(test, coords = c('lon', 'lat'))

# ggplot(my_sf) + 
#    geom_sf(aes(color = my_sf$status))

 st_write(my_sf, "2332PATest.shp", crs = 4326)

getwd()
