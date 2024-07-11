#Need to write DPS algorithm for each country
#Write R script with build including packages that need to be installed
#Need to convert to R Script
#Break each cell into functions and then call them in one at a time
#Use relative paths
#add script for loading in libraries and functions
#Start first with running from command line
#comment out parallelism, use in serial for now to make sure code is working
#R script conversion: 
#Terra

# install.packages("raster")
# install.packages("RItools")

#!/usr/bin/env Rscript

# This global processing script PART II is derived from the global processing notebook 
#the input can be the iso3 code (3-character) for one or multiple countries 

#mamba install -c conda-forge r-terra r-optmatch r-sp r-rgdal r-sf r-rgeos r-dplyr r-plyr r-ggplot2 r-raster r-mapview r-stringr r-maptools r-gridExtra r-lattice r-MASS r-foreach r-doParallel r-rlang r-tidyr r-magrittr r-viridis r-ggmap r-Hmisc r-hrbrthemes r-spatialEco r-bit64 r-randomForest r-modelr r-ranger r-caret r-rgeos r-ritools

options(warn=-1)
options(dplyr.summarise.inform = FALSE)

packages <- c("sp","rgdal","sf","rgeos","dplyr","plyr","ggplot2","raster","mapview","stringr",
              "maptools","gridExtra","lattice","MASS","foreach","optmatch","doParallel","RItools","gdalUtils",
              "rlang","tidyr","magrittr","viridis","ggmap","Hmisc","hrbrthemes","spatialEco","bit64","randomForest", "modelr","ranger","caret")
package.check <- lapply(packages, FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = TRUE))
})

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  
  iso3 <- args[1]  #country to process
  gediwk <- args[2]   #the # of weeks GEDI data to use
  mproc <- as.integer(args[3])#the number of cores to use for macthing 
}

#Get Info for GPKG Table Name
# ogrInfo(paste("~/GEDI_PA/Matching_Layers/SEN/SEN_Tiles/N14.30728831742391W-13.52951992770953.geojson_L2A.gpkg",sep=""))

# Subset GEDI data to first 3 features for running troubleshooting
#Run only once!
# sf<- st_read("~/GEDI_PA/Matching_Layers/SEN/WDPA_gedi_l2a+l2b_clean2_SEN/SEN_admin_L2A.gpkg",query = "select * from SEN_admin_L2A limit 3;")
# print(sf)
# st_write(sf,"~/GEDI_PA/Matching_Layers/SEN/WDPA_gedi_l2a+l2b_clean2_SEN/SEN_L2A_subset.gpkg")

# iso3 <- "SEN"
# gediwk<-24
# mproc <-1

cat("Step 0: Loading global variables to process country", iso3,"with GEDI data until week", gediwk, "\n")

f.path <- cat("~/GEDI_PA/Matching_Layers/SEN/")
ecoreg_key <- read.csv(paste("~/GEDI_PA/wwf_ecoregions_key.csv",sep=""))
allPAs <- readRDS(paste("~/GEDI_PA/Matching_Layers/SEN/SEN_PA_poly.rds",sep=""))
MCD12Q1 <- raster(paste("~/GEDI_PA/GEDI_ANCI_PFT_r1000m_EASE2.0_UMD_v1_projection_defined_6933.tif",sep=""))
projection(MCD12Q1) <- sp::CRS(paste("+init=epsg:",6933,sep=""))
world_region <- raster(paste("~/GEDI_PA/GEDI_ANCI_CONTINENT_r1000m_EASE2.0_UMD_v1_revised_projection_defined_6933.tif",sep=""))
projection(world_region) <- sp::CRS(paste("+init=epsg:",6933,sep=""))
adm <- readOGR(paste("~/GEDI_PA/Matching_Layers/SEN/SEN_admin.geojson"),verbose=FALSE)
adm_prj <- spTransform(adm, "+init=epsg:6933") 
load("~/GEDI_PA/rf_noclimate.rdata")
source("~/GEDI_PA/matching_func-Copy1.r")
# flag <- "don't ran extraction"
flag <- "run all"
# flag <- "run remaining"

#---------------STEP5. GEDI PROCESSING - using GEDI shots to extract the treatment/control status, also extract the MODIS PFT for AGB prediction---------------- 
# if (file.exists(paste(f.path,"WDPA_GEDI_extract/",iso3,"_wk",gediwk,"/",iso3,"_gedi_extracted_matching_wk",gediwk,".RDS", sep=""))){
cat(paste("Step 5: Performing WK ",gediwk,"GEDI extraction for", iso3,"\n"))
#matched_all <-read.csv(paste(f.path,"WDPA_extract4_residual_PAs/", iso3, "_wk_", gediwk, "_null_matches_rerun.csv",sep="")) 
matched_all<-list.files(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_points/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = FALSE)
# readr::parse_number() %>% unique()
#list.files(paste(f.path,"WDPA_matching_results/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = FALSE)
registerDoParallel(3)
matched_PAs <- foreach(this_rds=matched_all, .combine = c, .packages=c('sp','magrittr', 'dplyr','tidyr','raster')) %dopar% {   #non-NA matched results
  matched_PAs=c()
  print(this_rds)
  if(nchar(iso3)>3){
    id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[4]  
  } else {
    id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[3]
  }
  matched <- readRDS(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/SEN_wk24/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep=""))
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

# UPDATED FOR ALL OF SEN
cat(paste("Step 5: Performing WK ",gediwk,"GEDI extraction for", iso3,"\n"))
#matched_all <-read.csv(paste(f.path,"WDPA_extract4_residual_PAs/", iso3, "_wk_", gediwk, "_null_matches_rerun.csv",sep="")) 
matched_all<-list.files(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/",iso3,"_wk",gediwk,sep=""),full.names = FALSE)
# readr::parse_number() %>% unique()
#list.files(paste(f.path,"WDPA_matching_results/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = FALSE)
registerDoParallel(3)
matched_PAs <- foreach(this_rds=matched_all, .combine = c, .packages=c('sp','magrittr', 'dplyr','tidyr','raster')) %dopar% {   #non-NA matched results
  matched_PAs=c()
  print(this_rds)
  if(nchar(iso3)>3){
    id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[4]  
  } else {
    id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[3]
  }
  matched <- readRDS(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/SEN_wk24/",iso3,"_pa_", id_pa,"_matching_results_wk",gediwk,".RDS", sep=""))
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

if(flag=="run all"){  #determine how many PAs to run the extraction process
  matched_PAs <- matched_PAs
  cat("Step 5: runing extraction on all", length(matched_PAs),"of non-NA matched results in", iso3,"\n")
} else if (flag=="run remaining"){
  pattern1 = c(paste("wk",gediwk,sep=""),"RDS")
  extracted_PAid <- list.files(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/",iso3,"_wk",gediwk,"/",sep=""), full.names = F, pattern=paste0(pattern1, collapse="|"))%>%
    readr::parse_number() %>% unique()
  matched_PA_id <- matched_PAs %>% readr::parse_number()
  runPA_id <- matched_PA_id[!(matched_PA_id %in% extracted_PAid)]
  if (length(runPA_id)>0){
    Pattern2 <-  paste(runPA_id, collapse="|")
    runPA <-  matched_PAs[grepl(Pattern2,matched_PAs)]
    # runPA_ind <- str_detect(matched_PAs, paste(runPA_id, collapse = "|"))
    matched_PAs <-runPA
  } else {
    matched_PAs <- NULL
    cat("Step 5 already done for", iso3, "\n")
  }
}

source("~/GEDI_PA/matching_func-Copy1.r")

# SUB VERSION
# mproc=2
# registerDoParallel(cores=round(mproc))
# getDoParWorkers()
# startTime <- Sys.time()
# foreach(this_rds=sub, .combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr','raster')) %dopar% {

#     cat("Extracting for no. ", match(this_rds,sub),"pa out of", length(sub),"\n")
#     if(nchar(iso3)>3){
#         id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[4]  
#     } else {
#         id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[3]
#     }
#     print(id_pa)
#     matched <- readRDS(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/SEN_wk24/",iso3,"_pa_",id_pa,"_matching_results_wk24.RDS", sep=""))
#     if (is.null(matched)==TRUE  | nrow(matched)==0) {
#         cat("Matched result is null for PA", id_pa, "quitting...\n")
#     } else if (!is.null(matched)==TRUE){
#     mras  <- tryCatch(matched2ras(matched),
#                       error=function(cond){
#                         message(cond)
#                         cat("Matched result is likely null for country", iso3,"pa", id_pa, "dimension of the match is", dim(matched),"\n")
#                         # writeLines("Matched results is likely null for country", paste(f.path,"WDPA_log/",iso3,"_log_matching.txt", sep=""))
#                         return(NULL)}) #convert the macthed df to a raster stack 
#     print(table(mras$status[]))
#     if(table(mras$status[])[2]==0 | table(mras$status[])[1]==0 | is.null(mras)){
#       cat("Rasterized results unbalanced for PA", id_pa, "quitting...\n")
#     } else {
#       startTime <- Sys.time()
#       iso_matched_gedi<- extract_gedi(matched=matched, mras = mras)#}  #run filtered csvs on mras for extarction

#       tElapsed <- Sys.time()-startTime
#       cat(tElapsed, "for extracting all PAs in", iso3,"\n")
#       cat("Done GEDI for no. ",grep(unique(matched$pa_id), matched_PAs),"pa out of", length(matched_PAs),"\n")

#       iso_matched_gedi <-  iso_matched_gedi %>%
# #         #ASK AMBER ABOUT OG COLUMNS
# #         #THIS NEEDS TO BE ADJUSTED!!!!
#         dplyr::select("pa_id","status",
#                       "wwfbiom","wwfecoreg","shot_number","lon_lowestmode", 
#                       "lat_lowestmode", "lon_lowestmode", 
#                       "lat_lowestmode","rh25", "rh50", "rh75","rh90", "rh98")  #write to individual country folder
#       if (length(unique(iso_matched_gedi$wwfbiom)) >1){
#         pabiome <- iso_matched_gedi$wwfbiom %>% unique() %>% gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1',.,perl = TRUE)%>% str_c( collapse = "+")
#       } else if (length(unique(iso_matched_gedi$wwfbiom))==1){
#         pabiome <- iso_matched_gedi$wwfbiom %>% unique() %>% gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1',.,perl = TRUE)
#       } else {
#         pabiome <- iso_matched_gedi$wwfbiom %>% unique()
#       }
#       # papaddd <- unique(iso_matched_gedi$PADDD) %>% getmode()
#       continent <- unique(iso_matched_gedi$region) %>% getmode()
#       print(paste('output df',dim(iso_matched_gedi)))

#       dir.create(file.path(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/",iso3,"_wk",gediwk,"/",sep="")))
#       saveRDS(iso_matched_gedi, file=paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/SEN_wk24/",iso3,"_pa_", id_pa,"_gedi_wk_",gediwk,"_conti_","biome_",pabiome,".RDS", sep=""))
#       cat(id_pa,"in",iso3,"results is written to dir\n")
#       write.csv(iso_matched_gedi, file=paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/SEN_wk24/","SEN_pa_", id_pa,"_iso_matched_gedi_sub_wk_",gediwk,".csv", sep=""))
#     }
    
#     }
    
#     return(NULL)
# }
# stopImplicitCluster()
# tElapsed <- Sys.time()-startTime
# cat(tElapsed, "for extracting all PAs in", iso3,"\n")
# cat("Done GEDI extraction for pa in ",iso3,"\n")    
    


getwd()

writeLines(c(""), "SENlogJul11.txt")  

# MATCHED_PA VERSION
# NOTE: PA 866 failed since it was over the water so it's been removed from "matching_results"
# mproc=3
# registerDoParallel(cores=round(mproc))
# getDoParWorkers()
# startTime <- Sys.time()

foreach(this_rds=matched_PAs, .combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr','raster')) %dopar% {
    # sink("SENlogJul5_1127.txt", append=TRUE)
    log.text <- paste0(Sys.time(), " processing loop run ", this_rds, length(matched_PAs))
    write.table(log.text, "SENlogJul1.txt", append = TRUE, row.names = FALSE, col.names = FALSE)

    cat("Extracting for no. ", match(this_rds,matched_PAs),"pa out of", length(matched_PAs),"\n")
    # flush.console()
    if(nchar(iso3)>3){
        id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[4]  
    } else {
        id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[3]
    }
    matched <- readRDS(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/SEN_wk24/",iso3,"_pa_",id_pa,"_matching_results_wk24.RDS", sep=""))
    if (is.null(matched)==TRUE  | nrow(matched)==0) {
        cat("Matched result is null for PA", id_pa, "quitting...\n")
    } else if (!is.null(matched)==TRUE){
    mras  <- tryCatch(matched2ras(matched),
                      error=function(cond){
                        message(cond)
                        cat("Matched result is likely null for country", iso3,"pa", id_pa, "dimension of the match is", dim(matched),"\n")
                        # writeLines("Matched results is likely null for country", paste(f.path,"WDPA_log/",iso3,"_log_matching.txt", sep=""))
                        return(NULL)}) #convert the macthed df to a raster stack 
    print(table(mras$status[]))
    if(table(mras$status[])[2]==0 | table(mras$status[])[1]==0 | is.null(mras)){
      cat("Rasterized results unbalanced for PA", id_pa, "quitting...\n")
    } else {
      startTime <- Sys.time()
      iso_matched_gedi<- extract_gedi(matched=matched, mras = mras)#}  #run filtered csvs on mras for extarction
        if (is.null(iso_matched_gedi)) {
        cat("Matched result is null for PA", id_pa, "quitting...\n")}

      tElapsed <- Sys.time()-startTime
      cat(tElapsed, "for extracting all PAs in", iso3,"\n")
      cat("Done GEDI for no. ",grep(unique(matched$pa_id), matched_PAs),"pa out of", length(matched_PAs),"\n")

      iso_matched_gedi <-  iso_matched_gedi %>%
#         #ASK AMBER ABOUT OG COLUMNS
#         #THIS NEEDS TO BE ADJUSTED!!!!
        dplyr::select("pa_id","status",
                      "wwfbiom","wwfecoreg","shot_number","lon_lowestmode", 
                      "lat_lowestmode", "lon_lowestmode", 
                      "lat_lowestmode","rh25", "rh50", "rh75","rh90", "rh98")  #write to individual country folder
      if (length(unique(iso_matched_gedi$wwfbiom)) >1){
        pabiome <- iso_matched_gedi$wwfbiom %>% unique() %>% gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1',.,perl = TRUE)%>% str_c( collapse = "+")
      } else if (length(unique(iso_matched_gedi$wwfbiom))==1){
        pabiome <- iso_matched_gedi$wwfbiom %>% unique() %>% gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1',.,perl = TRUE)
      } else {
        pabiome <- iso_matched_gedi$wwfbiom %>% unique()
      }
      # papaddd <- unique(iso_matched_gedi$PADDD) %>% getmode()
      continent <- unique(iso_matched_gedi$region) %>% getmode()
      print(paste('output df',dim(iso_matched_gedi)))

      dir.create(file.path(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/",iso3,"_wk",gediwk,"/",sep="")))
      saveRDS(iso_matched_gedi, file=paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/SEN_wk24/",iso3,"_pa_", id_pa,"_gedi_wk_",gediwk,"_conti_","biome_",pabiome,".RDS", sep=""))
      write.csv(iso_matched_gedi, file=paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/SEN_wk24/","SEN_pa_", id_pa,"_iso_matched_gedi_sub_wk_",gediwk,".csv", sep=""))
      cat(id_pa,"in",iso3,"results is written to dir\n")
    }
    }
    
    return(NULL)
}


# stopImplicitCluster()
# tElapsed <- Sys.time()-startTime
# cat(tElapsed, "for extracting all PAs in", iso3,"\n")
cat("Done GEDI extraction for pa in ",iso3,"\n")


length(list.files("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/SEN_wk24",pattern = ".csv"))

# fileA<-read.csv("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/SEN_wk24/SEN_pa_352673_iso_matched_gedi_sub_wk_24.csv")
# dim(fileA)

# fileB<-read.csv("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/SEN_wk24/SEN_pa_352642_iso_matched_gedi_sub_wk_24.csv")
# dim(fileB)

########***not used***##########
#SUB VERSION UPDATING FOR ALL RDS
#Pull the CRS from one of the GEDI files so results use 
#Reprojection is clunky and needs to be updated, and overall the projections of above layers do not match the WGS84 projection
# tile<- st_read("~/shared-buckets/abarenblitt/SEN_Tiles/N14.30728831742391W-13.52951992770953.geojson")
#WORKS THROUGH MRAS
# registerDoParallel(cores=round(mproc))
# getDoParWorkers()
# startTime <- Sys.time()

# foreach(this_rds=sub, .combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr','raster')) %dopar% {
#   cat("Extracting for no. ", match(this_rds,sub),"pa out of", length(sub),"\n")
#   if(nchar(iso3)>3){
#     id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[4]  
#   } else {
#     id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[3]
#   }
#   matched <- readRDS(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/SEN_wk24/",iso3,"_pa_",id_pa,"_matching_results_wk24.RDS", sep=""))
#   if (is.null(matched)==TRUE  | nrow(matched)==0) {
#     cat("Matched result is null for PA", id_pa, "quitting...\n")
#   } else if (!is.null(matched)==TRUE){
#     mras  <- tryCatch(matched2ras(matched),
#                       error=function(cond){
#                         message(cond)
#                         cat("Matched result is likely null for country", iso3,"pa", id_pa, "dimension of the match is", dim(matched),"\n")
#                         # writeLines("Matched results is likely null for country", paste(f.path,"WDPA_log/",iso3,"_log_matching.txt", sep=""))
#                         return(NULL)}) #convert the macthed df to a raster stack 
#     if(table(mras$status[])[2]==0 | table(mras$status[])[1]==0 | is.null(mras)){
#       cat("Rasterized results unbalanced for PA", id_pa, "quitting...\n")
#     } else {
#       startTime <- Sys.time()
#       iso_matched_gedi<- extract_gedi(matched=matched, mras = mras)#}  #run filtered csvs on mras for extarction
#       tElapsed <- Sys.time()-startTime

#       cat(tElapsed, "for extracting all PAs in", iso3,"\n")
#       iso_matched_gedi_sub <- iso_matched_gedi %>%
# #         #ASK AMBER ABOUT OG COLUMNS
# #         #THIS NEEDS TO BE ADJUSTED!!!!
#         dplyr::select("pa_id","status",
#                       "wwfbiom","wwfecoreg","shot_number","lon_lowestmode", 
#                       "lat_lowestmode", "lon_lowestmode", 
#                       "lat_lowestmode","rh25", "rh50", "rh75","rh90", "rh98")  #write to individual country folder
#       if (length(unique(iso_matched_gedi$wwfbiom)) >1){
#         pabiome <- iso_matched_gedi$wwfbiom %>% unique() %>% gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1',.,perl = TRUE)%>% str_c( collapse = "+")
#       } else if (length(unique(iso_matched_gedi$wwfbiom))==1){
#         pabiome <- iso_matched_gedi$wwfbiom %>% unique() %>% gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1',.,perl = TRUE)
#       } else {
#         pabiome <- iso_matched_gedi$wwfbiom %>% unique()
#       }
#       # papaddd <- unique(iso_matched_gedi$PADDD) %>% getmode()
#       continent <- unique(iso_matched_gedi$region) %>% getmode()

#       dir.create(file.path(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/",iso3,"_wk",gediwk,"/",sep="")))
#       saveRDS(iso_matched_gedi, file=paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/SEN_wk24/",iso3,"_pa_", id_pa,"_gedi_wk_",gediwk,"_conti_","biome_",pabiome,".RDS", sep=""))
#       cat(id_pa,"in",iso3,"results is written to dir\n")
#       write.csv(iso_matched_gedi, file=paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/SEN_wk24/","SEN_pa_", id_pa,"_iso_matched_gedi_sub_wk_",gediwk,".csv", sep=""))
#     }
#   }
#   return(NULL)
# }

# stopImplicitCluster()
# tElapsed <- Sys.time()-startTime
# cat(tElapsed, "for extracting all PAs in", iso3,"\n")
# cat("Done GEDI extraction for pa in ",iso3,"\n")    

#FULL VERSION UPDATING FOR ALL RDS
#Pull the CRS from one of the GEDI files so results use 
#Reprojection is clunky and needs to be updated, and overall the projections of above layers do not match the WGS84 projection
# tile<- st_read("~/shared-buckets/abarenblitt/SEN_Tiles/N14.30728831742391W-13.52951992770953.geojson")
#WORKS THROUGH MRAS
# registerDoParallel(cores=round(mproc))
# getDoParWorkers()
# startTime <- Sys.time()

# foreach(this_rds=matched_PAs, .combine = foreach_rbind, .packages=c('sp','magrittr', 'dplyr','tidyr','raster')) %dopar% {
#   cat("Extracting for no. ", match(this_rds,matched_PAs),"pa out of", length(matched_PAs),"\n")
#   if(nchar(iso3)>3){
#     id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[4]  
#   } else {
#     id_pa <- this_rds %>% str_split("_") %>% unlist %>% .[3]
#   }
#   matched <- readRDS(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_matching_results/SEN_wk24/",iso3,"_pa_",id_pa,"_matching_results_wk24.RDS", sep=""))
#   if (is.null(matched)==TRUE  | nrow(matched)==0) {
#     cat("Matched result is null for PA", id_pa, "quitting...\n")
#   } else if (!is.null(matched)==TRUE){
#     mras  <- tryCatch(matched2ras(matched),
#                       error=function(cond){
#                         message(cond)
#                         cat("Matched result is likely null for country", iso3,"pa", id_pa, "dimension of the match is", dim(matched),"\n")
#                         # writeLines("Matched results is likely null for country", paste(f.path,"WDPA_log/",iso3,"_log_matching.txt", sep=""))
#                         return(NULL)}) #convert the macthed df to a raster stack 
#     if(table(mras$status[])[2]==0 | table(mras$status[])[1]==0 | is.null(mras)){
#       cat("Rasterized results unbalanced for PA", id_pa, "quitting...\n")
#     } else {
#       startTime <- Sys.time()
#       iso_matched_gedi<- extract_gedi(matched=matched, mras = mras)#}  #run filtered csvs on mras for extarction
#       tElapsed <- Sys.time()-startTime

#       cat(tElapsed, "for extracting all PAs in", iso3,"\n")
#       iso_matched_gedi_sub <- iso_matched_gedi %>%
# #         #ASK AMBER ABOUT OG COLUMNS
# #         #THIS NEEDS TO BE ADJUSTED!!!!
#         dplyr::select("pa_id","status",
#                       "wwfbiom","wwfecoreg","shot_number","lon_lowestmode", 
#                       "lat_lowestmode", "lon_lowestmode", 
#                       "lat_lowestmode","rh25", "rh50", "rh75","rh90", "rh98")  #write to individual country folder
#       if (length(unique(iso_matched_gedi$wwfbiom)) >1){
#         pabiome <- iso_matched_gedi$wwfbiom %>% unique() %>% gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1',.,perl = TRUE)%>% str_c( collapse = "+")
#       } else if (length(unique(iso_matched_gedi$wwfbiom))==1){
#         pabiome <- iso_matched_gedi$wwfbiom %>% unique() %>% gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1',.,perl = TRUE)
#       } else {
#         pabiome <- iso_matched_gedi$wwfbiom %>% unique()
#       }
#       # papaddd <- unique(iso_matched_gedi$PADDD) %>% getmode()
#       continent <- unique(iso_matched_gedi$region) %>% getmode()

#       dir.create(file.path(paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/",iso3,"_wk",gediwk,"/",sep="")))
#       saveRDS(iso_matched_gedi, file=paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/SEN_wk24/",iso3,"_pa_", id_pa,"_gedi_wk_",gediwk,"_conti_","biome_",pabiome,".RDS", sep=""))
#       cat(id_pa,"in",iso3,"results is written to dir\n")
#       write.csv(iso_matched_gedi, file=paste("~/GEDI_PA/Matching_Layers/SEN/WDPA_GEDI_extract4/SEN_wk24/","SEN_pa_", id_pa,"_iso_matched_gedi_sub_wk_",gediwk,".csv", sep=""))
#     }
#   }
#   return(NULL)
# }

# stopImplicitCluster()
# tElapsed <- Sys.time()-startTime
# cat(tElapsed, "for extracting all PAs in", iso3,"\n")
# cat("Done GEDI extraction for pa in ",iso3,"\n")    

# #---------------STEP6: [FIGURE 4B] Calculating per pa summary stats, 1 pa per row, contain shot#/PA---------------------------- 
# f.path<-"~/GEDI_PA/Matching_Layers/SEN/"
# gedi_paf <-list.files(paste(f.path,"WDPA_GEDI_extract4/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = TRUE)
# cat(paste("Step 6: calculating per pa summary stats for", iso3,"\n"))
# # if (file.exists(paste(f.path,"WDPA_GEDI_extract3/pa_stats/",iso3,"_pa_stats_summary_wk",gediwk,".csv", sep=""))) {
# #   #Delete existing files in exists to avoid duplicate appending
# #   cat("old version for", iso3,"exists, removing...\n")
# #   file.remove(paste(f.path,"WDPA_GEDI_extract3/pa_stats/",iso3,"_pa_stats_summary_wk",gediwk,".csv", sep=""))
# # }

# pa_metrics <-readRDS(gedi_paf[3])

# length(table(pa_metrics$status))



# if (flag =="run all"){
#   if(file.exists(paste(f.path,"WDPA_GEDI_extract4/pa_stats/",iso3,"_pa_stats_summary_wk",gediwk,".csv", sep=""))){  #if flag indicates running all steps, we remove the pa_stats for the iso from previous iteration and prompt creation of a new file
#     file.remove(paste(f.path,"WDPA_GEDI_extract4/pa_stats/",iso3,"_pa_stats_summary_wk",gediwk,".csv", sep=""))
#   }
  
#   for (this_paf in gedi_paf){
#     cat(this_paf,"\n")
#     pa_metrics <-readRDS(this_paf)
#     if (length(table(pa_metrics$status))<2) {
#       cat(iso3, this_paf, "has 0 protected or treatment \n")
#     } else if (table(pa_metrics$status)[1]!=0 && table(pa_metrics$status)[2]!=0) {
#       #filter the datafrme by number of shots in each cell, which is euivalent of the number of occurence of each unique UID code
#       tt <- table(pa_metrics$UID)
#       qcellid <- table(pa_metrics$UID)[tt>5] %>% names()
#       pa_metrics_filtered <- pa_metrics %>% dplyr::filter(UID %in% qcellid)
      
#       #calc summary stats for each country 
#       if(nrow(pa_metrics_filtered)>0){
#         pa_stats_summary <- pa_metrics_filtered %>%
#           group_by(status) %>% 
#           dplyr::mutate(pa_id=as.character(pa_id)) %>%
#           dplyr::summarise(pa_id=na.omit(unique(pa_id)),
#                            count=length(rh98),meanrh98=mean(rh98, na.rm=TRUE), sdrh98=sd(rh98, na.rm = TRUE),medrh98=median(rh98, na.rm = TRUE),
#                            meanrh75=mean(rh75,na.rm=TRUE), sdrh75=sd(rh75,na.rm=TRUE), medrh75=median(rh75,na.rm=TRUE),
#                            meanrh50=mean(rh50,na.rm=TRUE), sdrh50=sd(rh50,na.rm=TRUE), medrh50=median(rh50,na.rm=TRUE),
#                            meanrh25=mean(rh25,na.rm=TRUE ), sdrh25=sd(rh25, na.rm=TRUE),medrh25=median(rh25, na.rm=TRUE),
#                            # meanpai=mean(pai, na.rm=TRUE), sdpai=sd(pai, na.rm=TRUE), medpai=median(pai, na.rm=TRUE),
#                            # meancov=mean(cover, na.rm=TRUE), sdcov=sd(cover,na.rm=TRUE),  medcov=median(cover,na.rm=TRUE),
#                            meanagbd=mean(agbd, na.rm=TRUE), sdagbd=sd(agbd, na.rm=TRUE),medagbd=median(agbd, na.rm=TRUE),
#                            wwfecoreg=getmode(wwfecoreg),REGION=getmode(region))%>% 
#           tidyr::pivot_wider(names_from=status, values_from= setdiff(names(.),c("pa_id", "status"))) #writeLine to a large txt file where world pas stats are
#         pa_stats_summary$iso3 <- iso3
#         print(ncol(pa_stats_summary))
#         if(ncol(pa_stats_summary)>30){
#           if(!file.exists(paste(f.path,"WDPA_GEDI_extract4/pa_stats/",iso3,"_pa_stats_summary_wk",gediwk,".csv", sep=""))){
#             print("not exists")
#             write.csv(pa_stats_summary, file=paste(f.path,"WDPA_GEDI_extract4/pa_stats/",iso3,"_pa_stats_summary_wk",gediwk,".csv", sep=""), row.names = FALSE)
#           } else if (file.exists(paste(f.path,"WDPA_GEDI_extract4/pa_stats/",iso3,"_pa_stats_summary_wk",gediwk,".csv", sep=""))){
#             print("exists so appending to existing file")
#             write.table(pa_stats_summary, file=paste(f.path,"WDPA_GEDI_extract4/pa_stats/",iso3,"_pa_stats_summary_wk",gediwk,".csv", sep=""),
#                         sep=",", append=TRUE , row.names=FALSE, col.names=FALSE)   #will not overwrite but append to existing files
#           }
          
#         }
#       }
#     }
#   }
#   cat("Done summarizing pa-level stats for country",iso3,"\n")  
# }

#---------------STEP7: [Figure 4A] Removing dup gedi shots in overlapping region, count shot#/PA w/o dups, results in [extract4/iso_full_nodup]-------------------
# f.path<-"~/GEDI_PA/Matching_Layers/SEN/"
# gedi_paf <-list.files(paste(f.path,"WDPA_GEDI_extract4/",iso3,"_wk",gediwk,sep=""), pattern=".RDS", full.names = TRUE)
# fullds <- data.frame()
# cat ("Step7a: Compiling country wide dataframe for",iso3, "and removing duplicates...\n")

# gedi_paf
# for (c in gedi_paf){  #removing dups based on shot#
#   cat(match(c, gedi_paf),"in total of", length(gedi_paf),"files\n")
#   tmpds <- tryCatch(readRDS(c),error=function(cond){return(NULL)})
#   if(!is.null(tmpds)){
#     tmpds$shot_number=as.character(tmpds$shot_number)
#     tt <- table(tmpds$UID)
#     qcellid <- table(tmpds$UID)[tt>5] %>% names()
#     tmpds_filtered <- tmpds %>% dplyr::filter(UID %in% qcellid)
#     fullds <- rbind(fullds, tmpds_filtered)
#     print(dim(fullds))
#     fullds <- fullds[!duplicated(fullds$shot_number), ]
#     print(dim(fullds))
#   } else {
#     fullds <- fullds
#     print(dim(fullds))
    
#   }
  
# }
# fullds$iso3 <- iso3
# write.csv(fullds, file=paste(f.path,"WDPA_GEDI_extract4/iso_full_nodup/",iso3,"_country_full_nodup_wk",gediwk,".csv", sep=""), row.names = F)
# cat(iso3,"Dup removed df is exported to /iso_full_nodup/ \n")

#---------------STEP7b: Rasterize the non-dup GEDI results to return the shots-per-1km results, results in [extract4/cell_stats/] & [extract4/matched_raster_stack]---------------------------------
# cat("Step 7b: Summarizing #of GEDI shots per 1km pixel\n ")
# # fullds0 <- read.csv(paste(f.path,"WDPA_GEDI_extract3/iso_full_nodup/",iso3,"_country_full_nodup_wk",gediwk,".csv", sep=""))
# iso_gedi_spdf <- SpatialPointsDataFrame(coords=fullds[,c("lon_lowestmode","lat_lowestmode")],
#                                                 proj4string=CRS("+init=epsg:4326"), data=fullds) %>%spTransform(., CRS("+init=epsg:6933"))
# ras <- crop(MCD12Q1, extent(iso_gedi_spdf)) #a little slow with buffer 
# gcount_ras <- rasterize(coordinates(iso_gedi_spdf),ras, fun="count",background=NA)
# names(gcount_ras) <- "gshot_counts"
# gpid_ras <- rasterize(coordinates(iso_gedi_spdf),ras, fun=getmode,field=iso_gedi_spdf$pa_id,background=NA)
# names(gpid_ras) <- "pid"
# gattr_ras <- rasterize(iso_gedi_spdf@coords, ras, fun=getmode, field=iso_gedi_spdf$status, background=NA)
# names(gattr_ras) <- "status"
# gstack <- stack(gcount_ras,gpid_ras,gattr_ras)
# g1km_sp <- as(gstack, 'SpatialPointsDataFrame')
# g1km <- cbind(g1km_sp@coords, g1km_sp@data, country=iso3)
# # dir.create(paste(f.path,"WDPA_GEDI_extract3/cell_stats/",sep=""))
# # writeRaster(gstack,paste(f.path,"WDPA_GEDI_extract4/matched_raster_stack/",iso3,"_cell_shots_pa_wk",gediwk,".tif", sep=""))  #just for later checks as needed
# write.csv(g1km, file=paste(f.path,"WDPA_GEDI_extract4/cell_stats/",iso3,"_cell_shots_wk",gediwk,".csv", sep=""))
# cat(iso3,"1km pixel level shot count df is exported to /cell_stats/ \n")
# rm(ras, gstack)

#---------------STEP8: Calculating per country summary stats, 1 country per row, summarize key stats for the country ---------------------    
# cat("Step 8: Calculating country level summary stats for ", iso3,"\n ")
# # fullds <- read.csv(paste(f.path,"WDPA_GEDI_extract3/iso_full_nodup/",iso3,"_country_full_nodup_wk",gediwk,".csv", sep=""))
# getmode <- function(v) {
#   uniqv <- na.omit(unique(v))
#   uniqv[which.max(tabulate(match(v, uniqv)))]
# }

# # fullds <- read.csv(paste(f.path,"WDPA_GEDI_extract4/iso_full_nodup/",iso3,"_country_full_nodup_wk",gediwk,".csv", sep=""))

# #Currently commented out since L2B needed
# # fullds$pai[!is.finite(fullds$pai)] <- NA

# iso_sum <- fullds %>%
#   group_by(status) %>%  
#   dplyr::summarise(count_ttl=length(rh98),
#                    meanrh98=mean(rh98, na.rm=TRUE), sdrh98=sd(rh98, na.rm = TRUE), medrh98=median(rh98, na.rm = TRUE),msrh98=sum(is.na(rh98)),
#                    # meanpai=mean(pai, na.rm=TRUE), sdpai=sd(pai, na.rm=TRUE),  medpai=median(pai, na.rm=TRUE),mspai=sum(is.na(pai)),
#                    # meancov=mean(cover, na.rm=TRUE), sdcov=sd(cover,na.rm=TRUE),  medcov=median(cover,na.rm=TRUE),mscov=sum(is.na(cov)),
#                    meanagbd=mean(agbd, na.rm=TRUE), sdagbd=sd(agbd, na.rm=TRUE), medagbd=median(agbd, na.rm=TRUE), msagbd=sum(is.na(agbd)))%>% 
#   tidyr::pivot_wider(names_from=status, values_from= setdiff(names(.),c("pa_id", "status")))#writeLine to a large txt file where world pas stats are
# iso_sum$iso3 <- iso3
# continent <- fullds$region %>% unique() %>% getmode()
# iso_sum$continent <- continent

# write.csv(iso_sum, file=paste(f.path,"WDPA_GEDI_extract4/iso_stats/",iso3,"_country_stats_summary_wk",gediwk,"2.csv", sep=""), row.names = F)
# cat(iso3,"country level summary stats is exported to /iso_stats/ \n")

