#### Use this code for PCS by Cell ####
# devtools::install_github("dinilu/paleoCLMs-package")
library(paleoCLMs)
library(reshape2)
library(raster)
library(rgdal)
library(picante)
library(phytools)

# Define parameters
INDIR <- "../Data/"
PERIODS <- seq(0, 21000, by=1000) # That generate a numerical vector with time slides. I use this afterwards to name files and objects
CLIMVARS <- c("tmin_low_month", "tmax_high_month", "prcp_low_month", "prcp_high_month","aet_year_ave", "etr_year_ave", "wdei_year_ave")
CLIMVARS_NEW_NAMES <- c("Tmin", "Tmax", "Pmin", "Pmax", "AET", "ETR", "WDI")
VARS <- c(CLIMVARS, "Deglac.")
VARS_NEW_NAMES <- c(CLIMVARS_NEW_NAMES, "Deglac.")
DISTANCES <- c(120, 360, 480)
PCS_METRICS <- c("nri", "nti")

# Load pollen data and remove low quality data
pollen.original <- lapply(PERIODS, loadPollen, INDIR)
pollen.hq <- lapply(pollen.original, removeLowQualitySamples, 0.75)

# Average values for those pollen sites in the same grid cell
raster.tmp <- raster("../Data/Climate/CCSM/1000BP/tmin_low_month.tif")
pollen.abun <- lapply(pollen.hq, reduceDuplicated, raster.tmp, weighted=F)
rm(raster.tmp)

# Load climate data
clim <- mapply(loadClim, period=PERIODS, pollen=pollen.abun, MoreArgs=list(clim_model="CCSM", indir=INDIR, vars=CLIMVARS))
clim <- lapply(clim, as.data.frame)
names(clim) <- paste(PERIODS, "BP", sep="")
  
# Remove incomplete cases from the pollen and the climate datasets
complete.index <- lapply(clim, FUN=function(x){which(complete.cases(x))})
pollen.abun <- mapply(FUN=function(x,i){x[i,]}, pollen.abun, complete.index, SIMPLIFY=F)
clim <- lapply(clim, FUN=function(x){x[which(complete.cases(x)),]})

# Extract coordinates to a new object and remove them from the pollen data.frame
coord <- lapply(pollen.abun, FUN=function(x){x[,which(colnames(x) %in% c("x","y"))]})
names(coord) <- paste0(PERIODS, "BP")
pollen.abun <- lapply(pollen.abun, FUN=function(x){x[,-which(colnames(x) %in% c("x","y"))]})

# Convert pollen concentrations to presences/absences matrix given a specific pollen threshold(Community matrix)
# var05 <- lapply(pollen.abun, pollenThreshold, 0.05)
# pollen <- mapply(applyThreshold, pollen.abun, threshold=var05, SIMPLIFY=F)
# rm(var05)

pollen <- pollen.abun

# Load ice-sheet data and extract ice-sheet info for pollen-sites
directory <- paste0(INDIR, "Dyke2004")
files <- list.files(directory, full.names=T, pattern=".tif")

ice <- stack(files)
ice <- list(ice[["ice1k"]], ice[["ice1k"]], 
                ice[["ice2k"]], ice[["ice3k"]], 
                ice[["ice4k"]], ice[["ice4k"]], 
                ice[["ice55k"]],  ice[["ice6k"]],
                ice[["ice72k"]], ice[["ice8k"]],
                ice[["ice9k"]], ice[["ice96k"]], 
                ice[["ice1025k"]], ice[["ice11k"]],
                ice[["ice12k"]], ice[["ice125k"]], 
                ice[["ice135k"]], ice[["ice14k"]], 
                ice[["ice15k"]], ice[["ice155k"]],
                ice[["ice165k"]], ice[["ice175k"]])
names(ice) <- paste(PERIODS, "BP", sep="")
ice <- stack(ice)
rm(directory, files)

sites.ice <- lapply(coord, FUN=raster::extract, x=ice)


# Calculate age since deglaciation
get_deglaciation_time <- function(period, list) {
  period_bp <- paste0(period, "BP")
  i <- which(colnames(list[[period_bp]]) == paste0("X", period_bp))
  if(i == 1){
    mat <- list[[period_bp]]
  }else{
    mat <- subset(list[[period_bp]], select=-c(1:i-1))
  }
  mat <- cbind(mat, X22000=1)
  v <- apply(mat, MARGIN=1, FUN=function(x){
    x[which(is.na(x))] <- 0
    x <- max(which(x != 1))
    return(x)})
  v <- v + 1
  txt <- colnames(mat)[v]
  txt <- gsub("X", "", txt) 
  txt <- gsub("BP", "", txt) 
  num <- as.numeric(txt)
  deglac <- num - period
  deglac <- data.frame(Deglac.=deglac)
  return(deglac)
}

deglac <- lapply(PERIODS, get_deglaciation_time, sites.ice)
names(deglac) <- paste0(PERIODS, "BP")

  
# Load phylogenetic tree 
genus.tree <- read.newick("Results/Data/angiosperm_pollen_tree.nwk")


# NRI & NTI
# Net Relatedness Index (Webb) (nri = -1 * mpd.obs.z)
# Nearest Taxon Index (nti = -1 * mntd.obs.z)
# nri < 0 and mpd.obs.p >= 0.95 indicates phylogenetic evenness
# nri > 0 and mpd.obs.p <= 0.05 indicates phylogeneic clustering
# nti < 0 and mntd.obs.p >= 0.95 indicates phylogenetic evenness
# nti > 0 and mntd.obs.p <= 0.05 indicates phylogeneic clustering

# Create a function to do the analysis
calculate_pcs <- function(comm, phy, null.model="taxa.labels", abundance.weighted=FALSE, runs = 999, iterations = 1000) {
  # Prepare data
  commData <- as.data.frame(comm)
  
  # pruning phylogenetic tree based on species present in the community dataset 
  tempData <- match.phylo.comm(phy, commData)
  prunedTree <- tempData$phy
  prunedComm <- tempData$comm 
  
  # Calculate cophenetic distance matrix
  phydist <- cophenetic(phy)
  
  # Calculating NRI and NTI and storing results
  output.nri <- ses.mpd(prunedComm, phydist, null.model = null.model, abundance.weighted = abundance.weighted, runs = runs, iterations = iterations)
  output.nti <- ses.mntd(prunedComm, phydist, null.model = null.model, abundance.weighted = abundance.weighted, runs = runs, iterations = iterations)
  output <- cbind(output.nri, output.nti)
  return(output)
}

# Calculate PCS metrics for all pollen sites and all ages.
pcs.metrics <- lapply(pollen, FUN=calculate_pcs, genus.tree, null.model = "independentswap", abundance.weighted = TRUE)
names(pcs.metrics) <- paste0(PERIODS, "BP")

# save(pcs.metrics, file="Results/PCS/pcs_metrics.RData")
# load("Results/PCS/pcs_metrics.RData")


# Arrange data in a single dataframe for subsequent analysis
PCS_COLNAMES <- colnames(pcs.metrics[[1]])
COORD_COLNAMES <- colnames(coord[[1]])

# Combine dataframes with pcs, variables (climate and deglaciation), and coordinates
pcs.data <- mapply(cbind, pcs.metrics, clim, deglac, coord, SIMPLIFY = FALSE)

# Melt all dataframes in a single dataframe for further analysis and plots 
pcs.data <- melt(pcs.data, id.vars=c(PCS_COLNAMES, VARS, COORD_COLNAMES))
colnames(pcs.data) <- c(PCS_COLNAMES, VARS_NEW_NAMES, COORD_COLNAMES, "year")
pcs.data$year <- factor(pcs.data$year, levels=paste(PERIODS, "BP", sep=""))

# Calculate PCS metrics in the "standard" metric 
pcs.data$nri <- -1 * pcs.data$mpd.obs.z
pcs.data$nti <- -1 * pcs.data$mntd.obs.z

# Calculate PCS p-values 
pcs.data$nri_p <- pcs.data$mpd.obs.p
pcs.data$nti_p <- pcs.data$mntd.obs.p

# Calculate significant PCS metrics 
pcs.data$nri_significant <- pcs.data$nri_p < 0.05 | pcs.data$nri_p > 0.95 
pcs.data$nti_significant <- pcs.data$nti_p < 0.05 | pcs.data$nti_p > 0.95 

# Remove NAs 
pcs.data <- na.omit(pcs.data)

# Save final dataframe 
write.csv(pcs.data, file="Results/PCS/PCS_Clim_Coords.csv")
