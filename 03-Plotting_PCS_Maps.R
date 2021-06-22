# load all packages
library(ggplot2)
library(gridExtra)
library(reshape2)
library(raster)
library(ggpubr)
library(dismo)

# Functions
load_all_clim <- function (period, indir, vars){
  require(raster)
  directory <- paste(indir, "Climate/CCSM/", period, sep = "") 
  file.list <- paste(directory, "/", vars, ".tif", sep = "")
  clim.stack <- stack(file.list)
  return(clim.stack)
}

shp2df <- function(shp){
  require(plyr)
  require(rgeos)
  require(maptools)
  require(broom)
  shp@data$id <- rownames(shp@data)
  shp.points <- tidy(shp, region="id")
  shp.df <- join(shp.points, shp@data, by="id")
  colnames(shp.df) <- c("x", "y", "order", "hole", "piece", "group", "id", "layer")
  return(shp.df)
}

# Subset dataframe for plotting
data.cels <- pcs.data
data.cels <- data.cels[,c("nri", "nti", "nri_p", "nti_p", "year", "x", "y", "nti_significant", "nri_significant")]

# Load climate rasters 
clim.raster <- lapply(paste(PERIODS, "BP", sep=""), FUN=load_all_clim, indir=INDIR, vars=CLIMVARS)
names(clim.raster) <- paste(PERIODS, "BP", sep="")
rm(load_all_clim)


# LOAD EXTRA DATA FOR PLOTTING 

# Create shapefiles maps with land boundaries
map.raster <- lapply(clim.raster, FUN=function(x){x[[2]]/x[[2]]})
map.raster[[1]] <- map.raster[[2]]
map.shapefile <- lapply(map.raster, FUN=function(x, y){rasterToPolygons(x[[1]], dissolve=y)}, y=TRUE)

# Prepare the data for ggplot
map.df <- lapply(map.shapefile, shp2df)
map.df.melt <- melt(map.df, id.vars=c("x", "y", "order", "hole", "piece", "group", "id", "layer"))
colnames(map.df.melt) <- c("x", "y", "order", "hole", "piece", "group", "id", "layer", "year")
map.df.melt$year <- factor(map.df.melt$year, levels=paste(PERIODS, "BP", sep=""))

# Load shapefiles for ice sheet - Try to replace from the previous in 02-script.R
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
rm(directory, files)

# Prepare the ice data for ggplot
ice.shapefile <- lapply(ice, FUN=function(x, y){rasterToPolygons(x[[1]], dissolve=y)}, y=T)
ice.df <- lapply(ice.shapefile, shp2df)
ice.df.melt <- melt(ice.df, id.vars=c("x", "y", "order", "hole", "piece", "group", "id", "layer"))
colnames(ice.df.melt) <- c("x", "y", "order", "hole", "piece", "group", "id", "layer", "year")
ice.df.melt$year <- factor(ice.df.melt$year, levels=paste(PERIODS, "BP", sep=""))

# Create subset for publication
selected.periods <- c(seq(1000, 10000, by=3000), seq(12000, 21000, by=2000))
selected.periods <- paste0(selected.periods, "BP")


# Create convex hull polygons for ploting results

library(sp)
library(gstat)
spatial.data.cels <- data.cels
coordinates(spatial.data.cels) <- ~ x + y
proj4string(spatial.data.cels) <- "+proj=longlat +ellps=WGS84"

get_grids <- function(raster){
  coords <- raster %>% coordinates %>% as.data.frame
  coordinates(coords) <- ~ x + y 
  values <- extract(raster, coords)
  grids <- cbind(as.data.frame(coords), values)
  grids <- subset(grids, !is.na(values))
  grids <- subset(grids, values == 1)
  coordinates(grids) <- ~ x + y  
  proj4string(grids) <- "+proj=longlat +ellps=WGS84"
  return(grids)
} 

# ch <- convHull(data.cels.subset[chull(data.cels.subset[, c("x", "y")] ), c("x", "y")] )
ch <- convHull(data.cels[chull(data.cels[, c("x", "y")] ), c("x", "y")] )
maps <- lapply(map.raster, function(x, y){predict(y, x)}, ch)
maps <- mapply(overlay, maps, map.raster, MoreArgs=list(fun=function(x,y){(x*y)}), SIMPLIFY=F)
interp.grids <- lapply(maps, get_grids)

# Calculate inverse distance weighted interpolation within convex hulls
run_idw <- function(period, grid, data, var){
  idw <- idw(formula(paste0(var, "~ 1")), subset(data, year=paste0(period, "BP")), grid)
  idw
} 

nri.idws <- mapply(run_idw, PERIODS, interp.grids, MoreArgs=list(data=spatial.data.cels, var="nri"), SIMPLIFY=FALSE)
nri.idws <- lapply(nri.idws, as.data.frame)

nti.idws <- mapply(run_idw, PERIODS, interp.grids, MoreArgs=list(data=spatial.data.cels, var="nti"), SIMPLIFY=FALSE)
nti.idws <- lapply(nti.idws, as.data.frame)

names(nri.idws) <- names(nti.idws) <- paste0(PERIODS, "BP")

nri.idws.melt <- melt(nri.idws, id.vars=c("x", "y", "var1.pred"))
nri.idws.melt <- nri.idws.melt[,c("x", "y", "var1.pred", "L1")] 
colnames(nri.idws.melt) <- c("x", "y", "nri", "year")
nri.idws.melt$year <- factor(nri.idws.melt$year, levels=paste(PERIODS, "BP", sep=""))

nti.idws.melt <- melt(nti.idws, id.vars=c("x", "y", "var1.pred"))
nti.idws.melt <- nti.idws.melt[,c("x", "y", "var1.pred", "L1")] 
colnames(nti.idws.melt) <- c("x", "y", "nti", "year")
nti.idws.melt$year <- factor(nti.idws.melt$year, levels=paste(PERIODS, "BP", sep=""))


# Create NRI maps
nri.map <- ggplot(nri.idws.melt, aes(x=x, y=y)) +
  geom_polygon(data=subset(map.df.melt, hole == FALSE), aes(group=group), fill="grey", color="black") +
  geom_raster(data=nri.idws.melt, aes(fill=nri)) +
  geom_polygon(data=subset(map.df.melt, hole == FALSE), aes(group=group), fill="000000AA", color="black") +
  geom_polygon(data=subset(map.df.melt, hole == TRUE), aes(group=group), fill="steelblue1", color="black") +
  geom_polygon(data=subset(ice.df.melt, layer == 1), aes(group=group), fill="white", color="black") +
  geom_polygon(data=subset(ice.df.melt, layer == 2), aes(group=group), fill="steelblue1", color="black") +
  scale_fill_gradient2(low="red", mid="#ffffbf", high="blue", name="", limits=c(-2.15, 2.15), breaks=c(-2, 0, 2)) +
  facet_wrap(~year, ncol=5) +
  coord_equal(ylim=c(25,61), xlim=c(-113, -53)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_rect(fill='steelblue1'))  

# Create NTI maps
nti.map <- ggplot(nti.idws.melt, aes(x=x, y=y)) +
  geom_polygon(data=subset(map.df.melt, hole == FALSE), aes(group=group), fill="grey", color="black") +
  geom_raster(data=nti.idws.melt, aes(fill=nti)) +
  geom_polygon(data=subset(map.df.melt, hole == FALSE), aes(group=group), fill="000000AA", color="black") +
  geom_polygon(data=subset(map.df.melt, hole == TRUE), aes(group=group), fill="steelblue1", color="black") +
  geom_polygon(data=subset(ice.df.melt, layer == 1), aes(group=group), fill="white", color="black") +
  geom_polygon(data=subset(ice.df.melt, layer == 2), aes(group=group), fill="steelblue1", color="black") +
  scale_fill_gradient2(low="red", mid="#ffffbf", high="blue", name="", limits=c(-2.15, 2.15), breaks=c(-2, 0, 2)) +
  facet_wrap(~year, ncol=5) +
  coord_equal(ylim=c(25,61), xlim=c(-113, -53)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_rect(fill='steelblue1'))  

# Create NRI subset maps 
nri.map.ss <- ggplot(subset(nri.idws.melt, year %in% selected.periods), aes(x=x, y=y)) +
  geom_polygon(data=subset(map.df.melt, hole == FALSE & year %in% selected.periods), aes(group=group), fill="grey", color="black") +
  geom_raster(data=subset(nri.idws.melt, year %in% selected.periods), aes(fill=nri)) +
  geom_polygon(data=subset(map.df.melt, hole == FALSE & year %in% selected.periods), aes(group=group), fill="000000AA", color="black") +
  geom_polygon(data=subset(map.df.melt, hole == TRUE & year %in% selected.periods), aes(group=group), fill="steelblue1", color="black") +
  geom_polygon(data=subset(ice.df.melt, layer == 1 & year %in% selected.periods), aes(group=group), fill="white", color="black") +
  geom_polygon(data=subset(ice.df.melt, layer == 2 & year %in% selected.periods), aes(group=group), fill="steelblue1", color="black") +
  scale_fill_gradient2(low="red", mid="#ffffbf", high="blue", name="", limits=c(-2.15, 2.15), breaks=c(-2, 0, 2)) +
  facet_wrap(~year, ncol=3) +
  coord_equal(ylim=c(25,61), xlim=c(-113, -53)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_rect(fill='steelblue1'))  

# Create NTI subset maps 
nti.map.ss <- ggplot(subset(nti.idws.melt, year %in% selected.periods), aes(x=x, y=y)) +
  geom_polygon(data=subset(map.df.melt, hole == FALSE & year %in% selected.periods), aes(group=group), fill="grey", color="black") +
  geom_raster(data=subset(nti.idws.melt, year %in% selected.periods), aes(fill=nti)) +
  geom_polygon(data=subset(map.df.melt, hole == FALSE & year %in% selected.periods), aes(group=group), fill="000000AA", color="black") +
  geom_polygon(data=subset(map.df.melt, hole == TRUE & year %in% selected.periods), aes(group=group), fill="steelblue1", color="black") +
  geom_polygon(data=subset(ice.df.melt, layer == 1 & year %in% selected.periods), aes(group=group), fill="white", color="black") +
  geom_polygon(data=subset(ice.df.melt, layer == 2 & year %in% selected.periods), aes(group=group), fill="steelblue1", color="black") +
  scale_fill_gradient2(low="red", mid="#ffffbf", high="blue", name="", limits=c(-2.15, 2.15), breaks=c(-2, 0, 2)) +
  facet_wrap(~year, ncol=3) +
  coord_equal(ylim=c(25,61), xlim=c(-113, -53)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_rect(fill='steelblue1'))  

# Arrange maps in composites for publication
pcs.maps <- ggarrange(nri.map, nti.map, ncol=1, labels=c("a)", "b)"), common.legend = TRUE, legend="right")
pcs.maps.ss <- ggarrange(nri.map.ss, nti.map.ss, ncol=1, labels=c("a)", "b)"), common.legend = TRUE, legend="right")

# Save figures as .pdf files 
ggexport(pcs.maps, filename="Results/Manuscript_figures_and_tables/Figure_A1.pdf", width=9, height=14)
ggexport(pcs.maps.ss, filename="Results/Manuscript_figures_and_tables/Figure_1.pdf", width=7, height=10)
