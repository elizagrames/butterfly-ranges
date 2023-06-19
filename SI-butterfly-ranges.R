# This script is to generate approximate species distributions as textured
# range maps by combining point occurrence data (e.g. museum records or 
# observations from participatory science programs) with expert range polygons.

# In cases where expert range polygons are not available, preliminary polygon
# estimates can be used to create a pseudo-expert range for fitting the final 
# range models. Because the input polygons are essentially priors (though not 
# formally -- see more below and in the associated manuscript), the estimated 
# range polygons can be regenerated with new observations to update our 
# expectations about where species occur. We include a comparison of estimated 
# ranges derived from using an expert range or a pseudo-expert polygon as the 
# basis for distance weighting the models. 

# Please direct any questions about the script to Eliza Grames (egrames@unr.edu)

# Functions --------------------------------------------------------------------

# This script uses the following packages:
  # raster
  # smoothr
  # sf
  # maptools
  # dismo
  # rgdal
  # geodata
  # rJava
  # lwgeom

# rJava is a package that interfaces with Java. Before instalilng it you will need
# Java installed. This software is free, and runs on Windows, Linux, and MacOS:
# https://www.oracle.com/java/technologies/downloads/

## to install, uncomment and run the following lines:
  # install.packages("smoothr")
  # install.packages("sf")
  # install.packages("raster")
  # install.packages("maptools")
  # install.packages("dismo")
  # install.packages("rgdal")
  # install.packages("geodata")
  # install.packages("rJava")
  # install.packages("lwgeom")

# This function takes input polygons (x) and converts that object to a uniform
# range polygon by smoothing edges and dropping 'crumbs' or small, isolated areas
# Crumb size (cs) is given in units of square kilometers
smooth_raster <- function(x, #input polygon
                          bw=5, # bandwidth parameter
                          cs=10000*10000){ # crumb size
  x <- smoothr::smooth(x, method="ksmooth", bandwidth=bw)
  x <- sf::st_as_sf(x)
  x <- sf::st_cast(x, "POLYGON")
  sf::sf_use_s2(FALSE)
  x <- smoothr::fill_holes(x, 10000000000)
  x <- smoothr::drop_crumbs(x, cs*100) # changed size of crumbs
  x <- sf::as_Spatial(x)
  x <- raster::aggregate(x, dissolve=T)
  return(x)
}

# This function generates random points over a probability surface
generate_points <- function(x, # a raster of a probability surface
                            n=1000, # number of points to generate
                            crs = NULL){
  if(is.null(crs)){crs = raster::crs(x)}            
  pt.samples <- as.data.frame(dismo::randomPoints(mask=x, n, prob=T))
  sp::coordinates(pt.samples) <- ~x + y
  sp::proj4string(pt.samples) <- crs
  return(pt.samples)
}

# This function streamlines converting a suitability map into a polygon based
# on a threshold of minimum suitability
generate_threshold_poly <- function(predicted.suitability, # raster of suitability predictions from model
                                    threshold, # cutoff for threshold
                                    cs, # crumb size
                                    bw){ # bandwidth
  predicted.suitability[predicted.suitability < threshold] <- NA
  predicted.suitability <- raster::aggregate(raster::rasterToPolygons(predicted.suitability),
                                             dissolve=T)
  predicted.suitability <- smooth_raster(predicted.suitability, bw = bw, cs=cs)
  return(predicted.suitability)
}

# Set up spatial placeholders --------------------------------------------------
data("wrld_simpl", package="maptools")

# Subset to only North America and remove internal country borders
northamerica <- wrld_simpl[wrld_simpl$NAME%in%c("Mexico", "Canada", "United States"),]
northamerica <- raster::aggregate(northamerica, dissolve=T) 

# Set the extent of the bounding box. This is mostly important because the tail
# of the Aleutians is in the eastern hemisphere, but we do not want to model 
# that large of an extent
bb <- raster::extent(northamerica)
bb[2] <- c(-48) # Ignore the tail of the Aleutians

# Covariate layers -------------------------------------------------------------

# International Geosphere Biosphere Program vegetation layer
landcover <- raster::raster("./covariates/igbp")
landcover <- raster::crop(landcover, bb)

# Bioclim variables
clim <- geodata::worldclim_global(var = "bio", res = 10, path = "./covariates/worldclim")
clim <- raster::stack(clim)
# Elevation
dem <- raster::raster("./covariates/dem")
dem[is.na(dem)] <- 0 # Tends to cause issues with coastal species in Florida

# Make our predictor layers have the same extent so we can stack them
dem <- raster::crop(dem, raster::extent(landcover))
clim <- raster::crop(clim, raster::extent(landcover))

# Define a base map of gridcells covering the area of interest
gridcells <- raster::raster(array(dim=dim(dem)[1:2]))
gridcells[is.na(gridcells)] <- 1
raster::values(gridcells) <- 1:raster::ncell(gridcells)
raster::extent(gridcells) <- raster::extent(landcover)
raster::crs(gridcells) <- raster::crs(northamerica)

# Species point occurrence data and expert polygons ----------------------------

# Read in species point occurrence data. For this example species, we have 
# already downloaded and cleaned records, however, users could also call the
# dismo::gbif() function to pull in species records and then clean coordinates.
# Alternatively, users could read in a .csv of coordinates or data stored in any
# other format and convert it to a SpatialPointsDataFrame.
sp.obs <- rgdal::readOGR("./hesperia-juba-gbif.kml") 

# Read in the expert range polygon
expert <- rgdal::readOGR("./hesperia-juba-expert-range.kml")

# Optional: generate a pseudo-expert polygon for the 'two step' approach -------

# If you are working on a taxa or region that lacks expert ranges, you can still
# fit a model that incorporates 'distance from expert range' by creating a 
# pseudo-expert range to serve as the underlying polygon. The logic is rather
# circular and we don't recommend this unless there is no expert opinion to rely
# on, but, this process does produce qualitatively similar results

# First, set up the raster stack of predictor layers
cov.stage1 <- raster::stack(dem, clim, landcover)

# Fit a simple MaxEnt model to the data, specifying that landcover is a factor
mod1 <- dismo::maxent(cov.stage1, p=sp.obs, factors=21) 

# Users will typically get a warning message at this stage, usually due to 
# observations near the coast that have NA values in cells that are mostly ocean

# Make predictions from the model
p1 <- dismo::predict(cov.stage1, mod1)

# Next, we want to create a 'pseudo-expert' polygon from this model by using
# an arbitrary cutoff to subset the predictions just to the 'best' areas and
# converting the resulting raster to a polygon. This is by no means a true 
# substitute for the expertise associated with expert ranges, and is just a 
# quick shortcut to still be able to generate a final estimated range when 
# expert ranges are not available to use as a prior.

# We use 0.5 as our arbitrary cutoff to generate fairly restricted polygons
pseudo <- p1
pseudo[pseudo<0.5] <- NA

# Convert the resulting raster to a polygon
pseudo.outline <- raster::aggregate(raster::rasterToPolygons(pseudo), dissolve=T)
pseudo.smooth <- smooth_raster(pseudo.outline, bw=1)

# Fit species distribution model: with expert range ----------------------------

# First, we calculate distance from the expert range for all cells in the extent
exp.dist <- gridcells
expert.boundaries <- raster::aggregate(expert, dissolve=T)
dists <- rgeos::gDistance(expert.boundaries, as(exp.dist, "SpatialPoints"), 
                       byid=TRUE, hausdorff = F)

# Input those distance values to the expert distance raster
exp.dist[!is.na(exp.dist)] <- dists/max(dists) # scale so max value is 1

# Next, we have to define the distance decay in a way that makes sense for 
# the geographic scope of the species being modeled and their biology (i.e.
# accounting for things like dispersal capability). In the case of butterflies
# of North America, we opted for log(distance+1)^2 because it provided a shallow enough
# rate of decay near the expert range to allow for fuzzy boundaries (the expert
# ranges we rely on are for breeding range only) while forcing steep rates of
# decay at greater distances to minimize the effects of misidentifications and
# to limit predictions within feasible ranges (e.g. to prevent predictions of
# suitable habitat in Florida for species whose ranges are primarily in Mexico
# and the southwestern US because of similar climatic conditions in the Yucatan
# and southern Florida). The choice of decay depends on the context, and may be
# more specific for single species models as opposed to whole taxonomic groups.
d.decay <- max(raster::values(log(exp.dist+1)), na.rm=T)-(log(exp.dist+1))

# Rather than weighting our observations by distance from the expert range, we
# modify the covariate layers such that far away places are dissimilar to
# conditions near the expert range for continuous predictors.
cov.stage2 <- raster::stack(dem*d.decay*d.decay, clim*d.decay*d.decay, landcover)

# Now we fit a MaxEnt model with the modified covariates and predict from that
# model across the extent of our covariates
mod2 <- dismo::maxent(cov.stage2, p=sp.obs, factors=21) 
p2 <- dismo::predict(cov.stage2, mod2)

# Integrated range polygons ----------------------------------------------------

# First we set up the bandwidth and crumb size values depending on the area
# of the expert range. We use larger bandwidths for larger species to smooth 
# out wrinkles, whereas we want to preserve those wobbles for small ranges. We
# also want to use different crumb sizes because species with small ranges
# (e.g. restricted to the Florida Keys) may have such small ranges they would
# be below the threshold to drop from a larger range.
if(raster::area(expert)/10000000000 < 1){
  cs <- 100000000
}else{cs <- 10000000000}

bw <- raster::area(expert)/100000000000
if(bw<1){bw <- 1} # Lower bound so edges aren't too wobbly
if(bw>5){bw <- 5} # Upper bound so ranges still preserve their shape

# Set up three different thresholds for clipping the suitability model
# We are naming these 'ir' for integrated range, to distinguish them from
# range polygons based on no boundaries (nb) and the two-step (ts) pseudo-expert
# approach


## Threshold of 0.05
ir.05 = generate_threshold_poly(predicted.suitability =p2, 
                                  threshold = 0.05, 
                                  cs=cs, bw=bw)
# Next threshold is 0.10 suitability minimum
ir.1 = generate_threshold_poly(predicted.suitability =p2, 
                                  threshold = 0.10, 
                                  cs=cs, bw=bw)
# And finally 0.20 suitability minimum
ir.2 = generate_threshold_poly(predicted.suitability =p2,
                                  threshold = 0.20, 
                                  cs=cs, bw=bw)

# Plot output ------------------------------------------------------------------
colpal <- colorRampPalette(rev(c("#f94144ff", "#f3722cff", "#f8961eff", "#f9c74fff", "#90be6dff", 
                             "#43aa8bff", "#577590ff", "grey95")))(256)
raster::plot(p2, col=colpal, box=F, axes=F)
raster::plot(ir.05, add=T, border="grey80")
raster::plot(ir.1, add=T, border="grey10")
raster::plot(ir.2, add=T, border="black")

# Two-step process -------------------------------------------------------------

# Earlier, we created a pseudo-expert range that could be substituted instead
# of an actual expert range. Let's repeat the modeling process with that to
# see how different the final estimates are. Code is the same as above but with
# comments removed for space. We overwrite intermediary objects (e.g. exp.dist).

pseudo.dist <- gridcells
pseudo.boundaries <- raster::aggregate(pseudo.smooth, dissolve=T)
p.dists <- rgeos::gDistance(pseudo.boundaries, as(exp.dist, "SpatialPoints"), 
                          byid=TRUE, hausdorff = F)

pseudo.dist[!is.na(pseudo.dist)] <- p.dists/max(p.dists) # scale so max value is 1
d.decay <- max(raster::values(log(pseudo.dist+1)), na.rm=T)-(log(pseudo.dist+1))

cov.stage3 <- raster::stack(dem*d.decay*d.decay, clim*d.decay*d.decay, landcover)

mod3 <- dismo::maxent(cov.stage3, p=sp.obs, factors=21) 
p3 <- dismo::predict(cov.stage3, mod3)

if(raster::area(pseudo.smooth)/10000000000 < 1){
  cs <- 100000000
}else{cs <- 10000000000}

bw <- raster::area(pseudo.smooth)/100000000000
if(bw<1){bw <- 1} 
if(bw>5){bw <- 5} 


## Threshold of 0.05
ts.05 = generate_threshold_poly(predicted.suitability = p3, 
                               threshold = 0.05, 
                               cs=cs, bw=bw)
# Next threshold is 0.10 suitability minimum
ts.1 = generate_threshold_poly(predicted.suitability = p3, 
                               threshold = 0.1, 
                               cs=cs, bw=bw)
# And finally 0.20 suitability minimum
ts.2 = generate_threshold_poly(predicted.suitability = p3, 
                               threshold = 0.2, 
                               cs=cs, bw=bw)

# No boundaries or distance decay polygons -------------------------------------

# Earlier, we created a model with no boundaries or distance decay whatsoever
# as the first stage of the two-stage process, even though we know MaxEnt models
# with no expert range polygons will overpredict suitability in areas where
# species do not occur. Let's see how different the estimates are for these
# models from our two-stage and expert-based approaches.

if(raster::area(pseudo.smooth)/10000000000 < 1){
  cs <- 100000000
}else{cs <- 10000000000}

bw <- raster::area(pseudo.smooth)/100000000000
if(bw<1){bw <- 1} 
if(bw>5){bw <- 5} 

## Threshold of 0.05
nb.05 = generate_threshold_poly(predicted.suitability = p1, 
                                threshold = 0.05, 
                                cs=cs, bw=bw)
# Next threshold is 0.10 suitability minimum
nb.1 = generate_threshold_poly(predicted.suitability = p1, 
                               threshold = 0.1, 
                               cs=cs, bw=bw)
# And finally 0.20 suitability minimum
nb.2 = generate_threshold_poly(predicted.suitability = p1, 
                               threshold = 0.2, 
                               cs=cs, bw=bw)

# Comparison of extent of occurrence across approaches -------------------------

par(mfrow=c(3,1))

# First, the expert-based approach
raster::plot(p2, col=colpal, box=F, axes=F, main = "Integrated map",
             legend.args = list(text = "Predicted suitability"))
raster::plot(ir.05, add=T, border="grey80")
raster::plot(ir.1, add=T, border="grey20")
raster::plot(ir.2, add=T, border="black")

# Then the two-stage approach
raster::plot(p3, col=colpal, box=F, axes=F, main ="Two-stage approach",
             legend.args = list(text = "Predicted suitability"))
raster::plot(ts.05, add=T, border="grey80")
raster::plot(ts.1, add=T, border="grey20")
raster::plot(ts.2, add=T, border="black")

# And finally, no polygon / distance decay
raster::plot(p1, col=colpal, box=F, axes=F, main = "No decay / range limitaitons",
             legend.args = list(text = "Predicted suitability"))
raster::plot(nb.05, add=T, border="grey80")
raster::plot(nb.1, add=T, border="grey20")
raster::plot(nb.2, add=T, border="black")

# We get similar, but not identical, range polygons
# across the three approaches for the three different
# cutoff stages, but notice that there are suitable areas
# selected in far-flung places in the no-boundaries version

# Comparison of predicted occurrence points across approaches ------------------

# Another way we can compare performance is the proportion
# of occurrences that would be predicted outside the expert
# range boundaries using the suitability models generated
# from the three approaches as probability surfaces.

# To get distance estimates, we want to use a projected expert range
albers <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83"
expert <- sp::spTransform(expert, CRSobj=albers)

ir.pts <- sp::spTransform(generate_points(p2, crs = raster::crs(p2)), albers)
ts.pts <- sp::spTransform(generate_points(p3, crs = raster::crs(p3)), albers)
nb.pts <- sp::spTransform(generate_points(p1, crs = raster::crs(p1)), albers)

expert$name <- "in.range"

par(mfrow=c(3,1), mar=c(5,4,4,2)+0.1)

# What proportion of 1000 points from the integrated suitability model
# fall inside the original expert range?
table(sp::over(ir.pts, expert)$name)/1000
ir.dists <- rgeos::gDistance(expert, ir.pts, byid = T)
hist(ir.dists[ir.dists>0]/1000, 50, border=F, col=colpal[1], 
     main="", xlim=c(0,4000),
     xlab="Distance from range boundary for points outside expert range (km)")
(med.cur = median(ir.dists[ir.dists>0]/1000))
abline(v = med.cur, col = "blue", lty =2 )
legend(x = "topright", lty = 2, col = 'blue', legend = "median", cex = 2)

# For the two-stage approach?
table(sp::over(ts.pts, expert)$name)/1000
ts.dists <- rgeos::gDistance(expert, ts.pts, byid = T)
hist(ts.dists[ts.dists>0]/1000, 50, border=F, col=colpal[1], 
     main="", xlim=c(0,4000),
     xlab="Distance from range boundary for points outside expert range (km)")
(med.cur = median(ts.dists[ts.dists>0]/1000))
abline(v = med.cur, col = "blue", lty =2 )

# For no restrictions? 
table(sp::over(nb.pts, expert)$name)/1000
nb.dists <- rgeos::gDistance(expert, nb.pts, byid = T)
hist(nb.dists[nb.dists>0]/1000, 100, border=F, col=colpal[1], 
     main="", xlim=c(0,4000),
     xlab="Distance from range boundary for points outside expert range (km)")
mean(nb.dists[nb.dists>0]/1000)
(med.cur = median(nb.dists[nb.dists>0]/1000))
abline(v = med.cur, col = "blue", lty =2 )


# There are many other ways that the different approaches could be compared, 
# these are just a few examples. For applications to single species, 
# researchers should carefully consider thresholds used, the distance decay
# function, and what makes sense for that species. The methods demonstrated
# in this script were used for the integrated ranges across 750 species because
# they, in general, worked well for that set of species. 

# Write outputs ----------------------------------------------------------------

if(!dir.exists("./output")){dir.create("./output")}

raster::writeRaster(p1, file="./output/hesperia-juba-unrestricted-range.tiff")
raster::writeRaster(p2, file="./output/hesperia-juba-integrated-range.tiff")
raster::writeRaster(p3, file="./output/hesperia-juba-twostage-range.tiff")

integrated.polygon <- sf::st_as_sf(ir.2)
sf::write_sf(integrated.polygon, 
             "./output/hesperia-juba-integrated-polygon.kml",
             driver="KML")


