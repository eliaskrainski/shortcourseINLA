library(sf)
library(rnaturalearth)
library(fmesher)
library(INLA)
library(ggplot2)
library(fields)

setwd(here::here())

### CRS PROJ4 string to project the BRAZIL map
##  https://epsg.io/29101
crsKM <- st_crs("+proj=poly +lat_0=0 +lon_0=-54 +x_0=5000000 +y_0=10000000 +ellps=aust_SA +towgs84=-67.35,3.88,-38.22,0,0,0,0 +units=km +no_defs +type=crs")

## get the world map
worldMapll <- ne_countries(returnclass = "sf")

plot(worldMapll[,"sov_a3"])

## project the world map
worldMap <- st_transform(
    worldMapll,
    crsKM)

## read the  map
map <- rnaturalearth::ne_countries(country = "Brazil")

boundary_longlat <- st_geometry(map)
st_bbox(boundary_longlat)

boundary <- st_transform(boundary_longlat, crsKM)

bb <- st_bbox(boundary)
bb

rxy <- c(x = diff(bb[c(1,3)]), y = diff(bb[c(2,4)]))
rxy

## project it
map2km <- st_transform(map, crsKM)

tavg <- read.csv("data/tavg_month.csv")

## have a look at the first lines
head(tavg, 2)

## fit a linear regression to check the coeffient for elevation
summary(lm(X202401 ~ elevation, data = tavg))

## again, but with elevation/1000
summary(lm(X202406 ~ I(elevation/1000), data = tavg))

## doing the same with INLA
inla(X202406 ~ I(elevation/1000),
     data = tavg)$summary.fixed

## convert the dataset into spatial dataset and project it
Tavg <- st_transform(
    st_as_sf(
        x = tavg, 
        coords = c(3, 2),
        crs = st_crs(map)),
    crsKM
)

st_bbox(map)
st_bbox(map2km)

## bounding box
bb <- st_bbox(map2km)
bb

ggplot() + theme_minimal() + 
    geom_sf(data = map, fill = 'transparent') +
    geom_sf(aes(fill=X202401), Tavg) 

## visualize the map, and stations with
## size proportional to elevation and color to temperature in week 29
ggplot() + theme_minimal() +
    geom_sf(data = worldMap,
            fill = rgb(165/256,142/256,42/256))+
    geom_sf(data = map, fill = "transparent")+
    geom_sf(aes(color = X202401,
                size = elevation), Tavg) +
    xlim(c(2500, 8000)) + 
    ylim(c(6200, 12000)) + 
    scale_color_distiller(
        palette = "RdBu"
    )

## unite polygons to create a boundary
bnd <- st_union(map2km)

## define the mesh around the boundary
hexpoints <- fm_hexagon_lattice(
    st_buffer(boundary, 100),
    edge_len = 100)
str(hexpoints)

plot(boundary)
points(hexpoints, pch = 4)

mesh <- fm_mesh_2d(
    loc = hexpoints,
    offset = 1000,
    max.edge = 300 
)

mesh$n

plot(mesh)
plot(boundary, add = TRUE, border = "red", lwd = 3)

### visualize 
ggplot() + theme_minimal() +
    geom_sf(data = worldMap,
            fill = rgb(165/256,142/256,42/256))+
    geom_sf(data = map2km) +
    inlabru::gg(mesh) +
    geom_sf(aes(color = X202401,
                size = elevation), Tavg) +
    scale_color_distiller(
        palette = "RdBu"
    ) +
    xlim(c(2500, 8000)) + 
    ylim(c(6200, 12000)) 
    
## define the SPDE model
spdeModel <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(30, 0.05), ## P(range < 30) = 0.05
    prior.sigma = c(1, 0.05)   ## P(sigma > 1) = 0.05
)

## create the projector/mapper from mesh nodes to data locations
Amap <- inla.spde.make.A(
    mesh = mesh,
    loc = st_coordinates(Tavg)
)

dim(Amap)

## eta = \beta_0 + \beta_1 * (Elevation/1000) + efeito.espacial
## y = eta + erro, 

## model formula
smodel <- X202401 ~ I(elevation/1000) +
    f(spatial, model = spdeModel, A.local = Amap)

## prepare the data for INLA: 3 coluns
## 'X202401': temperature at month 1
## 'elevation':
## 'spatial': just to have a column in the dataset
##     and to name it in the outputs,
##     but with NA because the SPDE model will use A.local

dataf <- data.frame(
    tavg[c("X202401", "elevation")],
    spatial = NA ## no index in the data (because it is on the mesh nodes)
)

head(dataf)

## fit the model
sfit <- inla(
    formula = smodel,
    data = dataf
)
sfit$cpu.used

## summary of the fitted values for first data
head(sfit$summary.fitted.values, 2)

## summary of the fixed effects: intercept and elevation/1000
sfit$summary.fixed

## summary of the error noise and the SPDE model parameters
sfit$summary.hyperpar

## plot posterior mean of the fitted values against observed
plot(sfit$summary.fitted.values$mean, dataf$X202401, pch = 8)
cor(sfit$summary.fitted.values$mean, dataf$X202401, use = 'pair')^2

bb

## create a grid/projector/mapper to visualization of the
## spatial effect, and later the fitted values, in a map
bb
bb0 <- list(x = c(2700, 7200), y = c(6200, 10700))
sapply(bb0, diff)/10
projGrid <- fm_evaluator(##)inla.mesh.projector(
    mesh = mesh,
    xlim = bb0$x,
    ylim = bb0$y,
    dims = sapply(bb0, diff)/10
)

str(projGrid)

## conver the grid locations into an sf object
sfGrid <- st_as_sf(
    as.data.frame(projGrid$lattice$loc),
    coords = 1:2,
    crs = crsKM
)

## find which grid points are outside the Ethiopia boundary
igrid.out <- which(sapply(st_within(sfGrid, bnd), length)==0)

## project the posterior mean of the spatial effect into the grid
gridSmean <- fm_evaluate(##)inla.mesh.project(
    projGrid,
    sfit$summary.random$spatial$mean)
gridSmean[igrid.out] <- NA ## make NA grid locations outise the boundary

## visualize the spatial effect posterior mean
image.plot(
    x = projGrid$x,
    y = projGrid$y,
    gridSmean,
    asp = 1
)
plot(st_geometry(map2km), add = TRUE)
plot(st_geometry(worldMap), add = TRUE)\


### ADDD ELEVATION

## download the TOPO dataset: elevation (on land) or depth (on ocean)
if(!file.exists("ETOPO2.RData")) {
    download.file(
        url = paste0(
            "http://leesj.sites.oasis.unc.edu/",
            "FETCH/GRAB/RPACKAGES/",
            "ETOPO2.RData"
        ),
        destfile = "ETOPO2.RData"
    )
}

### load the TOPO dataset
load("ETOPO2.RData")

## extract the longitude (and fixit), and latitude
elon <- attr(ETOPO2, "lon")
elon[elon >= 180] <- 180 - rev(elon[elon >= 180])
elat <- attr(ETOPO2, "lat")

### fix the order of the columns (swap)
ETOPO2 <- ETOPO2[, ncol(ETOPO2):1]

## transform the grid coordinates into latlong
sfGrid.ll <- st_transform(sfGrid, st_crs(map))

ggplot(sfGrid.ll) + geom_sf()

## extract it as a matrix
glocs <- st_coordinates(sfGrid.ll)
summary(glocs)

## find which pixels in the long/lat of TOPO data belong each grid location
ij <- list(
  i = findInterval(glocs[, 1], c(-180, elon + 1 / 60)),
  j = findInterval(glocs[, 2], elat)
)

sapply(ij, summary)

## extract the TOPO data values at the grid locations
etopoll <- sapply(1:nrow(glocs), function(i) ETOPO2[ij$i[i], ij$j[i]])
summary(etopoll)

### make NA the grid locations outside the boundary
etopoll[igrid.out] <- NA

summary(etopoll)

## visualize
ggplot() +
    geom_sf(aes(color = etopoll), data = sfGrid, cex = 2) +
    scale_color_distiller(palette = "RdBu", direction = 1)

## compute the fitted values (mean) for the grid
## by sum of each linear predictor term
gridfitt <- sfit$summary.fixed$mean[1] +
    (etopoll/1000) * sfit$summary.fixed$mean[2] +
    gridSmean

str(gridfitt)

### visualise the computed fitted value
par(mfrow = c(1,1), mar = c(0,0,0,0))
image.plot(
    x = projGrid$x,
    y = projGrid$y,
    gridfitt,
    asp = 1
)
plot(st_geometry(map2km), add = TRUE)
plot(st_geometry(worldMap), add = TRUE)

## in order to compute the uncertainty for the predictions
## make a prediction scenario and re-evaluate the model

## make a prediction dataset
preddf <- data.frame(
  X202401 = NA, ## response set to NA so it will be predicted
  elevation = sapply(
      1:nrow(glocs),
      function(i) ETOPO2[ij$i[i], ij$j[i]]), ## elevation 
  spatial = NA ## the 'spatial' column
)

dim(projGrid$proj$A)

## A.local now need to include the data and prediction projector matrices
smodel_pred <- X202401 ~ I(elevation/1000) +
    f(spatial, model = spdeModel,
      A.local = rbind(Amap, projGrid$proj$A))

## consider the fitted mode
sfit$mode$theta

## re-evaluate the model at the fitted mode
sfit.prd <- inla(
  formula = smodel_pred,
  data = rbind(dataf, preddf), 
  control.mode = list(
      theta = sfit$mode$theta,
      restart = FALSE)
)

## create an index for the prediction
i.pred <- nrow(dataf) + 1:nrow(preddf)

head(sfit.prd$summary.fitted.values[i.pred, ])

length(projGrid$x)

## organize the predicted values (posterior mean) into the matrix
## which represent the grid 
pred.grid <- matrix(
  sfit.prd$summary.fitted.values$mean[i.pred], 
  length(projGrid$x)
)
pred.grid[igrid.out] <- NA ## values outside boundary st to NA

## organize the standard error of the predictions to visualize
sd.pred.grid <- matrix(
  sfit.prd$summary.fitted.values$sd[i.pred], 
  length(projGrid$x)
)
sd.pred.grid[igrid.out] <- NA ## values outside boundary st to NA

### visualize the predicted (posterior mean) temperature
par(mfrow = c(1,2), mar = c(0,0,0,0))
image.plot(
  x = projGrid$x,
  y = projGrid$y,
  pred.grid,
  asp = 1,
  axes = FALSE
)
plot(st_geometry(map2km), add = TRUE)
plot(st_geometry(worldMap), add = TRUE)
## visualize the standard error of the prediction
image.plot(
  x = projGrid$x,
  y = projGrid$y,
  sd.pred.grid,
  asp = 1,
  axes = FALSE
)
plot(st_geometry(map2km), add = TRUE)
plot(st_geometry(worldMap), add = TRUE)
points(st_geometry(Tavg), pch = 8)

## NOTICE that the error (std dev of the prediction)
## near the stations is smaller 
