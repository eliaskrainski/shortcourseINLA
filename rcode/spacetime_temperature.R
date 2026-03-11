library(sf)
library(rnaturalearth)
library(fmesher)
library(INLA)
library(ggplot2)
library(fields)
library(INLAspacetime)

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

## convert the dataset into spatial dataset and project it
Tavg <- st_transform(
    st_as_sf(
        x = tavg, 
        coords = c(3, 2),
        crs = st_crs(map)),
    crsKM
)

## some visualization
jjd <- 5:ncol(tavg)
(nt <- length(jjd))
(ns <- nrow(tavg))

par(mfrow = c(1, 1), mar = c(0,0,0,0))
plot(st_geometry(Tavg))
stlines(t(tavg[,jjd]), as(Tavg, "Spatial"))
plot(worldMap, add = TRUE, col = 'transparent')

## define a temporal mesh
tmesh <- fm_mesh_1d(loc = 1:nt)

## unite polygons to create a boundary
bnd <- st_union(map2km)

## define the mesh around the boundary
hexpoints <- fm_hexagon_lattice(
    st_buffer(boundary, 100),
    edge_len = 120)
str(hexpoints)

plot(boundary)
points(hexpoints, pch = 4)

smesh <- fm_mesh_2d(
    loc = hexpoints,
    offset = 300,
    max.edge = 300 
)
smesh$n

plot(smesh)

## define the spacetime model
## https://github.com/eliaskrainski/INLAspacetime
stmodel <- stModel.define(
    smesh = smesh, ## spatial mesh
    tmesh = tmesh, ## temporal mesh
    model = '220', ## model, see the paper
    ## model = "102" ## separable spacetime of 1st order over time
    control.priors = list(
        prs = c(300, 0.05), ## P(spatial range < 300) = 0.05
        prt = c(5, 0.05), ## P(temporal range < 5) = 0.05
        psigma = c(1, 0.1) ## P(sigma > 1) = 0.1
    )
)

## prepare the data in the long format
colnames(tavg) ## data is in columns 5 to 16
longdf <- data.frame(
    xloc = rep(st_coordinates(Tavg)[, 1], nt),
    yloc = rep(st_coordinates(Tavg)[, 2], nt),
    time = rep(1:nt, each = ns),
    elevation = rep(tavg$elevation, nt),
    temp = unlist(tavg[, jjd])
)

## projection matrix
Aproj <- fm_basis(
  fm_tensor(list(space = smesh, time = tmesh)),
  loc = list(
    space = cbind(longdf$xloc, longdf$yloc),
    time = longdf$time
  )
)

## prior for the precision
pprec <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))

## model formula
fmodel <- temp ~ I(elevation/1000) +
##    f(time, model = 'rw2', scale.model = TRUE, hyper = pprec) +
    f(st, model = stmodel, A.local = Aproj)
longdf$st <- NA

## model fitting
stfit <- inla(
    formula = fmodel,
    data = longdf,
    verbose = TRUE,
    control.mode = list(
        theta = c(0.23, 7.3, 2.16, 0.97),
        restart = TRUE
    )
    control.family = list(hyper=pprec),
    num.threads = "8:1") ### 2 x num. hyperpar

stfit$mode$theta
stfit$cpu.used

## summary of the fixed effects: intercept and elevation/1000
stfit$summary.fixed

## summary of the error noise and the SPDE model parameters
stfit$summary.hyperpar

## plot posterior mean of the fitted values against observed
plot(stfit$summary.fitted.values$mean, longdf$temp, pch = 8)
cor(stfit$summary.fitted.values$mean, longdf$temp, use = 'pair')^2

## spatial range
srange <- inla.tmarginal(exp, stfit$internal.marginals.hyperpar[[2]])
inla.zmarginal(srange)

## temporal range
trange <- inla.tmarginal(exp, stfit$internal.marginals.hyperpar[[3]])
inla.zmarginal(trange)

if(FALSE) {
    ## plot time random effect
    plt <- function(s, ...) {
        n <- nrow(s)
        ry <- range(s$"0.025quant", s$"0.975quant")
        plot(s$ID, s$mean, type = "l", ylim = ry, ...)
        polygon(s$ID[c(1:n, n:1, 1)],
                c(s$"0.025quant", rev(s$"0.975quant"),
                  s$"0.025quant"[1]),
                col = gray(0.5, 0.5), border = gray(0.5, 0.5))
    }
    plt(stfit$summary.random$time)
}

## create a grid/projector/mapper to visualization of the
## spatial effect, and later the fitted values, in a map
bb0 <- list(x = c(2700, 7200), y = c(6200, 10700))
sapply(bb0, diff)/10
projGrid <- fm_evaluator(##)inla.mesh.projector(
    mesh = smesh,
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

## organize the spacetime effect at the mesh nodes in a matrix
st.mean <- matrix(stfit$summary.random$st$mean, ncol = tmesh$n)

## project and organize the spacetime posterior mean into an array
## which represent the grid at each time (an 3D array)
st.array <- array(
    NA,
    c(length(projGrid$x), length(projGrid$y), tmesh$n)
)
for(k in 1:tmesh$n) {
    st.array[, , k] <- inla.mesh.project(projGrid, st.mean[, k])
    st.array[,,k][igrid.out] <- NA ## values outside boundary st to NA
}

## visualize the spacetime effect
par(mfrow = c(3,4), mar = c(0,0,0,0))
for(k in 1:tmesh$n) {
    image.plot(
        x = projGrid$x,
        y = projGrid$y,
        st.array[, , k],
        asp = 1,
        axes = FALSE
    )
    plot(st_geometry(worldMap), add = TRUE)
}


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

## extract the longitude (and fix it), and latitude
elon <- attr(ETOPO2, "lon")
elon[elon >= 180] <- 180 - rev(elon[elon >= 180])
elat <- attr(ETOPO2, "lat")

### fix the order of the columns (swap)
ETOPO2 <- ETOPO2[, ncol(ETOPO2):1]

## transform the grid coordinates into latlong
sfGrid.ll <- st_transform(sfGrid, st_crs(worldMapll))
## extract it as a matrix
glocs <- st_coordinates(sfGrid.ll)

## find which pixels in the long/lat of TOPO data belong each grid location
ij <- list(
  i = findInterval(glocs[, 1], c(-180, elon + 1 / 60)),
  j = findInterval(glocs[, 2], elat)
)

## extract the TOPO data values at the grid locations
etopo.grid <- sapply(1:nrow(glocs), function(i) ETOPO2[ij$i[i], ij$j[i]])
summary(etopo.grid)

### make NA the grid locations outside the boundary
etopo.grid[igrid.out] <- NA
summary(etopo.grid)

## visualize
ggplot() +
    geom_sf(aes(color = etopo.grid), data = sfGrid, cex = 2) +
    scale_color_distiller(palette = "RdBu", direction = 1, na.value = 'transparent')

## in order to compute the uncertainty for the predictions
## make a prediction scenario and re-evaluate the model

## make a prediction dataset
head(longdf)
preddf <- data.frame(
    xloc = rep(projGrid$lattice$loc[, 1], tmesh$n),
    yloc = rep(projGrid$lattice$loc[, 2], tmesh$n),
    time = rep(tmesh$loc, each = nrow(projGrid$lattice$loc)),
    elevation = rep(etopo.grid, tmesh$n),
    temp = NA, ## set temp to NA so to predict it
    st = NA
)

## projector matrix 
Apred <- fm_basis(
   fm_tensor(list(space = smesh, time = tmesh)),
   loc = list(
     space = cbind(preddf$xloc, preddf$yloc),
     time = preddf$time
   )
)

## consider the (previous) fitted mode
stfit$mode$theta

fmodel.prd <- temp ~ I(elevation/1000) +
##    f(time, model = 'rw2', scale.model = TRUE, hyper = pprec) +
    f(st, model = stmodel,
      A.local = rbind(Aproj, Apred))

## re-evaluate the model at the fitted mode
stfit.pred <- inla(
    formula = fmodel.prd,
    data = rbind(longdf, preddf),
    control.family = list(hyper = pprec),
    control.mode = list(
        theta = stfit$mode$theta,
        fixed = TRUE
    )
)

stfit.pred$cpu.used

## create an index for the prediction
i.pred <- nrow(longdf) + 1:nrow(preddf)
length(i.pred)
length(projGrid$x) * length(projGrid$y) * tmesh$n

## organize the predicted values (posterior mean) into an array
## which represent the grid at each time (an 3D array)
pred.array <- array(
  stfit.pred$summary.fitted.values$mean[i.pred], 
  c(length(projGrid$x), length(projGrid$y), tmesh$n)
)
for(k in 1:tmesh$n)
    pred.array[,,k][igrid.out] <- NA ## values outside boundary st to NA

## organize the standard error of the predictions to visualize
sd.pred.array <- array(
  stfit.pred$summary.fitted.values$sd[i.pred], 
  c(length(projGrid$x), length(projGrid$y), tmesh$n)
)
for(k in 1:tmesh$n)
    sd.pred.array[,,k][igrid.out] <- NA ## values outside boundary st to NA

colors1 <- rgb(1:64/64, 0.7, 64:1/64)
bk <- seq(min(pred.array, na.rm = TRUE)-0.001,
          max(pred.array, na.rm = TRUE)+0.001,
          length = 65)

### visualize the predicted (posterior mean) temperature
par(mfrow = c(3,4), mar = c(0,0,0,0))
for(k in 1:tmesh$n) {
    image.plot(
        x = projGrid$x,
        y = projGrid$y,
        pred.array[, , k],
        asp = 1,
        axes = FALSE,
        colors = colors1,
        breaks = bk
    )
    plot(st_geometry(bnd), add = TRUE)
    plot(st_geometry(worldMap), add = TRUE)
}

## visualize the standard error of the prediction (for one time, other times are similar)
x11()
image.plot(
    x = projGrid$x,
    y = projGrid$y,
    sd.pred.array[, , 1],
    asp = 1,
    axes = FALSE
)
plot(st_geometry(bnd), add = TRUE)
plot(st_geometry(worldMap), add = TRUE)
points(st_geometry(Tavg), pch = 8)

## NOTICE that the error (std dev of the prediction)
## near the stations is smaller 
