
library(sf)
library(ggplot2)
library(INLAspacetime)

setwd(here::here(""))

downloadUtilFiles(data.dir = ".", year = 2024)

### read the informations on stations
all.stations.longlat <- read.fwf(
    file = "ghcnd-stations.txt", 
    widths = diff(c(0, 11, 20, 30, 37, 40, 71, 75, 79, 85)),
    comment.char = ""
)

### set colnames
colnames(all.stations.longlat) <- c(
  "station", "latitude", "longitude",
  "elevation", "state", "name", "gsn", "hcn/crn", "wmo"
)

allst <- st_as_sf(
    all.stations.longlat,
    coords = 3:2,
    crs = st_crs("+proj=longlat +datum=WGS84")
)

wd0 <- ghcndSelect("2024.csv.gz")/10

str(wd0)
dim(wd0)

map <- rnaturalearth::ne_countries(country = "Brazil")

plot(map)

## select stations with dist=5e5 around the map
ii <- which(sapply(st_is_within_distance(allst, map, dist = 1e5), length)>0)
str(ii)

ii.d <- allst$station[ii] %in% dimnames(wd0)[[1]]
table(ii.d)

ggplot() + theme_minimal() +
    geom_sf(data = map) +
    geom_sf(aes(size = elevation,
                color = ii.d),
            data = allst[ii, ]) + labs(color = "2024")


id0 <- pmatch(allst$station[ii], dimnames(wd0)[[1]])
id <- id0[complete.cases(id0)]

dates <- as.Date(dimnames(wd0)[[2]], "%Y%m%d")
range(dates)

plot(dates, apply(is.na(wd0[id,,]), 2, sum))

tt <- 1:dim(wd0)[2]##182:236
range(tt)
head(dates[tt], 3)
tail(dates[tt], 3)

plot(dates[tt], apply(is.na(wd0[id,tt,]), 2, sum))

apply(is.na(wd0[id, , ]), 3, table)
apply(is.na(wd0[id, tt, ]), 3, table)

summary(apply(is.na(wd0[id, tt, ]), c(1,3), sum))

j <- 6
par(mfrow = c(1, 1), mar = c(3, 3, 0.1, 0.1), mgp = c(2, 0.5, 0), bty = "n")
plot(dates, wd0[id[j], , 1], type = "l", col = "blue",
     xlab = '', ylab = "Temperature",
     ylim = range(wd0[id[j], , ], na.rm = TRUE))
points(dates, wd0[id[j], , 1], pch = 19, col = "blue")
points(dates, wd0[id[j], , 2], pch = 19, col = "green")
lines(dates, wd0[id[j], , 2], col = "green")
points(dates, wd0[id[j], , 3], pch = 19, col = "red")
lines(dates, wd0[id[j], , 3], col = "red")

##apply(is.na(wd0[id, , ]), c(1,3), sum)

par(mfrow = c(2, 2), mar = c(0.1,0.1,0.1,0.1), mgp = c(2, 0.5, 0), bty = "n")
for(v in 1:3) {
    plot(st_coordinates(allst[ii[ii.d], ]), pch = 19, asp = 1,
         axes = FALSE, xlab = '', ylab = '')
    box()
    plot(st_geometry(map), add = TRUE)
    stlines(t(wd0[id, tt, v]), as(allst[ii[ii.d],], "Spatial"), yscale = 0.5, cex = 0.1)
}

par(mfrow = c(1,1), mar = c(0.1,0.1,0.1,0.1), mgp = c(2, 0.5, 0), bty = "n")
plot(st_coordinates(allst[ii[ii.d], ]), pch = 19, asp = 1,
     axes = FALSE, xlab = '', ylab = '')
plot(st_geometry(map), add = TRUE)
stlines(t(wd0[id, tt, 2]), as(allst[ii[ii.d],], "Spatial"), yscale = 0.5, cex = 0.5)

tavg.day <- data.frame(
    all.stations.longlat[ii[ii.d], 1:4],
    round(wd0[id,tt,2], 2)
)

tavg.day[1:3, 1:9]
tavg.day[1:3, -3:0+ncol(tavg.day)]

write.csv(tavg.day, "data/tavg_day.csv", row.names = FALSE)

### average over months
tm <- apply(wd0[id,,2], 1, tapply, substr(dimnames(wd0)[[2]], 1, 6), mean, na.rm = TRUE)
dim(tm)

## organize and save
tavg.m <- data.frame(
    all.stations.longlat[ii[ii.d], 1:4],
    round(t(tm), 2)
)

dim(tavg.m)

write.csv(tavg.m, "data/tavg_month.csv", row.names = FALSE)

## average over weeks
format(dates, "%V")

tw <- apply(wd0[id,,2], 1, tapply, format(dates, "%V"), mean, na.rm = TRUE)

par(mfrow = c(10, 9), mar = c(1.5, 1.5, .2, 0.2))
for(i in 1:90) {
    plot.ts(tw[, i])
    points(tw[, i], pch = 19, cex = 0.3)
}

par(mfrow = c(1,1), mar = c(0.1,0.1,0.1,0.1), mgp = c(2, 0.5, 0), bty = "n")
plot(st_coordinates(allst[ii[ii.d], ]), pch = 19, asp = 1,
     axes = FALSE, xlab = '', ylab = '')
plot(st_geometry(map), add = TRUE)
stlines(tw[18:40, ], as(allst[ii[ii.d],], "Spatial"), yscale = 0.5, cex = 0.1)

format(as.Date("2025-05-01"), "%V")
format(as.Date("2025-10-01"), "%V")

dim(tw)

tavg.df <- data.frame(
    all.stations.longlat[ii[ii.d], 1:4],
    round(t(tw[21:34, ]), 2)
)

tavg.df

write.csv(tavg.df, "tavg.csv", row.names = FALSE)
