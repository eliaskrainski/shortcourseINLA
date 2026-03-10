
setwd(here::here(""))
getwd()

nasc <- read.csv2(
    file = "data/nasc.csv", skip = 3,
    nrows = 855, encoding = 'latin1',
    na.strings = '-')

head(nasc,2)
tail(nasc,2)

iobt <- read.csv2(
    file = "data/iobt.csv", skip = 3,
    nrows = 855, encoding = 'latin1',
    na.strings = '-')

head(iobt,2)
tail(iobt,2)

jj <- 2:(ncol(iobt)-1)
years <- as.integer(substr(
    colnames(iobt)[jj], 2, 5))
years

if(!file.exists("data/mInfant.rds")) {
    muns <- c("BETIM", "BELO HORIZONTE", "CONTAGEM",
              "NOVA LIMA", "RIBEIRAO DAS NEVES", "SABARA")
    (nm <- length(muns))
    (i1 <- sapply(muns, grep, nasc$Mun))
    (i2 <- sapply(muns, grep, iobt$Mun))
    (nyears <- length(years))
    
    mInfant <- data.frame(
        year = rep(years, nm),
        mun = rep(muns, each = nyears),
        nasc = as.vector(t(nasc[, jj])[, i1]),
        iobt = as.vector(t(iobt[, jj])[, i2])
    )
    saveRDS(mInfant, "data/mInfant.rds")
    
    str(mInfant)
    
    library(ggplot2)
    
    ggplot(mInfant) + theme_minimal() +
        geom_line(aes(x = year, y = iobt/nasc, color = mun))

}

## package to work work with spatial/maps

library(sf)

if(file.exists("maps/mg_mun.shp")) {
    mg_mun <- st_read("maps/mg_mun.shp")
} else {
    library(geobr)
    mg_mun <- read_municipality(
        code_muni = "MG", 
        year = 2020, 
        showProgress = TRUE
    )
    st_write(mg_mun, "maps/mg_mun.shp")
    dir("maps")
}

head(mg_mun)

o1 <- pmatch(substr(nasc$Mun, 1, 6),
             substr(mg_mun$code, 1, 6))
summary(o1)
nasc[is.na(o1),]

o2 <- pmatch(substr(iobt$Mun, 1, 6),
             substr(mg_mun$code, 1, 6))
summary(o2)
nasc[is.na(o2),]

mg_mun$nasc <- nasc$Total[o1[complete.cases(o1)]]
mg_mun$iobt <- iobt$Total[o2[complete.cases(o2)]]

bb <- st_bbox(mg_mun)
bb

rr <- c(x = diff(bb[c(1,3)]),
        y = diff(bb[c(2,4)]))
rr

rr * 200

png("iobt_map_mg.png", 2300, 1800, res = 300)
ggplot(mg_mun) + theme_minimal() +
    geom_sf(aes(fill = iobt/nasc), color = 'transparent') +
    scale_fill_distiller(
        palette = "RdBu", direction = -1
    )
dev.off()

system("eog iobt_map_mg.png &")
