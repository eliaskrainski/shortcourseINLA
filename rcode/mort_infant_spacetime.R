library(ggplot2)    

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

(jj <- 2:(ncol(iobt)-1))
years <- as.integer(substr(
    colnames(iobt)[jj], 2, 5))
years

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

names(nasc)


## data frame
jj2 <- 21:30
nt <- length(jj2)
nt
n_areas <- nrow(mg_mun)
n_areas

idat <- data.frame(
    area = rep(1:n_areas, nt),
    tempo = rep(1:nt, each = n_areas),
    N = unlist(nasc[o1[complete.cases(o1)], jj2]),
    y = unlist(iobt[o2[complete.cases(o2)], jj2])
)

lambda0 <- sum(idat$y, na.rm = TRUE) / sum(idat$N)
lambda0

idat$E <- lambda0 * idat$N

head(idat)

library(spdep)

viz <- poly2nb(mg_mun)

viz

xy <- st_coordinates(st_centroid(mg_mun))

par(mar = c(0,0,0,0))
plot(viz, xy)

nb2INLA(file = "grafo", nb = viz) ## 

### replicar o modelo 'besagproper' para areas
ff0_st <- y ~ f(area, model = 'besagproper', graph = 'grafo',
                replicate = tempo)

### combinar o modelo 'besagproper' para ares
### com o modelo 'ar1' para o tempo
ff1_st <- y ~ f(area, model = 'besagproper', graph = 'grafo',
                group = tempo, control.group = list(model = 'ar1')) 

library(INLA)

fit0 <- inla(
    formula = ff0_st,
    family = "poisson",
    data = idat, 
    E = E,
    control.predictor = list(link = 1)
)


fit1 <- inla(
    formula = ff1_st,
    family = "poisson",
    data = idat, 
    E = E,
    control.predictor = list(link = 1)
)

fit0$summary.fixed
fit1$summary.fixed

fit0$summary.hyperpar

fit1$summary.hyperpar


gcv0 <- inla.group.cv(fit0, num.level.sets = 3)
gcv1 <- inla.group.cv(fit1, num.level.sets = 3)

str(gcv0,1)

gcv0$groups[1:2]

mean(log(gcv0$cv), na.rm = TRUE)
mean(log(gcv1$cv), na.rm = TRUE)


