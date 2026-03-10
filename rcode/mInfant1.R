
mInfant <- readRDS("data/mInfant.rds")

head(mInfant)

library(ggplot2)

ggplot(mInfant) + theme_minimal() +
    geom_line(aes(x = year, y = iobt/nasc, color = mun))


## lambda_0 = sum_i y_i / sum_i N_i
lambda0 <- sum(mInfant$iobt) / sum(mInfant$nasc)
lambda0

## E_i = lambda_0 * N_i

mInfant$Ei <- lambda0 * mInfant$nasc

sum(mInfant$iobt)
sum(mInfant$Ei)

muns <- sort(unique(mInfant$mun))
muns

years <- sort(unique(mInfant$year))
years

mInfant$t <- mInfant$year -years[1] + 1
mInfant$m <- pmatch(mInfant$mun, muns, duplicates.ok = TRUE)

head(mInfant)

R <- crossprod(diff(diag(4)))
R

kronecker(diag(3), R)

f1 <- iobt ~ 1 + f(t, model = 'rw1', replicate = m,
                   scale.model = TRUE) ## tau marginal

library(INLA)

fit1 <- inla(
    formula = f1,
    family = 'poisson',
    data = mInfant,
    E = Ei,
    control.compute = list(waic = TRUE)
)

fit1$summary.fixed

fit1$summary.hy

1/sqrt(53)
1/sqrt(12)

plot(fit1)

## RW2
f2 <- iobt ~ 1 + f(t, model = 'rw2', replicate = m,
                   scale.model = TRUE) ## tau marginal

fit2 <- inla(
    formula = f2,
    family = 'poisson',
    data = mInfant,
    E = Ei,
    control.compute = list(waic = TRUE)
)

fit2$summary.fixed

fit2$summary.hy

1/sqrt(2325)
1/sqrt(56)

plot(fit2)

c(m1 = fit1$waic$waic,
  m2 = fit2$waic$waic)
