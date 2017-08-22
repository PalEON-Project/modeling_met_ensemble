library(ncdf4)

mod.dir <- "~/Desktop/phase3_met_ensembe/bias_correct_day/HARVARD/MIROC-ESM/day/HARVARD_MIROC-ESM_day_001"

files.mod <- dir(mod.dir)

ncT <- nc_open(file.path(mod.dir, files.mod[1]))

summary(ncT$var)

tmax <- ncvar_get(ncT, "tmax")
summary(tmax)
plot(tmax[1:(365*5)], type="b")
plot(tmax, type="l")

