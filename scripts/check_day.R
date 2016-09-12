library(ncdf4)

mod.dir <- "~/Desktop/phase3_met_ensembe/bias_correct_day/HARVARD/MIROC-ESM/day/HARVARD_MIROC-ESM_day_001"

files.mod <- dir(mod.dir)

ncT <- nc_open(file.path(mod.dir, files.mod[1]))

summary(ncT$var)

tmax <- ncvar_get(ncT, "tmax")
summary(tmax)
plot(tmax[1:365], type="b")

load('~/Desktop/phase3_met_ensembe/bias_correct_day/HARVARD/MIROC-ESM/day/MIROC-ESM_day_alldata.Rdata')


summary(dat.out.full)
summary(dat.out.full$tmax$sims)

# Sort the data
dat.out.full$tmax$sims <- dat.out.full$tmax$sims[order(dat.out.full$tmax$sims$time),] 

dat.out.full$tmax$sims[1:100,1:10]

plot(dat.out.full$tmax$sims[1:700,"X1"], type="b")