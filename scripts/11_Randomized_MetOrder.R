# Generate a randomized priority list for Forward Runs
set.seed(1535)
site.name = "GOOSE"
vers=".v1"
site.lat  = 43.068496
site.lon  = -73.287425

# met.all <- read.csv(file.path("/Volumes/GoogleDrive/My Drive/PalEON_Met_Ensembles/data/", paste(site, vers, sep="."), "/aggregated/month/PDSI_AllMembers.csv"))
met.all <- read.csv("~/met_ensemble/data/met_ensembles/", paste(site, vers, sep="."), "/aggregated/month/PDSI_AllMembers.csv")
ens.names <- names(met.all[,2:ncol(met.all)])

ens.out <- data.frame(order=1:length(ens.names), ensemble=sample(ens.names, replace=F))
# write.csv(ens.out, file.path("/Volumes/GoogleDrive/My Drive/PalEON_Met_Ensembles/data/", paste(site, vers, sep="."), "/EnsembleOrder.csv"), row.names=F)
write.csv(ens.out, file.path("~/met_ensemble/data/met_ensembles/", paste(site, vers, sep="."), "/EnsembleOrder.csv"), row.names=F)

head(ens.out)
