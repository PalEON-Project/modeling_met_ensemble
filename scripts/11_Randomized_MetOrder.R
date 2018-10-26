# Generate a randomized priority list for Forward Runs
set.seed(1535)
site.name = "PALMGHATT"
vers=".v1"
site.lat  = 41.703568
site.lon  = -74.253876

met.all <- read.csv(file.path("/Volumes/GoogleDrive/My Drive/PalEON_Met_Ensembles/data/", paste0(site.name, vers), "/aggregated/month/PDSI_AllMembers.csv"))
# met.all <- read.csv("~/met_ensemble/data/met_ensembles/", paste0(site.name, vers), "/aggregated/month/PDSI_AllMembers.csv")
ens.names <- names(met.all[,2:ncol(met.all)])

ens.out <- data.frame(order=1:length(ens.names), ensemble=sample(ens.names, replace=F))
write.csv(ens.out, file.path("/Volumes/GoogleDrive/My Drive/PalEON_Met_Ensembles/data", paste0(site.name, vers), "EnsembleOrder.csv"), row.names=F)
# write.csv(ens.out, file.path("~/met_ensemble/data/met_ensembles", paste(site.name, vers), "EnsembleOrder.csv"), row.names=F)

head(ens.out)
