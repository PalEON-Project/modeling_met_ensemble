# Generate a randomized priority list for Forward Runs
set.seed(1535)
site.name = "GLSP"
vers=".v1"
site.lat  = 45.54127
site.lon  = -95.5313

# met.all <- read.csv(file.path("/Volumes/GoogleDrive/My Drive/PalEON_Met_Ensembles/data/", paste(site, vers, sep="."), "/aggregated/month/PDSI_AllMembers.csv"))
met.all <- read.csv("~/met_ensemble/data/met_ensembles/", paste(site, vers, sep="."), "/aggregated/month/PDSI_AllMembers.csv")
ens.names <- names(met.all[,2:ncol(met.all)])

ens.out <- data.frame(order=1:length(ens.names), ensemble=sample(ens.names, replace=F))
# write.csv(ens.out, file.path("/Volumes/GoogleDrive/My Drive/PalEON_Met_Ensembles/data/", paste(site, vers, sep="."), "/EnsembleOrder.csv"), row.names=F)
write.csv(ens.out, file.path("~/met_ensemble/data/met_ensembles/", paste(site, vers, sep="."), "/EnsembleOrder.csv"), row.names=F)

head(ens.out)
