# Generate a randomized priority list for Forward Runs
set.seed(1535)
site <- "TENSIONZONE"
vers <- "v1"

met.all <- read.csv("~/Google Drive/PalEON_Met_Ensembles/data/", paste(site, vers, sep="."), "/month/PDSI_AllMembers.csv")
ens.names <- names(met.all[,2:ncol(met.all)])

ens.out <- data.frame(order=1:length(ens.names), ensemble=sample(ens.names, replace=F))
write.csv(ens.out, "~/Google Drive/PalEON_Met_Ensembles/data/", paste(site, vers, sep="."), "/EnsembleOrder.csv", row.names=F)
