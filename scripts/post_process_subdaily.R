ensemble.post.process <- function(sims.out, site.name, GCM, )

for(e in ens.day){
  # # Do the prediction
  # sims.out <- predict.subdaily(dat.mod=dat.ens[[paste0("X", e)]], n.ens=ens.hr, path.model=file.path(dat.base, "subday_models"), lags.init=lags.init[[paste0("X", e)]], dat.train=dat.train)
  
  
  # If this is one of our designated QAQC years, makes some graphs
  if(y %in% yrs.plot){
    day.name <- paste0(site.name, "_", GCM, "_1hr_", str_pad(e, 3, pad=0))
    fig.ens <- file.path(path.out, "subdaily_qaqc", day.name)
    if(!dir.exists(fig.ens)) dir.create(fig.ens, recursive=T)
    
    for(v in names(sims.out)){
      graph.predict(dat.mod=dat.ens[[paste0("X", e)]], dat.ens=sims.out, var=v, fig.dir=fig.ens)
    }
  }
  
  
  # Update the initial lags for next year
  for(v in names(sims.out)){
    lags.init[[paste0("X",e)]][[v]] <- data.frame(sims.out[[v]][length(sims.out[[v]]),])
  }
  
  # -----------------------------------
  # Write each year for each ensemble member into its own .nc file
  # -----------------------------------
  for(i in 1:ens.hr){
    ens.name <- paste0(site.name, "_", GCM, "_1hr_", str_pad(e, 3, pad=0), "-", str_pad(i, 3, pad=0))
    
    if(!dir.exists(file.path(path.out, ens.name))) dir.create(file.path(path.out, ens.name), recursive=T)
    path.out <- file.path(dat.base, GCM, "1hr")
    
    var.list <- list()
    dat.list <- list()
    for(v in names(sims.out)){
      var.cf = vars.info[vars.info$name==v, "name.cf"]
      var.list[[v]] <- ncvar_def(v, units=paste(vars.info[vars.info$name==v, "units"]), dim=list(dimX, dimY, dim.t), longname=paste(vars.info[vars.info$name==v, "longname"]))
      dat.list[[v]] <- array(sims.out[[v]][,i], dim=c(1,1,length(hrs.now)))
    }
    
    # Naming convention: [SITE]_[GCM]_1hr_[bias_ens_member]-[subday_ens_member]_[YEAR].nc
    nc <- nc_create(file.path(path.out, ens.name, paste0(ens.name,"_", str_pad(y, 4, pad=0), ".nc")), var.list)
    for(v in 1:length(var.list)) {
      ncvar_put(nc, var.list[[v]], dat.list[[v]])
    }
    nc_close(nc)    
  }
  # -----------------------------------
  
} # End ensemble member prediction for 1 year