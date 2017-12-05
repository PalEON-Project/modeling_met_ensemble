##' Aggregate met 
# ----------------------------------- 
# Description
# -----------------------------------
##' @title aggregate.met
##' @family Meteorology
##' @author Christy Rollinson
##' @description Take sub-daily met data and turn it into daily &/or monthly data; saving by year as we go
# ----------------------------------- 
# Parameters
# -----------------------------------
##' @param path.in File path to where the .nc files should be housed
##' @param years.agg number of betas to pull from; if NULL (default) does ALL files
##' @param save.day [TRUE]/FALSE, whether to save daily data
##' @param save.month [TRUE]/FALSE whether to calculate/save daily data
##' @param out.base root file path for where you want to write data
##' @param day.dir subdirectory you want to save daily output in (defaults to "day")
##' @param mo.dir subdirectory you want to save monthly output in (defaults to "month")
##' @param add.vars c("daylength", "temp.max", "temp.min") additional handy variables that we can compute other than just mean
##' @param parallel [TRUE]/FALSE whether to perform aggregation in parrallel (years in parallel)
##' @param n.cores number of cores to use for parallel processing
##' @param print.progress [TRUE]/FALSE print a progress bar when not using parallel processing?
##' @param verbose spit out extra information with ncdf4
##' @return Automatically saved files withthe same naming structure as previously, just in a different location (this is dangerous!)
##' @export
# -----------------------------------
#----------------------------------------------------------------------
# Begin Function
#----------------------------------------------------------------------
aggregate.met <- function(path.in, years.agg=NULL, save.day=T, save.month=T, 
                          out.base, day.dir="day", mo.dir="month", 
                          add.vars=c("daylength", "air_temperature_maximum", "air_temperature_minimum"),
                          parallel=F, n.cores=NULL, 
                          print.progress=T, verbose=FALSE){
  
  # Create our output folders
  if(save.day==TRUE) dir.create(file.path(out.base, day.dir), recursive=T, showWarnings = FALSE)
  if(save.month==TRUE) dir.create(file.path(out.base, mo.dir), recursive=T, showWarnings = FALSE)
  
  # Get a list of all the files available to do the aggregation for
  files.agg <- dir(path.in)
  yrs.files <- strsplit(files.agg, "[.]")
  yrs.files <- matrix(unlist(yrs.files), ncol=length(yrs.files[[1]]), byrow=T)
  yrs.files <- as.numeric(yrs.files[,ncol(yrs.files)-1])
  
  # Subset to just the files we want to work with
  if(is.null(years.agg)){
    years.agg <- yrs.files
  }
  files.agg <- files.agg[which(yrs.files %in% years.agg)]
  yrs.agg <- yrs.files[which(yrs.files %in% years.agg)]
  
  out <- list()
  if(parallel==T){
    files.agg <- as.list(files.agg)
    
    out <- parallel::mclapply(files.agg, aggregate.file, 
                              mc.cores=min(n.cores, length(files.agg)), mc.silent=TRUE,
                              path.in=path.in, 
                              save.day=save.day, save.month=save.month, 
                              out.base=out.base, day.dir=day.dir, mo.dir=mo.dir, 
                              add.vars=add.vars, verbose=verbose)
  } else {
    if(print.progress==T) pb <- txtProgressBar(min=0, max=length(files.agg), style=3)
    for(i in 1:length(files.agg)){
      out[[i]] <- aggregate.file(f.agg=files.agg[i], path.in=path.in, 
                                 save.day=save.day, save.month=save.month, 
                                 out.base=out.base, day.dir=day.dir, mo.dir=mo.dir, 
                                 add.vars=add.vars, verbose=verbose)
      if(print.progress==T) setTxtProgressBar(pb, i)
    }
  }

}