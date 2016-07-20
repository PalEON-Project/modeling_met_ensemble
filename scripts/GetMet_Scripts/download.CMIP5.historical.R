# Script to extract a single point from a CMIP5 historical run (1850-2005)

library(httr)
library(RCurl)
library(XML)
library(lubridate)
library(ncdf4)
library(stringr)
creds <- scan("~/Desktop/ceda_creds.txt", what="character")

ceda_usr=creds[1]
ceda_pw=creds[2]

openid="https://ceda.ac.uk/openid/Christine.Rollinson"
out_base="/projectnb/dietzelab/paleon/met_ensemble/data/full_raw/bcc-csm1-1/"
cmip5_base="/badc/cmip5/data/cmip5/output1/"

model="BCC"
variant = "bcc-csm1-1"
experiment="historical"
time="day"
ensemble="r1i1p1"

path_histo_day=file.path(cmip5_base, model, variant, experiment, time, "atmos", time, ensemble, "latest/")


# Getting the base file paths from ceda
ceda_base=paste("ftp://ftp.ceda.ac.uk")
ceda_now=file.path(ceda_base, path_histo_day)

userpwd <- paste(ceda_usr, ceda_pw, sep=":")
vars <- strsplit(getURL(ceda_now, userpwd = userpwd,ftp.use.epsv = FALSE,dirlistonly = T), "\n")[[1]]
files <- strsplit(getURL(file.path(ceda_now, vars[1], ""), userpwd = userpwd,ftp.use.epsv = FALSE,dirlistonly = T), "\n")[[1]]



# Acutally accessing the files through thredds
dap_file = paste0(dap_base, "/",year, "/", doy, "/",dap.log[h,1],".ascii?")
latlon <- getURL(paste0(dap_file,"lat[0:1:599],lon[0:1:1439]"))

example <- "https://esgf1.dkrz.de/thredds/dodsC/cmip5/cmip5/output1/MPI-M/MPI-ESM-P/past1000/day/atmos/day/r1i1p1/v20111028/pr/pr_day_MPI-ESM-P_past1000_r1i1p1_18400101-18491231.nc?lat[0:1:95],lon[0:1:191]"
test <- "https://crollinson:PSU3c0l0gy@ceda.ac.uk/openid/Christine.Rollinson:esgf1.dkrz.de/thredds/dodsC/cmip5/cmip5/output1/MPI-M/MPI-ESM-P/past1000/day/atmos/day/r1i1p1/v20111028/pr/pr_day_MPI-ESM-P_past1000_r1i1p1_18400101-18491231.nc?lat[0:1:95],lon[0:1:191]"
https://esgf1.dkrz.de/thredds/dodsC/cmip5/cmip5/output1/MPI-M/MPI-ESM-P/past1000/day/atmos/day/r1i1p1/v20111028/pr/pr_day_MPI-ESM-P_past1000_r1i1p1_18400101-18491231.nc.ascii?pr[0:1:0][0:1:0][0:1:0]
part1 <- "https://crollinson:PSU3c0l0gy@ceda.ac.uk/openid/Christine.Rollinson"

test2 <- getURL(part1, curl=h)

# https://stat.ethz.ch/pipermail/r-help/2011-July/282914.html
h = getCurlHandle( cookiefile = "")
ans = getForm(openid, Username = ceda_usr, Password = ceda_pw, curl = h)

rdr = dynCurlReader(h)
ans2 = getForm(openid, Username = ceda_usr, Password = ceda_pw, curl = h, header = rdr$update)

test2 <- getURLContent(test, curl=h)
test2 <- getURLContent(test, curl=h)


library(httr)
html <- content(GET('https://jump.valueline.com/login.aspx'), "text")

viewstate <- as.character(sub('.*id="_VIEWSTATE" value="([0-9a-zA-Z+/=]*).*', '\\1', html))

params <- list(
  'ct100$ContentPlaceHolder$LoginControl$txtUserID' = 'MY USERNAME',
  'ct100$ContentPlaceHolder$LoginControl$txtUserPw' = 'MY PASSWORD',
  'ct100$ContentPlaceHolder$LoginControl$btnLogin' = 'Sign In',
  '_VIEWSTATE' = viewstate
)
POST('https://jump.valueline.com/login.aspx', body = params)

https://ceda.ac.uk/OpenID/Provider/server?openid.ns=http%3A%2F%2Fspecs.openid.net%2Fauth%2F2.0&openid.claimed_id=http%3A%2F%2Fspecs.openid.net%2Fauth%2F2.0%2Fidentifier_select&openid.identity=http%3A%2F%2Fspecs.openid.net%2Fauth%2F2.0%2Fidentifier_select&openid.return_to=https%3A%2F%2Fesgf1.dkrz.de%2Fesg-orp%2Fj_spring_openid_security_check.htm&openid.realm=https%3A%2F%2Fesgf1.dkrz.de%2F&openid.assoc_handle={HMAC-SHA256}{577eb2ea}{1uBpAQ%3D%3D}&openid.mode=checkid_setup&openid.ns.ext1=http%3A%2F%2Fopenid.net%2Fsrv%2Fax%2F1.0&openid.ext1.mode=fetch_request&openid.ext1.type.email=http%3A%2F%2Fopenid.net%2Fschema%2Fcontact%2Finternet%2Femail&openid.ext1.type.firstname=http%3A%2F%2Fopenid.net%2Fschema%2FnamePerson%2Ffirst&openid.ext1.type.lastname=http%3A%2F%2Fopenid.net%2Fschema%2FnamePerson%2Flast&openid.ext1.required=email%2Cfirstname%2Clastname

latlon <- getURL(example)
lat.ind <- gregexpr("lat", latlon)
lon.ind <- gregexpr("lon", latlon)
lats <- as.vector(read.table(con <- textConnection(substr(latlon, lat.ind[[1]][3], lon.ind[[1]][3]-1)), sep=",", fileEncoding="\n", skip=1))
lons <- as.vector(read.table(con <- textConnection(substr(latlon, lon.ind[[1]][3], nchar(latlon))), sep=",", fileEncoding="\n", skip=1))

lat.use <- which(lats-0.25/2<=lat.in & lats+0.25/2>=lat.in)
lon.use <- which(lons-0.25/2<=lon.in & lons+0.25/2>=lon.in)




# n 6/25/2011 6:16 PM, Duncan Temple Lang wrote:
#   > Hi Steve
# >
#   >   RCurl can help you when you need to have more control over Web requests.
# > The details vary from Web site to Web site and the different ways to specify
# > passwords, etc.
# >
#   > If the JSESSIONID and NCES_JSESSIONID are regular cookies and returned in the first
# > request as cookies, then you can just have RCurl handle the cookies
# > But the basics for your case are
# >
#   >    library(RCurl)
# >    h = getCurlHandle( cookiefile = "")
# >
#   > Then make your Web request using getURLContent(), getForm() or postForm()
# > but making certain to pass the curl handle  stored in h in each call, e.g.
# >
#   >    ans = getForm(yourURL, login = "bob", password = "jane", curl = h)
# >
#   >    txt = getURLContent(dataURL, curl = h)
# >
#   >
#   > If JSESSIONID and NCES_JSESSIONID are not returned as cookies but HTTP header fields, then you
# > need to process the header.
# > Something like
# >
#   >    rdr = dynCurlReader(h)
# >
#   >    ans = getForm(yourURL, login = "bob", password = "jane", curl = h, header = rdr$update)
# >
#   > Then the header  from the HTTP response is available as
# >    rdr$header()
# >
#   > and you can use parseHTTPHeader(rdr$header()) to convert it into a named vector.
# >
