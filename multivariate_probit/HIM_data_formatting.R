### Take data from the HIM dataset and get it in form for the hybrid metacom model ###-----------------------------------------------------
### Sylvia Ranjeva
### April 2016
###-------------------------------------------------------------------------------------------------------------------------------------------
library(reshape)
library(RSQLite)

## FUNCTIONS ################## ----------------------------------------------------------
## stupid function to convert visit dates to numeric times ###
get_timestep <- function(t){
 return(as.numeric(strsplit(t,"v")[[1]][2]))
}

## function to get infection data for one strain 
get_data_for_strain <- function(i, strain_list = test_strains){
  this_strain = strain_list[[i]]
  cat("strain is", this_strain, "\n")
  data <- inf_status_complete[inf_status_complete$strainId == this_strain & inf_status_complete$subjectId %in% test_pat,]
  data$time <- sapply(data$visitId, get_timestep)
  data <- data[order(data$time),]
  site <- as.numeric(data$subjectId)
  time <- as.numeric(data$time)
  y <- as.numeric(data$status)
  this_strain_df <- data.frame(site = site, time = time, y = y)
  return(this_strain_df)
}

## Create data for stan model #####-----------------------------------------------------------------
dbFilename <- "HIMdata_Mar_2015.sqlite"
db <- dbConnect(SQLite(), dbFilename)
tableNames <- dbListTables(db)
inf_status <- dbReadTable(db,"infectionStatusByVisit")
visit_dates <- dbReadTable(db,"visitDates")
dbDisconnect(db)


vis_dates <- reshape(visit_dates, idvar = "subjectId",  timevar = "visitId", direction = "wide")
vis_dates2 <- apply(vis_dates[,2:11],2, as.Date, origin = '1970-1-1')
vis_dates3 <- apply(vis_dates2,2, as.numeric)
complete_data_indices <- which(!is.na(rowSums(vis_dates3)))
complete_data_ids <- vis_dates$subjectId[complete_data_indices]
inf_status_complete_visits <- subset(inf_status, subjectId %in% complete_data_ids)
missing_data_indices_2 <- which(is.na(inf_status_complete_visits$status))
problem_patients <- unique(inf_status_complete_visits[missing_data_indices_2,]$subjectId)
inf_status_complete <- subset(inf_status_complete_visits, !(subjectId %in% problem_patients))


visitIds <- unique(inf_status$visitId)
strainIds <- unique(inf_status$strainId)
n_vis = length(visitIds)
n_strains =  length(strainIds)

# get people with data from all visits 
#n_obs_per_pat <- table(inf_status_complete$subjectId)
#complete_obs_pat <- as.numeric(names(n_obs_per_pat[which(n_obs_per_pat == n_strains * n_vis)]))
#inf_status_complete <- inf_status_complete[inf_status_complete$subjectId %in% complete_obs_pat,]

strainIds <- unique(inf_status_complete$strainId)
subjectIds <- unique(inf_status_complete$subjectId)
strainIds <- unique(inf_status_complete$strainId)
subjectIds <- unique(inf_status_complete$subjectId)
n_pat <- 200
test_pat <- sample(subjectIds, size = n_pat, replace = FALSE)
## choose two HPV types to test method
test_strains <- as.list(strainIds[c(1,3)])

## TODO: set this up with lapply and merge data frame results
#df_list <- lapply(c(1:length(test_strains)),get_data_for_strain)
df1 <- get_data_for_strain(1)
df2 <- get_data_for_strain(2)
df_all <- merge(df1,df2, by = c("site","time"))
df_all <- df_all[order( df_all$time),]
names(df_all) <- c("patient", "visit", test_strains)


## Format data for stan code 
n_visits <- length(visitIds)
n_patients <- n_pat
n <- n_patients * n_visits
y = as.matrix(df_all[,-c(1:2)])
patient <- rep(c(1:n_patients), times = n_visits)
visit <- df_all$visit
d <- ncol(y)

subjectIds <- data.frame(subjectId = df_all$patient,
                           patient_num = patient)

stan_d <- list(n = n, 
               d = d, 
               y = y,
               eta = 2,
               patient = patient, 
               n_patient = n_patients,
               n_visit = n_visits, 
               visit = visit,
               dir_prior = c(1, 1))

save(stan_d, subjectIds, file = "test_data_HIM.rda")


