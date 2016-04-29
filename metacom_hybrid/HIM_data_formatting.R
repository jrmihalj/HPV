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
dbDisconnect(db)
inf_status_complete <- inf_status[!is.na(inf_status$status),]
# who has all ten visits
visitIds <- unique(inf_status_complete$visitId)
strainIds <- unique(inf_status_complete$strainId)

n_vis = length(visitIds)
n_strains =  length(strainIds)

# get people with data from all visits 
n_obs_per_pat <- table(inf_status_complete$subjectId)
complete_obs_pat <- as.numeric(names(n_obs_per_pat[which(n_obs_per_pat == n_strains * n_vis)]))
inf_status_complete <- inf_status_complete[inf_status_complete$subjectId %in% complete_obs_pat,]

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

stan_d <- list(n_unit = n_pat,
               n_time = n_vis,
               n = n_pat * n_vis,
               m = length(test_strains),
               unit = rep(c(1:n_pat), times = n_vis),
               time = df_all$time,
               y = as.matrix(df_all[,3:ncol(df_all)]),
               y_mat = as.matrix(df_all[,3:ncol(df_all)])
               )
save(stan_d, file = "test_data_HIM.rda")


