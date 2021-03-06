### Take data from the HIM dataset and get it in form for the hybrid metacom model ###-----------------------------------------------------
### Sylvia Ranjeva
### April 2016
###-------------------------------------------------------------------------------------------------------------------------------------------
library(reshape)
library(RSQLite)
library(tidyr)
use_complete_data = FALSE
## FUNCTIONS ################## ----------------------------------------------------------
## stupid function to convert visit dates to numeric times ###
get_timestep <- function(t){
 return(as.numeric(strsplit(t,"v")[[1]][2]))
}

## function to get infection data for one strain 
get_data_for_strain <- function(i, inf = inf_status, strain_list = test_strains, test_pat = subjectIds){
  this_strain = strain_list[[i]]
  cat("strain is", this_strain, "\n")
  data <- inf[inf$strainId == this_strain & inf$subjectId %in% test_pat,]
  data$time <- sapply(data$visitId, get_timestep)
  data <- data[order(data$subjectId),]
  site <- as.numeric(data$subjectId)
  time <- as.numeric(data$time)
  y <- as.numeric(data$status)
  this_strain_df <- data.frame(site = site, time = time, y = y)
  return(this_strain_df)
}

## Create data for stan model #####-----------------------------------------------------------------
dbFilename <- "data/HIMdata_Mar_2015.sqlite"
db <- dbConnect(SQLite(), dbFilename)
tableNames <- dbListTables(db)
inf_status <- dbReadTable(db,"infectionStatusByVisit")
visit_dates <- dbReadTable(db,"visitDates")
dbDisconnect(db)

vis_dates <- reshape(visit_dates, idvar = "subjectId",  timevar = "visitId", direction = "wide")
vis_dates2 <- apply(vis_dates[,2:11],2, as.Date, origin = '1970-1-1')
vis_dates3 <- as.data.frame(apply(vis_dates2,2, as.numeric))
vis_dates3$subjectId <- vis_dates$subjectId

if(use_complete_data){
  complete_data_indices <- which(!is.na(rowSums(vis_dates3)))
  complete_data_ids <- vis_dates$subjectId[complete_data_indices]
  inf_status_complete_visits <- subset(inf_status, subjectId %in% complete_data_ids)
  missing_data_indices_2 <- which(is.na(inf_status_complete_visits$status))
  problem_patients <- unique(inf_status_complete_visits[missing_data_indices_2,]$subjectId)
  inf_status_complete <- subset(inf_status_complete_visits, !(subjectId %in% problem_patients))
  inf_status <- inf_status_complete
  vis_dates3 <- subset(vis_dates3, subjectId %in% inf_status$subjectId)
}

##  Specify strains/patients for analysis ------------------------------------------------------
strainIds <- unique(inf_status$strainId)
subjectIds <- unique(inf_status$subjectId)
strainIds <- unique(inf_status$strainId)
n_strains =  length(strainIds)
n_pat <- length(subjectIds)
visitIds <- unique(inf_status$visitId)
n_vis = length(visitIds)

test_pat <- sample(subjectIds, size = n_pat, replace = FALSE)
test_strains <- paste0("hpv",c(6,11,16,18,31,33,45,52,58,84))
n_test_strains <- length(test_strains)

## Construct matrix of time between visits # -----------
tbv <- matrix(NA, nrow(vis_dates3), n_vis)
tbv[,1] <- array(0, nrow(tbv))
for( i in 2:n_vis){
  tbv[,i] <- vis_dates3[,i] - vis_dates3[,(i-1)]
}

tbv_df <- data.frame(subjectId = vis_dates3$subjectId,
                     tbv)
names(tbv_df) <- c("subjectId",paste0("v",c(1:n_vis)))


## Get the data and format it # ------------------------------------------------------------
df1 <- get_data_for_strain(1, test_pat = test_pat)
df2 <- get_data_for_strain(2, test_pat = test_pat )
df_all <- merge(df1,df2, by = c("site","time"))
if(n_test_strains > 2){
  for( i in 3:n_test_strains){
    df <- get_data_for_strain(i, test_pat = test_pat)
    df_all <- merge(df_all, df, by = c("site","time"))
  }
}

df_all <- df_all[order( df_all$site, df_all$time),]
names(df_all) <- c("patient", "visit", test_strains)

patients <- data.frame(subjectId = unique(df_all$patient),
                       patient = c(1:length(unique(df_all$patient)))
                    )
df_all$pat <- 0
for( i in 1:nrow(df_all)){
  cat(" i is: ", i , "\n")
  df_all[i,]$pat<- patients[which(df_all[i,]$patient == patients$subjectId ),]$patient
}

# get the correct set of time between visit data according to test patients --------------
tbv <- melt(subset(tbv_df, subjectId %in% df_all$patient), id.vars = "subjectId")
tbv <- tbv[order(tbv$subjectId),]


## Format all data into stan code input --------------------------------------------------
## if variable numbers of visits are allowed, remove missing data 
if(!use_complete_data){
  missing_genotype <- which(is.na(rowSums(df_all[,-c(1:2)])))
  missing_tbv <- which(is.na(tbv$value))
  missing_any <- union(missing_genotype, missing_tbv)
  df_all <- df_all[-missing_any,]
  tbv <- tbv[-missing_any,]
}

stopifnot(tbv$subjectId == df_all$patient)
n_visits <- length(visitIds)
n_patients <- n_pat
n <- nrow(df_all)
y = as.matrix(df_all[,names(df_all) %in% test_strains])
patient <- df_all$pat
visit <- df_all$visit
time_between_visits <- tbv$value
d <- ncol(y)

subjectIds <- data.frame(subjectId = df_all$patient,
                           patient_num = patient)

stan_d <- list(n = n, 
               n_strains = d, 
               y = y,
               eta = 2,
               patient = patient, 
               n_patient = n_patients,
               n_visit_max = n_visits, 
               visit = visit,
               tbv = time_between_visits,
               dir_prior = c(1, 1))

save(stan_d, subjectIds, file = "data/test_data_HIM_full_10_strains.rda")

## NOT including 139,55,64###
