### Take data from the HIM dataset and get it in form for the hybrid metacom model ###-----------------------------------------------------
### Sylvia Ranjeva
### April 2016
###-------------------------------------------------------------------------------------------------------------------------------------------
library(reshape)
library(RSQLite)
use_complete_data = TRUE
## FUNCTIONS ################## ----------------------------------------------------------
## stupid function to convert visit dates to numeric times ###
get_timestep <- function(t){
 return(as.numeric(strsplit(t,"v")[[1]][2]))
}

## function to get infection data for one strain 
get_data_for_strain <- function(i, inf = inf_status, strain_list = test_strains){
  this_strain = strain_list[[i]]
  cat("strain is", this_strain, "\n")
  data <- inf[inf$strainId == this_strain & inf$subjectId %in% test_pat,]
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


## Construct matrix of time between visits # -----------
tbv <- matrix(NA, n_pat, n_vis)
tbv[,1] <- array(0, nrow(tbv))
for( i in 2:n_vis){
  tbv[,i] <- vis_dates3[,i] - vis_dates3[,(i-1)]
}

tbv_df <- data.frame(subjectId = vis_dates3$subjectId,
                     tbv)
names(tbv_df) <- c("subjectId",paste0("v",c(1:n_vis)))

##  Specify strains/patients for analysis ------------------------------------------------------
strainIds <- unique(inf_status$strainId)
subjectIds <- unique(inf_status$subjectId)
strainIds <- unique(inf_status$strainId)
n_strains =  length(strainIds)
n_pat <- length(subjectIds)
visitIds <- unique(inf_status$visitId)
n_vis = length(visitIds)

test_pat <- sample(subjectIds, size = n_pat, replace = FALSE)
test_strains <- paste0("hpv",c(6,11,16,18)) #,31,33,45,52,58,84))
n_test_strains <- length(test_strains)


## Get the data and format it # ------------------------------------------------------------
df1 <- get_data_for_strain(1)
df2 <- get_data_for_strain(2)
df_all <- merge(df1,df2, by = c("site","time"))
if(n_test_strains > 2){
  for( i in 3:n_test_strains){
    df <- get_data_for_strain(i)
    df_all <- merge(df_all, df, by = c("site","time"))
  }
}

df_all <- df_all[order( df_all$time),]
names(df_all) <- c("patient", "visit", test_strains)

# get the correct set of time between visit data according to test patients --------------
tbv <- melt(subset(tbv_df, subjectId %in% df_all$patient), id.vars = "subjectId")

## Format all data into stan code input --------------------------------------------------
n_visits <- length(visitIds)
n_patients <- n_pat
n <- n_patients * n_visits
y = as.matrix(df_all[,names(df_all) %in% test_strains])
patient <- rep(c(1:n_patients), times = n_visits)
visit <- df_all$visit
time_between_visits <- tbv$value
d <- ncol(y)

if(!use_complete_data){
  missing <- which(is.na(y[,1]))
  visit[missing] <- -1
}

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
               tbv = time_between_visits,
               dir_prior = c(1, 1))

save(stan_d, subjectIds, file = "test_data_HIM_full.rda")

## NOT including 139,55,64###
