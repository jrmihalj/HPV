## Raw data analysis and visualization 
## 12.2.2016 
## Sylvia Ranjeva 
## -----------------------------------------------------

###-------------------------------------------------------------------------------------------------------------------------------------------
library(reshape)
library(RSQLite)
library(tidyr)
library(ggplot2)
library(dplyr)
library(directlabels)
library(cooccur)
theme_set(theme_bw(base_size = 14))

## Step 1: Extract infection data and visit dates ---------------------------------------------
dbFilename <- "./data/HIMdata_Mar_2015.sqlite"
db <- dbConnect(SQLite(), dbFilename)
tableNames <- dbListTables(db)
inf_status <- dbReadTable(db,"infectionStatusByVisit")
visit_dates <- dbReadTable(db,"visitDates")
dbDisconnect(db)

# vis_dates <- reshape(visit_dates, idvar = "subjectId",  timevar = "visitId", direction = "wide")
# vis_dates2 <- apply(vis_dates[,2:11],2, as.Date, origin = '1970-1-1')
# vis_dates3 <- as.data.frame(apply(vis_dates2,2, as.numeric))
# vis_dates3$subjectId <- vis_dates$subjectId
HPV_types <- paste0("hpv",c(6,11,16,18,31,33,45,52,58,84))
n_test_strains <- length(HPV_types)


## Step2 ## Examine  patterns of cooccurance ---------------------------------------------------------------

get_cooccurance <- function(data_subset, vis_id, threshold = TRUE){
  inf_stat <- data_subset[data_subset$visitId == vis_id,]
  inf_stat$visitId <- NULL
  inf_v <-reshape(inf_stat, idvar = "subjectId", timevar = "strainId", direction = "wide")
  inf_v$subjectId <- NULL
  strainIds <- unique(data_subset$strainId)
  names(inf_v) <- strainIds
  inf_v <- as.data.frame(apply(inf_v,2,as.numeric))
  inf_v <- inf_v[!is.na(rowSums(inf_v)),]
  inf <- t(inf_v)
  cooccur.hpv <- cooccur(mat = inf, type = "spp_site", thresh = threshold, spp_names = TRUE)
  return(cooccur.hpv)
}

full_data_cooccur <- get_cooccurance(data_subset = inf_status[inf_status$strainId %in% HPV_types,], vis_id = "v1", threshold = TRUE)
pdf("./figs/cooccurance_patterns.pdf")
plot(full_data_cooccur)
dev.off()
pdf("./figs/coocurrance_trends.pdf")
obs.v.exp(mod = full_data_cooccur)
dev.off()
pdf("./figs/profile_cooccurance.pdf")
pair.profile(mod = full_data_cooccur)
dev.off()


## Step3 ## Examine prevalence of types across visits and across time ----------------------------
## Visit level prevalences #######################################################
visit_prev <- data.frame()
strains <- HPV_types
for( i in 1:length(strains)){
  this_strain <- strains[i]
  inf_df <- inf_status[inf_status$strainId == this_strain,]
  inf_df$strainId <- NULL
  inf <- reshape(inf_df, idvar = "subjectId", timevar = "visitId", direction = "wide" )
  prev <- as.numeric(colMeans(inf[names(inf)!="subjectId"], na.rm=T))
  visit_prev <- rbind(visit_prev,prev)
  visit_prev$type[i] <- this_strain
}

names(visit_prev) <- c(paste0("v", c(1:10)), "strain")
dfm <- melt(visit_prev, id.vars = "strain")

# Calculate average prevalence per type per visit 

avg_prev <- visit_prev %>%
  mutate(mean = select(., starts_with("v")) %>%
           rowMeans(na.rm = TRUE)) %>% select(.,c(strain,mean))

## Plot the prevalence ## 
pdf("./figs/visit_level_prev.pdf")
p <- ggplot(dfm, aes(x = variable, y = value, group = strain)) + geom_line(aes(color = strain)) + geom_text(aes(label = toupper(strain)))
p + xlab("") + ylab("Prevalence") 
dev.off()
  
## Get prevalence in real time ## ---------------------------------------------------------------------------------
dfm_dates <- visit_dates #melt(visit_dates, id.vars = "subjectId")
names(dfm_dates) <- c("subjectId", "visitId","date")
min_date <- min(as.numeric(as.Date(dfm_dates$date)), na.rm=T)
max_date <- max(as.numeric(as.Date(dfm_dates$date)), na.rm=T)

date_seq <- seq(to = max_date, from = min_date, by = .5*365)
get_prev <- function(type, st){
  inf <- inf_status %>% filter(strainId == type) %>% 
    select(-strainId)
  dfm <- merge(inf, dfm_dates, by = c("subjectId","visitId"))
  dfm$date <- as.numeric(as.Date(dfm$date))
  prev <- array(NA, length(st))
  for( i in 2:length(st)){
    dfm %>% filter(date > st[i-1] & date < st[i] & !is.na(status)) -> dfm_sub
    prev[i] <- sum(dfm_sub$status)/nrow(dfm_sub)
  }
  prev <- prev[-1]
  return(prev)
}

df_all <- data.frame()

for ( i in 1:length(HPV_types)){
  df <- data.frame(date = date_seq[-1],
                   type = HPV_types[i],
                   prev = get_prev(HPV_types[i], date_seq)
  )
  df_all <- rbind(df_all,df)
}
pdf("./figs/true_prevalence.pdf")
df_all$type = toupper(df_all$type)
p <- ggplot(df_all, aes(x = as.Date(date, origin = '1970-1-1'), y = prev, group = type)) + geom_line(aes(color = type)) 
p + xlab("") + ylab("Prevalence") 
dev.off()

# Get mean prevalence  ---------------------------------------------
prev <- reshape(df_all, timevar = "date", idvar = "type", direction = "wide")
avg_prev <- prev %>%
  mutate(mean = select(., starts_with("prev")) %>%
           rowMeans(na.rm = TRUE)) %>% select(.,c(type,mean))
avg_prev[order(avg_prev$mean),]

