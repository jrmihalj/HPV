# Analyze results 
library(rstan)
library(ggplot2)
library(reshape)
library(RSQLite)
source("plotting_functions.R")
dbResultsFilename <- "results_hybrid_test.sqlite"
db <- dbConnect(SQLite(), dbResultsFilename)
tables <- dbListTables(db)
cor_time <- dbReadTable(db, tables[[1]])
cor_patient <- dbReadTable(db, tables[[2]])
summary <- dbReadTable(db,tables[[5]])
dbDisconnect(db) 

### read in fixed params ###########
m <- 2                   # species
n_timesteps <- c(3,5,10,20,30)          # visits/repeat observations
n_site <- c(20,40,60,80,100,150,200,300,400,500) 
version_corr <- c(1,2,3)

sweep_df <- data.frame( n_time = rep(n_timesteps, each = length(n_site)), n_s = rep(n_site, each = length(version_corr)), version = rep(version_corr, times = length(n_timesteps)*length(n_site)))
sweep_df <- sweep_df[sweep_df$version == 1,]
sweep_df$trial <- 1:nrow(sweep_df)

Rp_filename <- paste0("Rp_",1,".csv" )
Ro_filename <- paste0("Ro_",1,".csv" )

Rp = as.matrix(read.csv(Rp_filename))[,-1]
Ro = as.matrix(read.csv(Ro_filename))[,-1]

sweep_df_1 <- sweep_df[sweep_df$version == 1,]
cor_patient_1 <- cor_p_total[cor_p_total$trial %in% sweep_df_1$trial,]
cor_time_1 <- cor_t_total[cor_t_total$trial %in% sweep_df_1$trial,]

sweep_df_2 <- sweep_df[sweep_df$version == 2,]
cor_patient_2 <- cor_p_total[cor_p_total$trial %in% sweep_df_2$trial,]
cor_time_2 <- cor_t_total[cor_t_total$trial %in% sweep_df_2$trial,]

sweep_df_3 <- sweep_df[sweep_df$version == 3,]
cor_patient_3 <- cor_p_total[cor_p_total$trial %in% sweep_df_3$trial,]
cor_time_3 <- cor_t_total[cor_t_total$trial %in% sweep_df_3$trial,]


###### Summary ##############################
names <- unique(summary$row_names)[1:8]
summary <- summary[summary$row_names %in% names,]
trials <- unique(summary$trial)
sum <- summary[!is.na(summary$Rhat),]
trials_NC <- unique(sum[sum$Rhat>1.1,]$trial)
trials_C <- trials[!(trials %in% trials_NC)]



#######
cor_patient <- cor_patient_1
cor_time <- cor_time_1


## Make traceplots 
trials <- unique(cor_patient$trial)
for( i in 1:length(trials)){
  cat("i is ", i)
  df_sub <- cor_patient[cor_patient$trial == trials[i],]
  t <- make_traceplot(df_sub, "X2.1" ,2000,4)
  t <- t + ggtitle(paste0("Rp n_pat = ", sweep_df[sweep_df$trial == trials[i],]$n_s, ", n_vis = ", sweep_df[sweep_df$trial == trials[i],]$n_time))
  assign(paste0("t_",trials[i]), t)
}


trials <- trials_C
res_df <- data.frame()
for( i in 1:length(trials)){
  #cor_p <- cor_patient[cor_patient$trial == trials[i],]
  #cor_p$trial <- NULL
  cor_t <- cor_time[cor_time$trial == trials[i],]
  cor_t$trial <- NULL
  means <- apply(cor_t,2,mean)
  LCI <- apply(cor_t,2,quantile,c(.025,.975))[1,]
  UCI <- apply(cor_t,2,quantile,c(.025,.975))[2,]
  df <- data.frame(var = names(means),
                   mean = means,
                   LCI = LCI,
                   UCI = UCI,
                   trial = sweep_df[sweep_df$trial == trials[i],]$trial,
                   n_times = sweep_df[sweep_df$trial == trials[i],]$n_time,
                   n_pat = sweep_df[sweep_df$trial == trials[i],]$n_s
                   )
res_df <- rbind(res_df, df)
}



# plotting 
cor_dim <- c(2,1)
var_name <- paste0("X",cor_dim[1],".",cor_dim[2])
test <- res_df[res_df$var == var_name,]
test$true_val <- Ro[cor_dim[1],cor_dim[2]]
p <- ggplot(test, aes(x = n_pat, y = mean)) + geom_point() + geom_hline(aes(yintercept =true_val, color = "red")) + geom_errorbar(aes(ymin = LCI, ymax = UCI, width = .1)) + facet_wrap(~n_times, scale = "free")
p + theme_bw() +  ggtitle(paste0("estimated vs. true correlation r_time for ",var_name)) + ggsave(paste0(var_name, "Ro_version_3.png"))

