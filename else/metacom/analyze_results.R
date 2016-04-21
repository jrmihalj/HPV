# Analyze results 
library(rstan)
library(ggplot2)
library(reshape)
library(RSQLite)
dbResultsFilename <- "results.sqlite"
db <- dbConnect(SQLite(), dbResultsFilename)
tables <- dbListTables(db)
cor_patient <- dbReadTable(db, tables[[1]])
cor_time <- dbReadTable(db, tables[[2]])
dbDisconnect(db)


### read in fixed params ###########
m <- 2                   # species
n_timesteps <- c(3,5,10,20,30)          # visits/repeat observations
n_site <- c(20,40,60,80,100,150,200,300,400,500)  
sweep_df <- data.frame( n_time = rep(n_timesteps, each = length(n_site)), n_s = rep(n_site, times = length(n_timesteps)))
sweep_df$trial <- 1:nrow(sweep_df)
Rp = as.matrix(read.csv("Rp.csv"))[,-1]
Ro = as.matrix(read.csv("Ro.csv"))[,-1]

#######
trials <- unique(cor_patient$trial)
res_df <- data.frame()
for( i in 1:length(trials)){
  cor_p <- cor_patient[cor_patient$trial == trials[i],]
  cor_p$trial <- NULL
  means <- apply(cor_p,2,mean)
  LCI <- apply(cor_p,2,quantile,c(.025,.975))[1,]
  UCI <- apply(cor_p,2,quantile,c(.025,.975))[2,]
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
test$true_val <- Rp[cor_dim[1],cor_dim[2]]
p <- ggplot(test, aes(x = n_pat, y = mean)) + geom_point() + geom_hline(aes(yintercept =true_val, color = "red")) + geom_errorbar(aes(ymin = LCI, ymax = UCI, width = .1)) + facet_wrap(~n_times, scale = "free")
p + theme_bw() +  ggtitle(paste0("estimated vs. true correlation r_patient for ",var_name)) + ggsave(paste0(var_name, ".png"))

