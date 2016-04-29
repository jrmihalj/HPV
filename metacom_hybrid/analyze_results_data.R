# Analyze results 
library(rstan)
library(ggplot2)
library(reshape)
library(RSQLite)

dbResultsFilename <- "results_hybrid_4_29_200_pat_LKJ10.sqlite"
db <- dbConnect(SQLite(), dbResultsFilename)
tables <- dbListTables(db)
cor_time <- dbReadTable(db, tables[[1]])
cor_patient <- dbReadTable(db, tables[[2]])
beta_gam <- dbReadTable(db, tables[[3]])
beta_phi <- dbReadTable(db, tables[[4]])
summary <- dbReadTable(db,tables[[5]])
dbDisconnect(db) 

make_traceplot <- function(results, param ,iter,nChains){
  test <- results[,names(results) == param]
  cut_point <- iter - 100
  df <- data.frame(
    matrix(test, ncol = nChains)
  )
  chain_names <- paste0("chain_",c(1:nChains))
  names(df) <- chain_names
  df$iter = c(1:nrow(df))
  dfm <- melt(df, id.vars = "iter")
  title <- paste("traceplot ", param, sep = "")
  traceplot <- ggplot(dfm, aes(x=iter, y = value, color = variable)) + geom_line()+ ggtitle(title) + theme_bw()
  return(traceplot)
}  



## Summary #####-------------------------------------------------------------------------------------------

sum <- summary[!is.na(summary$Rhat),]

#### -------------------------------------------------------------------------------------------------------

## Make traceplots 
t <- make_traceplot(cor_patient, "X2_1" ,2000,4)
t + ggtitle("cor_patient_sp1_sp2")


trials <- trials_C
res_df <- data.frame()
for( i in 1:length(trials)){
  cor_p <- cor_patient[cor_patient$trial == trials[i],]
  cor_p$trial <- NULL
  #cor_t <- cor_time[cor_time$trial == trials[i],]
  #cor_t$trial <- NULL
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
beta_phi_intxn_true <- as.matrix(beta_phi_intxn_list[[beta_phi_intxn_ind]])

cor_dim <- c(2,1)
var_name <- paste0("X",cor_dim[1],".",cor_dim[2])
test <- res_df[res_df$var == var_name,]
test$true_val <- Rp[cor_dim[1],cor_dim[2]]
p <- ggplot(test, aes(x = n_pat, y = mean)) + geom_point() + geom_hline(aes(yintercept =true_val, color = "red")) + geom_errorbar(aes(ymin = LCI, ymax = UCI, width = .1)) + facet_wrap(~n_times, scale = "free")
p + theme_bw() +  ggtitle(paste0("estimated vs. true correlation r_patient for ",var_name)) + ggsave(paste0(var_name, "Rp_version_1.png"))



