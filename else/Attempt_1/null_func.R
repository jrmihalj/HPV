# Null occurrence indices

null_func <- function(data, n.iter){
  
  require(tidyr)
  
  S1 <- data[data$Strain==1,]
  S2 <- data[data$Strain==2,]
  
  freq1 <- sum(S1$Y)/length(S1$Y)
  freq2 <- sum(S2$Y)/length(S2$Y)
  
  # Calculate the occurrence metrics:
  occur <- occur_func(data)
  
  co <- NULL
  
  for(i in 1:n.iter){
    
    null1 <- rbinom(nrow(S1), 1, freq1)
    null2 <- rbinom(nrow(S2), 1, freq2)
    
    mat <- data.frame(null1, null2, S1$Patient, S1$Visit)
    colnames(mat) <- c("1", "2", "Patient", "Visit.Pat")
    
    # Restructure to match the output of sim_func()
    # This is required to work with occur_func
    data.null <- gather(mat, key=Strain, value=Y, 1:2)
    
    null.oc <- occur_func(data.null)
    co[i] <- null.oc$coinf
  
  }
  
  co.mean.null <- mean(co)
  co.sd.null <- sd(co)
  
  # Estimate how different the real data is from the null:
  # The null ~ Normal (I checked...)
  co.dif <- pnorm(occur$coinf, mean = co.mean.null, sd = co.sd.null)
  
  L <- list(coinf = occur$coinf, 
            co.mean.null = co.mean.null, 
            coinf.sd.null = co.sd.null,
            co.dif = co.dif)
  
  return(L)
  
}


# NOT USED: WITH PERSIST AND COLON...

# 
# null_func <- function(data, n.iter){
#   
#   require(tidyr)
#   
#   S1 <- data[data$Strain==1,]
#   S2 <- data[data$Strain==2,]
#   
#   freq1 <- sum(S1$Y)/length(S1$Y)
#   freq2 <- sum(S2$Y)/length(S2$Y)
#   
#   # Calculate the occurrence metrics:
#   occur <- occur_func(data)
#   
#   co <- NULL
#   per <- NULL
#   col <- NULL
#   
#   for(i in 1:n.iter){
#     
#     null1 <- rbinom(nrow(S1), 1, freq1)
#     null2 <- rbinom(nrow(S2), 1, freq2)
#     
#     mat <- data.frame(null1, null2, S1$Patient, S1$Visit)
#     colnames(mat) <- c("1", "2", "Patient", "Visit.Pat")
#     
#     # Restructure to match the output of sim_func()
#     # This is required to work with occur_func
#     data.null <- gather(mat, key=Strain, value=Y, 1:2)
#     
#     null.oc <- occur_func(data.null)
#     co[i] <- null.oc$coinf
#     per[i] <- null.oc$persist
#     col[i] <- null.oc$colonize
#     
#   }
#   
#   co.mean.null <- mean(co)
#   co.sd.null <- sd(co)
#   per.mean.null <- mean(per)
#   per.sd.null <- sd(per)
#   col.mean.null <- mean(col)
#   col.sd.null <- sd(col)
#   
#   # Estimate how different the real data is from the null:
#   # The nulls all follow normals (I checked...)
#   co.dif <- pnorm(occur$coinf, mean = co.mean.null, sd = co.sd.null)
#   per.dif <- pnorm(occur$persist, mean = per.mean.null, sd = per.sd.null)
#   col.dif <- pnorm(occur$colonize, mean = col.mean.null, sd = col.sd.null)
#   
#   
#   L <- list(coinf = occur$coinf, 
#             persist = occur$persist, 
#             colonize = occur$colonize,
#             
#             co.mean.null = co.mean.null, 
#             coinf.sd.null = co.sd.null,
#             persist.mean.null = per.mean.null, 
#             persist.sd.null = per.sd.null,
#             colonize.mean.null = col.mean.null,
#             colonize.sd.null = col.sd.null,
#             
#             co.dif = co.dif,
#             per.dif = per.dif,
#             col.dif = col.dif)
#   
#   return(L)
#   
# }