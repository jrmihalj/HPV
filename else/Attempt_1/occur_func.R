
##############################
# Current measure of persist:
##############################

#     Strain 1   Strain 2
# t       1         1
# t+1     1

##############################
# Current measure of colonize:
##############################

#     Strain 1   Strain 2
# t       0         1
# t+1     1

occur_func <- function(data){
  
  # Storage
  co <- NULL
  
  coinf <- NULL 
  
  S1 <- data[data$Strain==1,]
  S2 <- data[data$Strain==2,]
  Mat <- cbind(S1$Y, S2$Y)
  
  ######################
  # Total co-infections:
  ######################
  
  co <- ifelse(Mat[,1] == 1 & Mat[, 2] == 1, 1, 0)
  coinf <- sum(co)
  
  L <- list(coinf = coinf)
 
  return(L)

}

# # NOT USED RIGHT NOW:
#  # Storage
# co <- NULL
# per <- NULL
# col <- NULL
# 
# coinf <- NULL 
# persist <- NULL
# colonize <- NULL
# ############################################
# #Persistence and Colonization of each strain
# ############################################
# 
# j <- 1
# 
# for(p in 1:max(data$Patient)){ # For each patient
#   
#   sub <- data[data$Patient == p, ]
#   mat <- cbind(sub[sub$Strain == 1, ]$Y, sub[sub$Strain == 2, ]$Y)
#   
#   t <- 1
#   while(t < nrow(mat)){
#     
#     ###### JUST FOR STRAIN 1 ##########
#     
#     # Persistence
#     # If at Visit 't', s1 & s2 == 1, & at Visit 't+1' s1 == 1, then persistence = 1
#     if( (mat[t, 1] == 1 & mat[t, 2] == 1) & mat[t+1, 1] == 1){
#       per[j] <- 1 }else{ per[j] <- 0 }
#     
#     # Colonization
#     # If at Visit 't' s1 == 0 & s2 == 1, & at Visit 't+1' s1 == 1, then colonize = 1
#     if( (mat[t, 1] == 0 & mat[t, 2] == 1) & mat[t+1, 1] == 1){
#       col[j] <- 1}else{ col[j] <- 0}     
#     
#     j <- j + 1
#     t <- t + 1
#   }
#   
# }
# 
# persist <- sum(per)
# colonize <- sum(col)
# 
# L <- list(coinf = coinf, 
# persist = persist, 
# colonize = colonize)
# 
