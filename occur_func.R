
occur_func <- function(data){
  
  # Storage
  co <- mat.or.vec(nr = nrow(data)/2, nc = 2) #dimensions = Visit.Pats, alt/null, strain 1/strain 2
  per <- array(0, dim=c(nrow(data)/2, 2, 2))
  col <- array(0, dim=c(nrow(data)/2, 2, 2))
  
  coinf <- NULL # just a vector of size 2
  persist <- mat.or.vec(nr = 2, nc = 2)
  colonize <- mat.or.vec(nr = 2, nc = 2)
  
  # Total co-infections
  j <- 1
  k <- 1
  while(j < nrow(data)){
    if (data$Y_Alt[j] == 1 && data$Y_Alt[j+1] == 1 ){
      co[k,1] <- 1
    } else{co[k,1] <- 0}
    
    if (data$Y_Null[j] == 1 && data$Y_Null[j+1] == 1 ){
      co[k,2] <- 1
    } else{co[k,2] <- 0}
    
    k <- k + 1
    j <- j + 2
  }
  coinf[1] <- sum(co[,1]) # total coinf in Alt
  coinf[2] <- sum(co[,2]) # total coinf in Null
  
  #Persistence and Colonization of each strain
  j <- 1
  for(p in 1:max(data$Patient)){ # For each patient
    sub <- data[data$Patient == p, ]
    for(t in 1:(max(sub$Visit.Pat)-1)){ # For each Visit.Pat
      
      ###### FOR STRAIN 1 ##########
      # If at Visit.Pat 't', s1 & s2 == 1, & at Visit.Pat 't+1' s1 == 1, then persistence = 1
      if( (sub[sub$Visit.Pat == t & sub$Strain == 1, ]$Y_Alt == 1 &
             sub[sub$Visit.Pat == t & sub$Strain == 2, ]$Y_Alt == 1)
          & (sub[sub$Visit.Pat == t+1 & sub$Strain == 1, ]$Y_Alt == 1)){
        per[j, 1, 1] <- 1 }else{ per[j, 1, 1] <- 0 }
      if( (sub[sub$Visit.Pat == t & sub$Strain == 1, ]$Y_Null == 1 &
             sub[sub$Visit.Pat == t & sub$Strain == 2, ]$Y_Null == 1)
          & (sub[sub$Visit.Pat == t+1 & sub$Strain == 1, ]$Y_Null == 1)){
        per[j, 2, 1] <- 1 }else{ per[j, 2, 1] <- 0 }
      # If at Visit.Pat 't' s1 == 0 & s2 == 1, & at Visit.Pat 't+1' s1 == 1, then colonize = 1
      if( (sub[sub$Visit.Pat == t & sub$Strain == 1, ]$Y_Alt == 0 &
             sub[sub$Visit.Pat == t & sub$Strain == 2, ]$Y_Alt == 1)
          & (sub[sub$Visit.Pat == t+1 & sub$Strain == 1, ]$Y_Alt == 1)){
        col[j, 1, 1] <- 1}else{ col[j, 1, 1] <- 0}
      if( (sub[sub$Visit.Pat == t & sub$Strain == 1, ]$Y_Null == 0 &
             sub[sub$Visit.Pat == t & sub$Strain == 2, ]$Y_Null == 1)
          & (sub[sub$Visit.Pat == t+1 & sub$Strain == 1, ]$Y_Null == 1)){
        col[j, 2, 1] <- 1}else{ col[j, 2, 1] <- 0}
      
      ###### FOR STRAIN 2 ##########
      # If at Visit.Pat 't', s1 & s2 == 1, & at Visit.Pat 't+1' s2 == 1, then persistence = 1
      if( (sub[sub$Visit.Pat == t & sub$Strain == 1, ]$Y_Alt == 1 &
             sub[sub$Visit.Pat == t & sub$Strain == 2, ]$Y_Alt == 1)
          & (sub[sub$Visit.Pat == t+1 & sub$Strain == 2, ]$Y_Alt == 1)){
        per[j, 1, 2] <- 1 }else{ per[j, 1, 2] <- 0 }
      if( (sub[sub$Visit.Pat == t & sub$Strain == 1, ]$Y_Null == 1 &
             sub[sub$Visit.Pat == t & sub$Strain == 2, ]$Y_Null == 1)
          & (sub[sub$Visit.Pat == t+1 & sub$Strain == 2, ]$Y_Null == 1)){
        per[j, 2, 2] <- 1 }else{ per[j, 2, 2] <- 0 }
      # If at Visit.Pat 't' s2 == 0 & s1 == 1, & at Visit.Pat 't+1' s2 == 1, then colonize = 1
      if( (sub[sub$Visit.Pat == t & sub$Strain == 1, ]$Y_Alt == 1 &
             sub[sub$Visit.Pat == t & sub$Strain == 2, ]$Y_Alt == 0)
          & (sub[sub$Visit.Pat == t+1 & sub$Strain == 2, ]$Y_Alt == 1)){
        col[j, 1, 2] <- 1}else{ col[j, 1, 2] <- 0}
      if( (sub[sub$Visit.Pat == t & sub$Strain == 1, ]$Y_Null == 1 &
             sub[sub$Visit.Pat == t & sub$Strain == 2, ]$Y_Null == 0)
          & (sub[sub$Visit.Pat == t+1 & sub$Strain == 2, ]$Y_Null == 1)){
        col[j, 2, 2] <- 1}else{ col[j, 2, 2] <- 0}
      
      
      j <- j + 1
      
    }
  }
  
  for(i in 1:2){
    for(z in 1:2){
      persist[i,z] <- sum(per[,i,z], na.rm=T)
      colonize[i,z] <- sum(col[,i,z], na.rm=T)
      
    }
  }
  
  
  L <- list(coinf = coinf, 
                persist = persist, 
                colonize = colonize)
  
 return(L)

}

