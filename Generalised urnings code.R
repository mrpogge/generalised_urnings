#### Code for the illustrative example

# load relevant packages
library(dplyr)
library(tidyr)
library(purrr)
library(Rcpp)
library(tictoc)
library(parallel)

# Load custom cpp functions
sourceCpp('urnings_empirical.cpp')

# R functions

# Function which runs the urnings algorithm through the data for a specific 'model' ('Rasch','DSRT','PCM','NRM','MNRM')
# and return either the MSE of prediction of the item score or the MSE of prediction of the binary accuracy (binary=FALSE or TRUE)
# the input is a vector of urn sizes, data (should have $Item, $Person and $X - 4-category item score), 
# argument 'scores' (s_1 and s_2 scores for the slow and fast correct responses) is only relevant for NRM

urnings_math_garden=function(n,data,scores=NULL,model,binary,MSE=TRUE,person_specific=FALSE,expected=FALSE){
  P=max(data$Pers) # number of students
  I=max(data$Item) # number of items
  J=length(n)
  if(model!='MNRM'){
    D=1
  }
  if(model!='DSRT'){
    # the matrices b are put next to each other other
    b=matrix(c(0,1,1,0),2,2)
    colnames(b)<-c(-1,1) # for each column of b indicate which b_c it belongs to
    rownames(b)<-c(-1,1)
  }
  if(model=="Rasch"){
    s=cbind(c(0,1),c(1,0))
    # indication that 1st item parameter is included in set w={1}
    paired=1
    Item_score=1*(data$X>1)
  }
  if(model=='DSRT'){
    s=cbind(c(0,1,2,3),c(3,2,1,0))
    b=cbind(c(0,0,0,0,0,1),c(0,0,0,1,1,0),c(0,0,0,3,0,0),c(0,0,0,0,1,0),c(0,0,0,2,0,0),c(0,0,0,1,0,0),c(0,0,1,0,0,0),c(0,1,0,0,0,0),c(0,0,2,0,0,0),c(1,0,0,0,0,0),c(0,1,1,0,0,0),c(0,0,3,0,0,0))
    colnames(b)=c(-3,-3,-3,-2,-2,-1,1,2,2,3,3,3)
    rownames(b)=c(-3,-2,-1,1,2,3)
    paired=1
    Item_score=data$X
  }
  if(model=='PCM'){
    s=cbind(c(0,1,2,3),c(1,0,0,0),c(1,1,0,0),c(1,1,1,0))
    #indication that the 1st and the 3rd item parameters are not included in set w, but the 2nd is (w={2})
    paired=c(0,1,0)
    Item_score<-data$X
  }
  if(model=='NRM'){
    s=cbind(c(0,scores[1],scores[2]),c(1,0,0),c(1,1,0))
    paired=c(0,1)
    Item_score<-data$X-1
    Item_score[Item_score<0]=0
  }
  if(model=='MNRM'){
    D=2
    s=cbind(c(0,1,2,3),c(1,0,0,1),c(1,0,0,0),c(1,1,0,0),c(1,1,1,0))
    paired=c(1,1,0)# i.e., w=c(1,2)
    Item_score<-data$X
  }
  M=max(Item_score)+1
  S=apply(s,2,max)
  tilde_t=as.numeric(rownames(b)) # vector of possible non-zero updates
  K=length(tilde_t) # number of possible non-zero updates
  
  # specify where in the combined matrix b the columns for each element of tilde_t start 
  cn=colnames(b)
  ucn=rownames(b)
  blocks=NULL
  for(j in 1:length(ucn)){
    blocks=c(blocks,which(cn==ucn[j])[1])
  }
  blocks=c(blocks-1,ncol(b))
  
  # Specifying starting values
  # person urnings; we start ability at 0.7
  if(model=='MNRM'){
    u=c(round(0.7*n[1]),round(0.5*n[2]))
  }
  if(model!='MNRM'){
    u=round(0.7*n[1:D])
  }
  # item urnings
  v=0.5*n[(D+1):J]
  
  res=urnings_empirical(X=as.integer(Item_score),
                        Pers=as.integer(data$Pers-1),
                        Item=as.integer(data$Item-1),
                        n_students=as.integer(P),
                        n_items=as.integer(I),
                        n_pars=as.integer(J), 
                        n_dims=as.integer(D), 
                        n_resp=as.integer(nrow(data)), 
                        n_categories=as.integer(M), 
                        n_scores=as.integer(K),
                        start_u=as.integer(u),
                        start_v=as.integer(v),
                        S=as.integer(s), 
                        Smax=as.integer(S),
                        urn_size=as.integer(n), 
                        blocks=as.integer(blocks), 
                        Upd=as.integer(b),
                        Score=as.integer(tilde_t), 
                        paired=as.integer(paired)) 
  if(MSE==TRUE){
    if(binary==TRUE){
      if(person_specific==TRUE){
        output=mean(res[[4]])
      }else{
      output=res[[6]]}
    }else{
      if(person_specific==TRUE){
        output=mean(res[[3]])
      }else{
      output=res[[5]]}
    }
  }else{ # here the urning values are returned
    output=matrix(res[[2]],ncol=J,byrow=TRUE)
  }
  if(expected==TRUE){output=res[[7]]}
  return(output)
}  
  
# Function which runs the urnings function multiple times across multiple cores and returns the mse
MSE_rep <- function(x,data,scores=NULL,model,binary, reps = 200, cores = 20) {
  
  # Run one of the MSE functions in parallel reps times
  results <- mclapply(
    X = 1:reps,
    FUN = function(i) urnings_math_garden(x,data,scores,model,binary),
    mc.cores = cores
  )
  
  # Convert list to numeric vector and average
  mean(unlist(results))
}

# Function for optimizing urn size via coordinate pattern search
optimize_urn_size <- function(data,
                              scores=NULL,
                              model,
                              binary,
                              start = 128, 
                              start_step = 64, 
                              max_iter = 50) {
  
  # specify the number of parameters (and urn sizes)
  if((model=='Rasch')|(model=='DSRT')){J=2}
  if((model=='NRM')){J=3}
  if(model=="PCM"){J=4}
  if(model=='MNRM'){J=5}
  
  # First start form non-informed starting values
  x=rep(start,J)
  
  # run coordinate-wise search where we run the urnings and compute MSE once
  for(k in 0:5){
    current_val <- urnings_math_garden(x,data=data,scores=scores,binary=binary,model=model)
    step=start_step/(2^k)
    
    for (iter in 1:max_iter) {
      
      moves <- list()
      for(i in 1:(2*J)){
        moves[[i]]=x
      }
      for(i in 1:J){
        moves[[i]][i]=x[i]+step
        moves[[i+J]][i]=x[i]-step
      }
      moves <- Filter(function(x) all(x >= 4), moves)
      
      # Evaluate neighbors
      vals <- sapply(moves, function(v)urnings_math_garden(v,scores=scores,data=data,model=model,binary=binary))
      
      # Best neighbor
      best_idx <- which.min(vals)
      best_move <- moves[[best_idx]]
      best_val <- vals[best_idx]
      
      # If no improvement, stop
      if (best_val >= current_val) {
        break
      }
      
      # Move to the improving neighbor
      x=best_move  
      current_val=best_val
      
      print(x)
      print(best_val)
    }
    
  }
  
  # From the resulting point, run coordinate-wise search with 200 replications
  
  for(k in 0:5){
    current_val <- MSE_rep(x,data=data,scores=scores,binary=binary,model=model)
    step=start_step/(2^k)
    
    for (iter in 1:max_iter) {
      
      moves <- list()
      for(i in 1:(2*J)){
        moves[[i]]=x
      }
      for(i in 1:J){
        moves[[i]][i]=x[i]+step
        moves[[i+J]][i]=x[i]-step
      }
      moves <- Filter(function(x) all(x >= 4), moves)
      
      # Evaluate neighbors
      vals <- sapply(moves, function(v) MSE_rep(v,scores=scores,data=data,model=model,binary=binary))
      
      # Best neighbor
      best_idx <- which.min(vals)
      best_move <- moves[[best_idx]]
      best_val <- vals[best_idx]
      
      # If no improvement, stop
      if (best_val >= current_val) {
        break
      }
      
      # Move to the improving neighbor
      x=best_move  
      current_val=best_val
      
      print(x)
      print(best_val)
    }
    
  }
  
  return(list(x=x, fval = current_val))
}

#================
# Load data
#================
load(file="~/research-collaboration/data-ro/data-old-new-scoring-rule/tables-2024-10-03.RData")

#===================================
# Prepare the data for the analysis
#===================================

# select data within one academic year, where the students chose to see the coins, and the time limit was 10 seconds
# select relevant variables (user ID, item ID, timestamp, response accuracy and response time)
logs <- logs %>%
  filter(created >= '2018-09-01', 
         created < '2019-07-01', 
         show_coins == 1,
         maximum_response_in_seconds.x==10) %>%
  select(user_id, created, response_in_milliseconds, correct_answered,item_id)

# remove persons with less than 20 responses
logs <- logs %>%
  group_by(user_id) %>% 
  filter(n() >= 20 )%>% 
  ungroup() %>%
  as.data.frame()

# count the number of persons
N=length(unique(logs$user_id))

# Divide the sample in two folds
# Specify a vector with fold values of length N
if(N%%2==0){
  folds=rep(c(1:2),N/2)
}
if(N%%2!=0){
  folds=c(rep(c(1:2),(N-N%%2)/2),c(1:(N%%2)))
}

# create the 4-category item score variable
logs<-logs%>%
  mutate(X=(response_in_milliseconds>5000)*(correct_answered==0)*1+
           (response_in_milliseconds>5000)*(correct_answered==1)*2+
           (response_in_milliseconds<=5000)*(correct_answered==1)*3)

# convert item IDs to numbers from 1 to 100
logs$Item <- match(logs$item_id, unique(logs$item_id))

# make sure that the item score is an integer
logs$X<-as.integer(logs$X)

logs<-logs%>%select(Item,X,user_id)

# randomly divide the persons into the folds
set.seed(108501865)
folds=sample(folds,N)
# write the userIDs for different folds in a list
folds_list=list()
for(i in 1:2){folds_list[[i]]=unique(logs$user_id)[which(folds==i)]}

# For saving the results across folds for the prediction of the 4-category item
MSE=matrix(ncol=length(folds_list),nrow=4)
Urn_sizes=list()

# for the prediction of the binary scores
MSE2=matrix(ncol=length(folds_list),nrow=9)
Urn_sizes2=list()

#Part 1: Models for 4-category items

Models=c('DSRT','PCM','MNRM')

load(file='polytomous_analysis.RData')
for(II in 1:2){
  
  data<-logs%>%
    filter(!user_id%in%folds_list[[II]])
  
  data$Pers <- match(data$user_id, unique(data$user_id))
  
  optimum=list()
  
  for(J in 1:length(Models)){
    optimum[[J]]<-optimize_urn_size(data=data,model=Models[J],binary=FALSE)[[1]]
  }
  
  data<-logs%>%
    filter(user_id%in%folds_list[[3-II]])
  
  data$Pers <- match(data$user_id, unique(data$user_id))
  
  for(J in 1:3){
    MSE[J,II]=MSE_rep(optimum[[J]],data=data,model=Models[J],binary=FALSE)
  }
  
  # baseline prediction
  P=table(data$X)
  P/sum(P)
  sum((P/sum(P))*c(0:3))->overall
  MSE[4,II]=mean((data$X-overall)^2)
  
  Urn_sizes[[II]]=optimum
  
  save(MSE,Urn_sizes,file='polytomous_analysis.RData')
  
}

#Part 2: Comparing models in terms of prediction of binary scores

Models=c("Rasch",'DSRT','PCM',"MNRM","NRM","NRM","NRM",'NRM')
Scores=list(NULL,NULL,NULL,NULL,c(1,2),c(1,3),c(2,3),c(3,4))

for(II in 1:2){
  
  data<-logs%>%
    filter(!user_id%in%folds_list[[II]])
  
  data$Pers <- match(data$user_id, unique(data$user_id))
  
  optimum=list()
  
  for(J in 1:length(Models)){
    optimum[[J]]<-optimize_urn_size(data=data,model=Models[J],binary=TRUE,scores=Scores[[J]])[[1]]
  }
  
  data<-logs%>%
    filter(user_id%in%folds_list[[3-II]])
  
  data$Pers <- match(data$user_id, unique(data$user_id))
  
  for(J in 1:8){
    MSE2[J,II]=MSE_rep(optimum[[J]],data=data,scores=Scores[[J]],model=Models[J],binary=TRUE)
  }
  
  # baseline prediction
  P=table(1*(data$X>1))
  sum((P/sum(P))*c(0:1))->overall
  MSE2[9,II]=mean((1*(data$X>1)-overall)^2)
  
  Urn_sizes2[[II]]=optimum
  
  save(MSE2,Urn_sizes2,file='binary_analysis.RData')
  
}


# Collecting relevant output from the models
P_fast<-P_slow<-dif_fast<-dif_slow<-0
Cors=matrix(0,ncol=3,nrow=3)

for(II in 1:2){
  
  data<-logs%>%
    filter(!user_id%in%folds[[II]])
  
  data$Pers <- match(data$user_id, unique(data$user_id))
  
  R=urnings_math_garden(Urn_sizes[[3-II]][[3]],data=data,binary=NULL,model='MNRM',MSE=FALSE)
  
  pp=NULL # positions in the dataset where each person's last response is
  for(i in 1:max(data$Pers)){
    pp[i]=rev(which(data$Pers==i))[1]
  }
  
  ii=NULL # positions in the dataset where each item's last response is
  for(j in 1:100){
    ii[j]=rev(which(data$Item==j))[1]
  }
  
  # final person urnings
  R_pers=(R[pp,1:2])
  # final item urnings
  R_item=(R[ii,3:5])
  
  # based on the final ratings define fast and slow students
  which(R_pers[,2]>=quantile(R_pers[,2],2/3))->fast
  which(R_pers[,2]<quantile(R_pers[,2],1/3))->slow
  
  # for each fast person compute their proportion of fast responses
  p_fast=NULL;for(i in fast){p_fast=c(p_fast,mean(data$X[data$Pers==i]==3)+mean(data$X[data$Pers==i]==0))}
  # for each slow person compute their proportion of fast responses
  p_slow=NULL;for(i in slow){p_slow=c(p_slow,mean(data$X[data$Pers==i]==3)+mean(data$X[data$Pers==i]==0))}
  
  # Average proportions of fast response for the fast and slow persons (Averaged across the folds)
  P_slow=P_slow+1/2*mean(p_slow)
  P_fast=P_fast+1/2*mean(p_fast)
  
  X <- as.matrix(R_item)
  n <- as.numeric(Urn_sizes[[3-II]][[3]][3:5])
  
  N <- nrow(X)
  K <- ncol(X)
  
  # Compute the correlations between the item parameters corrected for attenuation
  
  # Step 1: observed proportions
  Theta_hat <- sweep(X, 2, n, "/")   # divide each column by its n_k
  
  # Step 2: observed covariance of proportions
  Cov_hat <- cov(Theta_hat)
  
  # Step 3: estimate measurement-noise variances
  noise_var <- numeric(K)
  for (k in 1:K) {
    noise_var[k] <- mean(Theta_hat[, k] * (1 - Theta_hat[, k]) / n[k])
  }
  
  # Step 4: de-noised latent covariance matrix
  Latent_Cov <- Cov_hat
  diag(Latent_Cov) <- pmax(diag(Cov_hat) - noise_var, 0)
  
  # Step 5: convert to correlation matrix
  Latent_Corr <- cov2cor(Latent_Cov)
  print(Latent_Corr)
  
  # Correlation matrix between the item parameters, averaged across the folds
  Cors=Cors+1/2*Latent_Corr
  
  R_PCM=urnings_math_garden(Urn_sizes[[3-II]][[2]],data=data,binary=NULL,model='PCM',MSE=FALSE)
  
  # ability estimates under PCM
  R_pers_PCM=R_PCM[pp,1]
  
  # convert the final urnings under the MNRM and PCM into proportions and compute the differences between
  # the models' ability estimates for the slow and fast persons
  # average the results across the folds
  dif_slow=dif_slow+1/2*mean(R_pers[slow,1]/Urn_sizes[[3-II]][[3]][1]-R_pers_PCM[slow]/Urn_sizes[[3-II]][[2]][1])
  dif_fast=dif_fast+1/2*mean(R_pers[fast,1]/Urn_sizes[[3-II]][[3]][1]-R_pers_PCM[fast]/Urn_sizes[[3-II]][[2]][1])
  
  save(P_slow,P_fast,dif_slow,dif_fast,Cors,file='Model_results.RData')
  
}


