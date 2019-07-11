# load packages
require(poolSeq)


args = commandArgs(trailingOnly=TRUE)

L <- as.numeric(args[1])
dT <- as.numeric(args[2])
N <- as.numeric(args[3])
ns <- as.numeric(args[4])
T <- as.numeric(args[5])
nosim <- as.numeric(args[6])
test_name <- args[7]
dir_name <- args[8]
mu <- as.numeric(args[9])

flag<-1
ite_sim <-1

s_LLS_mat <- array(0, c(nosim, L));
array_LLS_time<- array(0,c(nosim,1)) # time to run LLS

s_MPL_mat <- array(0, c(nosim, L));
array_MPL_time<- array(0,c(nosim,1)) # time to run MPL

s_SL_mat <- array(0, c(nosim, L));
array_SL_time<- array(0,c(nosim,1)) # time to run SL
array_load_time<- array(0,c(nosim,1)) # time to load data
array_load_time_onetwo<- array(0,c(nosim,1)) # time to load data
array_load_time_one<- array(0,c(nosim,1)) # time to load data

while (ite_sim <= nosim){


  start_time <- Sys.time()

  dat_name <- paste("wfsim_",test_name,"_",ite_sim-1,"_T",T,"_ns",ns,"_dt",dT,".dat",sep="")
  pathname <- paste(dir_name,dat_name,sep="")


  input_data <- read.table(pathname, sep="\t")
  colnames(input_data) <- c("generation", "count","sequences")

  end_time <- Sys.time()


  array_load_time[ite_sim] <-end_time-start_time

  #################################3
  # obtain one/two points

  start_time <- Sys.time()

  data_gen <- input_data["generation"]
  unique_gen <- unlist(unique(data_gen))
  num_time_points <- length(unique_gen) # number of observed generations

  T_overwrite<- dT*num_time_points # overrides user input for T. Not including this causes errors if T is not properly specified

  q <- array(0,c(num_time_points,L))
  q11<- array(0, c(L, L,num_time_points))

  for (ite_T in 1:num_time_points){
    curr_gen <- unique_gen[ite_T]

    curr_data <- input_data[input_data$generation==curr_gen,]
    curr_seqs <- curr_data["sequences"]
    curr_count <- as.numeric(unname(unlist(curr_data["count"])))
    curr_total_count = sum(curr_count)

    num_unique_seqs = lengths(curr_seqs)

    # convert sequences to binary matrix
    curr_seqs_bin <- array(0, c(num_unique_seqs, L));
    for (ind_unique in 1:num_unique_seqs){
      curr_str<-toString(curr_seqs[ind_unique,])
      curr_str_no_space <- gsub(" ", "", curr_str, fixed = TRUE)
      for (ind_L in 1:L){
        curr_seqs_bin[ind_unique,ind_L] <- as.numeric(substring(curr_str_no_space,ind_L,ind_L))

      }
    }


    # calculate one/two points
    if (length(curr_count)==1){
      curr_mut <-t(curr_seqs_bin) %*% curr_seqs_bin
    } else {

      curr_mut <-t(curr_seqs_bin) %*% diag(curr_count) %*% curr_seqs_bin

    }
    curr_mut <- curr_mut/curr_total_count
    q[ite_T,] <- diag(curr_mut)
    q11[,,ite_T] <- curr_mut

  }

  end_time <- Sys.time()
  array_load_time_onetwo[ite_sim] <-end_time-start_time
  ####################################


  ########################3

  # LLS estimation

  start_time <- Sys.time()
  s_LLS <- array(0,c(1,L))

  for (ind_residue in 1:L){

    simTraj <- q[,ind_residue]
    s_list <- estimateSH(simTraj, Ne=N, t=seq(0, T_overwrite-1, by=dT), h=0.5,haploid=TRUE,N.ctraj=1000)
    s_LLS[ind_residue] <- s_list[["s"]]
  }
  s_LLS_mat[ite_sim,]= t(s_LLS)
  end_time <- Sys.time()
  array_LLS_time[ite_sim] <-end_time-start_time

  #############

  # MPL estimation
  start_time <- Sys.time()

  integrated_covar <- array(0, c(L, L));
  for (ite_T in 1:(num_time_points-1)){

    curr_q = q[ite_T,]
    curr_q11 = q11[,,ite_T]

    covar_mat = curr_q11 - curr_q %o% curr_q # covariance matrix
    integrated_covar = integrated_covar + covar_mat*dT

  }

  if (L==1) {
    b <- q[num_time_points]  - q[1] - mu*dT*(num_time_points-2*sum(q[1:num_time_points-1]))
    s_MPL <- b/(integrated_covar+1)
  } else{
    b <- q[num_time_points,]  - q[1,] - mu*dT*(array(num_time_points-1,c(1,L))-2*apply(q[1:num_time_points-1,],2,sum))

    s_MPL <- solve(integrated_covar+diag(L),t(b))
  }
  s_MPL_mat[ite_sim,]= t(s_MPL)

  end_time <- Sys.time()
  array_MPL_time[ite_sim] <-end_time-start_time

  #############

  # SL estimation
  start_time <- Sys.time()

  integrated_var <- array(0, c(L, 1));
  for (ite_T in 1:(num_time_points-1)){

    curr_q = q[ite_T,]

    integrated_var = integrated_var + curr_q*(array(1, c(L, 1))-curr_q)*dT

  }

  if (L==1) {
    b <- q[num_time_points]  - q[1] - mu*dT*(num_time_points-2*sum(q[1:num_time_points-1]))
    s_SL <- b/(integrated_var+1)
  } else{
    b <- q[num_time_points,]  - q[1,] - mu*dT*(array(num_time_points-1,c(1,L))-2*apply(q[1:num_time_points-1,],2,sum))

    s_SL <- (1/(integrated_var+1))*t(b)
  }
  s_SL_mat[ite_sim,]= t(s_SL)

  end_time <- Sys.time()
  array_SL_time[ite_sim] <-end_time-start_time

  ################

  # The below code repeats the first iteration due to weird timing problem. Uncomment out if timing problem persists.
  if (flag==1 & ite_sim==1){
    flag<-0
    ite_sim <- 0
  }

  ite_sim <- ite_sim+1
  print(ite_sim)


}

stext <- array(0,c(1,L))
for (ind_residue in 1:L){
  stext[ind_residue] <-paste("s",ind_residue-1,sep="")

}

info_combine <- cbind(array_load_time,array_load_time_onetwo+array_LLS_time,s_LLS_mat)
info_combine_df <-  data.frame(info_combine)
colnames(info_combine_df) <- c("loadtime","runtime",stext)
write.csv(info_combine_df, file = paste(paste("LLS", test_name, sep="_"), ".csv", sep=""))
#
#
#
info_combine <- cbind(array_load_time,array_load_time_onetwo+array_MPL_time,s_MPL_mat)
info_combine_df <-  data.frame(info_combine)
colnames(info_combine_df) <- c("loadtime","runtime",stext)
write.csv(info_combine_df, file = paste(paste("MPL", test_name, sep="_"), ".csv", sep=""))
#
#
info_combine <- cbind(array_load_time,array_load_time_onetwo+array_SL_time,s_SL_mat)
info_combine_df <-  data.frame(info_combine)
colnames(info_combine_df) <- c("loadtime","runtime",stext)
write.csv(info_combine_df, file = paste(paste("SL", test_name, sep="_"), ".csv", sep=""))

# info_combine_df <-  data.frame(array_load_time)
# colnames(info_combine_df) <- paste("loadtime",stext,sep=" ")
# write.csv(info_combine_df, file = paste(test_name,"loadtime.csv",sep=""))
