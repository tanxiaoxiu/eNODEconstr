library(deSolve)
library(gtools)
library(tidyverse)
library(readxl)
library(pheatmap)
library(reshape2)
library(rhdf5)
library(dplyr)
library(tidyr)


#Definition consumer_resource_model
consumer_resource_model <- function(t, state, parms) {
  C <- state[1:parms$N]
  R <- state[(parms$N + 1):(parms$N + parms$M)]
  mu <- matrix(parms$mu, nrow = parms$N, ncol = parms$M)
  m <- parms$m
  rho <- parms$rho
  omega <- parms$omega
  l <- matrix(parms$l, nrow = parms$M, ncol = parms$M)
  lambda <- parms$lambda
  dCdt <- rep(0, parms$N)
  dRdt <- rep(0, parms$M)
  dCdt <- C * (R %*% t(mu) * (1 - lambda)) - C * m
  dRdt <- rho - R * omega - (C %*% mu) * R + lambda * ((C %*% mu) * R) %*% t(l)
  list(c(dCdt, dRdt))
}

#mu
create_consumption_matrix <- function(N, M, nonzero_consumption) {
  total_elements <- N * M
  num_nonzero <- floor(total_elements * nonzero_consumption)
  all_elements <- c(runif(num_nonzero, min = 0, max = 5), rep(0, total_elements - num_nonzero))
  all_elements <- sample(all_elements)
  consumption_matrix <- matrix(all_elements, nrow = N, ncol = M, byrow = TRUE)
  for (i in 1:N) {
    non_zero_elements <- which(consumption_matrix[i, ] != 0)
    non_zero_count <- length(non_zero_elements)
    if (non_zero_count > 0) {
      consumption_matrix[i, non_zero_elements] <- consumption_matrix[i, non_zero_elements] / non_zero_count
    }
  }

  return(consumption_matrix)
}

#l
generate_byproduct_matrix <- function(M, nonzero_byproduct) {
  num_elements <- M * M
  num_nonzeros <- floor(num_elements * nonzero_byproduct)
  values <- c(runif(num_nonzeros, 0, 1), rep(0, num_elements - num_nonzeros))
  D <- matrix(sample(values), nrow = M, ncol = M)
  col_sums <- colSums(D)
  col_sums[col_sums == 0] <- 1
  D_norm <- sweep(D, 2, col_sums, FUN = "/")
  return(D_norm)
}


run_simulation_and_save <- function(S, N, M) {
  base_work_dir <- sprintf("~/eNODEconstr/simulation/generation_data/n%dm%d/s%d", N, M, S)
  
  if (!dir.exists(base_work_dir)) {
    dir.create(base_work_dir, recursive = TRUE)
  }
  
  timepoints <- list(
    t5 = c(0, 2, 6, 14, 24),
    t10 = c(0, 1, 2, 6, 9, 12, 15, 18, 21, 24),
    t15 = c(0, 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 19, 21, 24),
    t20 = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24),
    t25 = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)
  )
  
  run_simulation_for_timepoint_and_iter <- function(timepoint_name, timepoint_values) {
    timepoint_dir <- sprintf("%s/%s", base_work_dir, timepoint_name)
    if (!dir.exists(timepoint_dir)) {
      dir.create(timepoint_dir, recursive = TRUE)
    }
    
    for (iter in 1:10) {
      iter_dir <- sprintf("%s/iter%d", timepoint_dir, iter)
      if (!dir.exists(iter_dir)) {
        dir.create(iter_dir, recursive = TRUE)
      }
      
      set.seed(iter)
      lambda <- 0.5
      m <- rep(0.2, N)
      rho <- runif(M)
      omega <- rep(0.2, M)
      nonzero_consumption <- 0.2
      mu <- create_consumption_matrix(N, M, nonzero_consumption)
      nonzero_byproduct <- 0.5
      l <- generate_byproduct_matrix(M, nonzero_byproduct)
      parameter_matrix <- list(N = N, M = M, mu = mu, m = m, rho = rho, omega = omega, l = l, lambda = lambda)
      
      seeds <- 1:S
      results <- lapply(seeds, function(seed) {
        set.seed(iter * 1000 + seed) 
        state <- c(runif(N), runif(M))
        parameters <- list(N = N, M = M, mu = mu, m = m, rho = rho, omega = omega, l = l, lambda = lambda)
        times <- seq(0, 50, by = 1)
        
        out <- ode(y = state, times = times, func = consumer_resource_model, parms = parameters, method = "ode45")
        times_matrix <- out[, 1]
        
        sparse_indices <- which(times_matrix %in% timepoint_values)
        out_sparse <- out[sparse_indices, ]
        
        new_out_list <- list()
        new_out_sparse_list <- list()
        
        for (i in 1:N) {
          new_N <- N
          new_mu <- parameters$mu
          new_m <- parameters$m
          new_state <- state
          new_state[i] <- 0 
          
          new_parameters <- list(N = new_N, M = M, mu = new_mu, m = new_m, 
                                 rho = parameters$rho, omega = parameters$omega, 
                                 l = parameters$l, lambda = parameters$lambda)
          
          new_out <- ode(y = new_state, times = times, func = consumer_resource_model, parms = new_parameters, method = "ode45")
          
          new_out_sparse <- new_out[sparse_indices, ]
          
          new_out_list[[i]] <- new_out
          new_out_sparse_list[[i]] <- new_out_sparse
          
        }
        
        return(list(out = out, out_sparse = out_sparse,new_out_list = new_out_list,new_out_sparse_list=new_out_sparse_list,parameter_matrix=parameter_matrix))
      })   
      
      saveRDS(results, sprintf("%s/results.rds", iter_dir))
    }
  }
  
  for (timepoint_name in names(timepoints)) {
    run_simulation_for_timepoint_and_iter(timepoint_name, timepoints[[timepoint_name]])
  }
  
  return(TRUE)
}


S_values <- c(10, 15, 20)
N <- 10
M <- 10 
lapply(S_values, function(S) run_simulation_and_save(S, N, M))



###################################################Save data
library(tidyr) 
library(rhdf5)

NM_combinations <- list(c(10, 10))
S_values <- c(10, 15, 20)
T_values <- c(5, 10, 15, 20, 25)
iter_values <- 1:10

for (NM in NM_combinations) {
  N <- NM[1]
  M <- NM[2]
  
  for (S in S_values) {
    for (t in T_values) {
      for (iter in iter_values) {
        work_dir <- sprintf("~/eNODEconstr/simulation/generation_data/n%dm%d/s%d/t%d/iter%d", N, M, S, t, iter)
        setwd(work_dir)
        
        if (file.exists("results.rds")) {
          results <- readRDS("results.rds")
          
          sparse_out_list <- lapply(results, function(x) x$out_sparse)
          sparse_out_list <- lapply(sparse_out_list, function(x) {
            x <- x[, -which(colnames(x) == "time"), drop = FALSE]
            return(x)
          })
          saveRDS(sparse_out_list, "true_initial_state.rds")
          
          true_initial_state <- simplify2array(sparse_out_list)
          h5save(file = "true_initial_state.h5", true_initial_state)
          
          sparse_new_out_list <- lapply(results, function(simulation) {
            lapply(simulation$new_out_sparse_list, function(data) {
              data[ , -which(colnames(data) == "time"), drop = FALSE]
            })
          })
          saveRDS(sparse_new_out_list, "true_new_state.rds")
          
          mu <- lapply(results, function(x) x$parameter_matrix$mu)
          mu <- mu[[1]]
          write.table(mu, "mu_true.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
          
          mu_binary <- mu
          mu_binary[mu_binary > 0] <- 1
          write.table(mu_binary, "mu_binary_matrix.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
          
          l <- lapply(results, function(x) x$parameter_matrix$l)
          l <- l[[1]]
          write.table(l, "l_true.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
          
          rho <- lapply(results, function(x) x$parameter_matrix$rho)
          rho <- rho[[1]]
          write.table(rho, "rho_true.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
          
          m <- lapply(results, function(x) x$parameter_matrix$m)
          m <- m[[1]]
          write.table(m, "m_true.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
          
          omega <- lapply(results, function(x) x$parameter_matrix$omega)
          omega <- omega[[1]]
          write.table(omega, "omega_true.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
          
          lambda <- lapply(results, function(x) x$parameter_matrix$lambda)
          lambda <- lambda[[1]]
          write.table(lambda, "lambda_true.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
          
          
          ##################Score
          #true_initial_state
          true_initial_state <- readRDS("true_initial_state.rds")
          
          #true_new_state
          true_new_state <- readRDS("true_new_state.rds")
          
          calculate_L2_Distances <- function(true_initial_state, true_new_state, S, N, M) {
            
            results <- vector("list", S)
            
            for (s in 1:S) {
              true_initial_state_del <- true_initial_state[[s]]
              
              C0 <- true_initial_state_del[, 1:N] 
              R0 <- true_initial_state_del[, c((N+1):(N+M))] 
              
              true_new_state_subject <- true_new_state[[s]]
              
              C_distance <- numeric(N)
              R_distance <- numeric(M)
              
              for (n in 1:N) {
                
                true_new_state_subject_del <- true_new_state_subject[[n]]
                C1 <- true_new_state_subject_del[, 1:N] 
                R1 <- true_new_state_subject_del[, c((N+1):(N+M))] 
                
                # Calculate L2 distances
                C_distance[n] <- sqrt(sum((C0[, -n] - C1[,-n])^2))
                R_distance[n] <- sqrt(sum((R0 - R1)^2))
                
              }
              
              results[[s]] <- list(C_distance = C_distance, R_distance = R_distance)
            }
            
            return(results)
          }
          
          calculate_l2_results <-calculate_L2_Distances(true_initial_state, true_new_state, S, N, M)
          
          # Calculate average distances across all runs
          C_l2_d_avg <- colMeans(do.call(rbind, lapply(calculate_l2_results, `[[`, "C_distance")))
          R_l2_d_avg <- colMeans(do.call(rbind, lapply(calculate_l2_results, `[[`, "R_distance")))

          microbes <- paste0("Microbe", 1:N)
          C_l2_d_avg <- data.frame(Regulate = microbes, Score = C_l2_d_avg, Group = "Microbe")
          C_l2_d_avg$Score_normalized <- C_l2_d_avg$Score / sum(C_l2_d_avg$Score)
          
          R_l2_d_avg <- data.frame(Regulate = microbes, Score = R_l2_d_avg, Group = "Metabolite")
          R_l2_d_avg$Score_normalized <- R_l2_d_avg$Score / sum(R_l2_d_avg$Score)
          
          combined_df <- rbind(C_l2_d_avg, R_l2_d_avg)
          write.table(combined_df,file ="score_true_mean_l2.txt",row.names = F,col.names = TRUE, sep = "\t",quote = F)
          
          ###Microbe
          C_l2_d_list <- do.call(rbind, lapply(calculate_l2_results, `[[`, "C_distance"))
          microbes <- paste0("Microbe", 1:N)
          subjects <- paste0("Subject", 1:S)
          
          C_l2_normalized <- as.data.frame(sweep(C_l2_d_list, 1, rowSums(C_l2_d_list), FUN = "/"))
          colnames(C_l2_normalized) <- microbes
          C_l2_normalized <- cbind(Subject = subjects, C_l2_normalized)
          C_l2_normalized_long <- pivot_longer(C_l2_normalized, cols = microbes, names_to = "Regulate", values_to = "Score")
          C_l2_normalized_long$Group <- 'Microbe'
          
          C_l2 <- data.frame(C_l2_d_list)
          colnames(C_l2) <- microbes
          C_l2 <- cbind(Subject = subjects, C_l2)
          C_l2_long <- pivot_longer(C_l2, cols = microbes, names_to = "Regulate", values_to = "Score")
          C_l2_long$Group <- 'Microbe'
          
          ###Metabolite
          R_l2_d_list <- do.call(rbind, lapply(calculate_l2_results, `[[`, "R_distance"))
          R_l2_normalized <- as.data.frame(sweep(R_l2_d_list, 1, rowSums(R_l2_d_list), FUN = "/"))
          colnames(R_l2_normalized) <- microbes
          R_l2_normalized <- cbind(Subject = subjects, R_l2_normalized)
          R_l2_normalized_long <- pivot_longer(R_l2_normalized, cols = microbes, names_to = "Regulate", values_to = "Score")
          R_l2_normalized_long$Group <- 'Metabolite'
          
          R_l2 <- data.frame(R_l2_d_list)
          colnames(R_l2) <- microbes
          R_l2 <- cbind(Subject = subjects, R_l2)
          R_l2_long <- pivot_longer(R_l2, cols = microbes, names_to = "Regulate", values_to = "Score")
          R_l2_long$Group <- 'Metabolite'
          
          combined_score_l2 <- rbind(C_l2_long, R_l2_long)
          write.table(combined_score_l2,file ="score_true_subject_l2.txt",row.names = F,col.names = TRUE, sep = "\t",quote = F)
          
          combined_score_normalized_l2 <- rbind(C_l2_normalized_long, R_l2_normalized_long)
          write.table(combined_score_normalized_l2,file ="score_normalized_true_subject_l2.txt",row.names = F,col.names = TRUE, sep = "\t",quote = F)
          
          calculate_L1_Distances <- function(true_initial_state, true_new_state, S, N, M) {
            
            results <- vector("list", S)
            
            for (s in 1:S) {
              true_initial_state_del <- true_initial_state[[s]]
              
              C0 <- true_initial_state_del[, 1:N] 
              R0 <- true_initial_state_del[, c((N+1):(N+M))] 
              
              true_new_state_subject <- true_new_state[[s]]
              
              C_distance <- numeric(N)
              R_distance <- numeric(M)
              
              for (n in 1:N) {
                
                true_new_state_subject_del <- true_new_state_subject[[n]]
                C1 <- true_new_state_subject_del[, 1:N] 
                R1 <- true_new_state_subject_del[, c((N+1):(N+M))] 
                
                # Calculate L1 distances
                C_distance[n] <- sum(abs(C0[, -n] - C1[, -n]))
                R_distance[n] <- sum(abs(R0 - R1))
              }
              
              results[[s]] <- list(C_distance = C_distance, R_distance = R_distance)
            }
            
            return(results)
          }
          
          calculate_l1_results <-calculate_L1_Distances(true_initial_state, true_new_state, S, N, M)
  
          C_l1_d_avg <- colMeans(do.call(rbind, lapply(calculate_l1_results, `[[`, "C_distance")))
          R_l1_d_avg <- colMeans(do.call(rbind, lapply(calculate_l1_results, `[[`, "R_distance")))
          
          microbes <- paste0("Microbe", 1:N)
          C_l1_d_avg <- data.frame(Regulate = microbes, Score = C_l1_d_avg, Group = "Microbe")
          C_l1_d_avg$Score_normalized <- C_l1_d_avg$Score / sum(C_l1_d_avg$Score)
          
          R_l1_d_avg <- data.frame(Regulate = microbes, Score = R_l1_d_avg, Group = "Metabolite")
          R_l1_d_avg$Score_normalized <- R_l1_d_avg$Score / sum(R_l1_d_avg$Score)
          
          combined_df <- rbind(C_l1_d_avg, R_l1_d_avg)
          write.table(combined_df,file ="score_true_mean_l1.txt",row.names = F,col.names = TRUE, sep = "\t",quote = F)
          
          ##################
          ###Microbe
          C_l1_d_list <- do.call(rbind, lapply(calculate_l1_results, `[[`, "C_distance"))
          microbes <- paste0("Microbe", 1:N)
          subjects <- paste0("Subject", 1:S)
          
          C_l1_normalized <- as.data.frame(sweep(C_l1_d_list, 1, rowSums(C_l1_d_list), FUN = "/"))
          colnames(C_l1_normalized) <- microbes
          C_l1_normalized <- cbind(Subject = subjects, C_l1_normalized)
          C_l1_normalized_long <- pivot_longer(C_l1_normalized, cols = microbes, names_to = "Regulate", values_to = "Score")
          C_l1_normalized_long$Group <- 'Microbe'
          
          C_l1 <- data.frame(C_l1_d_list)
          colnames(C_l1) <- microbes
          C_l1 <- cbind(Subject = subjects, C_l1)
          C_l1_long <- pivot_longer(C_l1, cols = microbes, names_to = "Regulate", values_to = "Score")
          C_l1_long$Group <- 'Microbe'
          
          ###Metabolite
          R_l1_d_list <- do.call(rbind, lapply(calculate_l1_results, `[[`, "R_distance"))
          R_l1_normalized <- as.data.frame(sweep(R_l1_d_list, 1, rowSums(R_l1_d_list), FUN = "/"))
          colnames(R_l1_normalized) <- microbes
          R_l1_normalized <- cbind(Subject = subjects, R_l1_normalized)
          R_l1_normalized_long <- pivot_longer(R_l1_normalized, cols = microbes, names_to = "Regulate", values_to = "Score")
          R_l1_normalized_long$Group <- 'Metabolite'
          
          R_l1 <- data.frame(R_l1_d_list)
          colnames(R_l1) <- microbes
          R_l1 <- cbind(Subject = subjects, R_l1)
          R_l1_long <- pivot_longer(R_l1, cols = microbes, names_to = "Regulate", values_to = "Score")
          R_l1_long$Group <- 'Metabolite'
          
          combined_score_l1 <- rbind(C_l1_long, R_l1_long)
          write.table(combined_score_l1,file ="score_true_subject_l1.txt",row.names = F,col.names = TRUE, sep = "\t",quote = F)
          
          combined_score_normalized_l1 <- rbind(C_l1_normalized_long, R_l1_normalized_long)
          write.table(combined_score_normalized_l1,file ="score_normalized_true_subject_l1.txt",row.names = F,col.names = TRUE, sep = "\t",quote = F)
          
          save.image(file = "Save_data.RData")
          
        } else {
          cat("File 'results.rds' does not exist in", work_dir, "\n")
        }
      }
    }
  }
}
