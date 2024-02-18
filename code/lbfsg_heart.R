library(dplyr)
library(lbfgs)
library(fastDummies)
# The "goal" field refers to the presence of heart disease
# in the patient.  It is integer valued from 0 (no presence) to 4.
# Experiments with the Cleveland database have concentrated on simply
# attempting to distinguish presence (values 1,2,3,4) from absence (value 0). 
heart.dat <- read.csv("~/Desktop/Research/HeartDisease/processed.cleveland.data", header = FALSE)
head(heart.dat)

col_names <- c("age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
               "thalach", "exang", "oldpeak", "slope", "ca", "thal", "target")

colnames(heart.dat) <- col_names
head(heart.dat)

table(heart.dat$target)

dim(heart.dat)

str(heart.dat)

# Create the histogram
ggplot(heart.dat, aes(x = target)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of Target",
       x = "Target",
       y = "Frequency") +
  theme_minimal()


# Use mutate and recode to group "1", "2", "3", and "4" as "1"
# heart_data <- heart.dat %>%
#   mutate(target = case_when(
#     target %in% c("2", "3", "4") ~ "1",
#     TRUE ~ as.character(target)
#   ))

heart_data <- heart.dat
round( (table(heart_data$target) / length(heart_data$target))*100, 2)

# Assuming your data frame is named 'heart_data'
# heart_data$sex <- factor(heart_data$sex)
heart_data$cp <- factor(heart_data$cp)
heart_data$fbs <- factor(heart_data$fbs)
# heart_data$exang <- factor(heart_data$exang)
heart_data$slope <- factor(heart_data$slope)
heart_data$thal <- factor(heart_data$thal)
heart_data$ca <- as.numeric(heart_data$ca)

str(heart_data)
colSums(is.na(heart_data))

# cat_names <- names(heart_data)
# heartcat[cat_names] <- lapply(heartcat[cat_names] , as.factor)
# 
# heartcat1 <- dummy_cols(heartcat, remove_first_dummy = TRUE)


# check MSE for all 100 dataset dim(replicated_data_train)[1]
kpr_simulation4both_1 <- function(data){ 
  # MSE for replicates
  replicates_min1 <- matrix(nrow = 0, ncol = 4)
  # Assign column names
  colnames(replicates_min1) <- c("MSE", "LLhood", "Sigma", "lambda")
  
  # MSE for replicates
  replicates_max1 <- matrix(nrow = 0, ncol = 4)
  # Assign column names
  colnames(replicates_max1) <- c("MSE", "LLhood", "Sigma", "lambda")
  
  train_mat1 <- data
  
  # Using poisson glm for alpha
  {
    train_df1 <- as.data.frame(train_mat1)
    
    # Fit a Poisson regression model
    poisson_model1 <- glm(y ~ x, data = train_df1, family = "poisson")
    
    # using intercept as b
    kpr_b01 <- as.numeric(poisson_model1$coefficients[1])
    
    # using intercept as b
    kpr_b11 <- as.numeric(poisson_model1$coefficients[2])
    
    # X * beta1 for alpha
    kpr_alpha1 <- train_mat1[, "x"] * kpr_b11
    
    kpr_b1 <- kpr_b01
  }
  
  ## Check for range of both lambda and sigma.
  {
    # check for both sigma and lambda values
    # Create an empty dataframe to store the results
    param_mat1 <- matrix(nrow = 0, ncol = 2)
    # Assign column names
    colnames(param_mat1) <- c("sigma", "lambda")
    
    #### Using GLM
    # X * beta1 for alpha
    kpr_alpha1 <- train_mat1[, "x"]  * kpr_b11
    kpr_b1 <- kpr_b01
    
    # sequence for lambda
    # exp(seq(-7, 2, by = 0.1))
    lambda_seq1 <- exp(seq(-6, 2, by = 0.1))
    
    # sequence for sigma
    # exp(seq(-8, 1, by = 0.01))
    kpr_sigma_seq1 <- exp(seq(-6, 1, by = 0.01))
    
    #for (kpr_lambda1 in lambda_seq1) {
    for (kpr_sigma1 in kpr_sigma_seq1) {
      for (kpr_lambda1 in lambda_seq1) {
        tryCatch({
          # kernel matrix
          kpr_k_mat1 <- rbf_kernel(train_mat1[, "x"] , train_mat1[, "x"] , sigma = kpr_sigma1)
          kpr_alpha_mat1 <- kpr_alpha1
          kpr_b_mat1 <- matrix(kpr_b1, nrow = dim(kpr_k_mat1)[1])
          
          # predicted mean for the ith fold data
          kpr_mu_hat1 <- exp( kpr_eta_funct(train_mat1[, "x"] , kpr_alpha1, kpr_b1, sigma = kpr_sigma1) )
          
          # first derivative
          kpr_frst_deri_mat11 <- kpr_k_mat1 %*%
            ( kpr_mu_hat1 - train_mat1[, "y"] +
                kpr_lambda1*kpr_alpha_mat1 )
          
          # for scaled y
          # kpr_frst_deri_mat11 <- kpr_k_mat1 %*%
          #   ( kpr_mu_hat1 - as.matrix(scale(train_mat1[, "y"])) +
          #       kpr_lambda1*kpr_alpha_mat1 )
          
          kpr_frst_deri_mat21 <- sum(kpr_mu_hat1 - train_mat1[, "y"])
          
          # for scaled y
          # kpr_frst_deri_mat21 <- sum(kpr_mu_hat1 - as.matrix(scale(train_mat1[, "y"])) )
          
          # kpr_frst_deri_mat21 <- sum(kpr_mu_hat1 - as.matrix( train_mat1[, "y"] ))
          
          # score matrix
          c11 <- dim(kpr_k_mat1)[1]
          c21 <- c11 + 1
          kpr_score_mat1 <- matrix(NA, nrow = c21, ncol = 1)
          kpr_score_mat1[1:c11, 1] <- kpr_frst_deri_mat11
          kpr_score_mat1[c21, 1] <- kpr_frst_deri_mat21
          
          # second derivative
          # first term (1,1)
          kpr_scnd_deri_mat11 <- kpr_k_mat1 %*% diag(c(kpr_mu_hat1)) %*% kpr_k_mat1 + kpr_lambda1*kpr_k_mat1
          # second term (1,2)
          kpr_scnd_deri_mat21 <- kpr_k_mat1 %*% kpr_mu_hat1
          # (2,1)
          kpr_scnd_deri_mat31 <- t(kpr_mu_hat1) %*% kpr_k_mat1
          # (2,2)
          kpr_scnd_deri_mat41 <- sum(kpr_mu_hat1)
          
          # hessian matrix
          kpr_hess_mat1 <- matrix(NA, nrow = c21, ncol = c21)
          kpr_hess_mat1[1:c11, 1:c11] <- kpr_scnd_deri_mat11
          kpr_hess_mat1[1:c11, c21] <- kpr_scnd_deri_mat21
          kpr_hess_mat1[c21, 1:c11] <- kpr_scnd_deri_mat31
          kpr_hess_mat1[c21, c21] <- c(kpr_scnd_deri_mat41)
          
          # predicted parameters
          kpr_alpha_b_mat1 <- matrix(c(kpr_alpha_mat1, kpr_b1), ncol = 1)
          kpr_alpha_b_mat_new1 <- kpr_alpha_b_mat1 - solve(kpr_hess_mat1) %*% kpr_score_mat1
          
          # Check if the Hessian matrix is invertible
          if (!is.null(solve(kpr_hess_mat1))) {
            
            # Store the values in the dataframe
            param_mat1 <- rbind(param_mat1, c(kpr_sigma1, kpr_lambda1))
          }
          
        }, error = function(e) {
          # print(paste("Error for sigma =", kpr_sigma1, "& lambda =", kpr_lambda1,  ":", e$message))
        })
      }
    }
    param_mat1
  }
  
  # Skip to the next iteration if condition is met
  if (nrow(param_mat1) > 0) {
    
    ### For loop for all combination of sigma and lambda
    results_list1 <- list()
    # list for final parameters
    parameters_list1 <- list()
    
    # Iterate over different pairs of sigma and lambda values
    for (k in 1:nrow(param_mat1)) {
      tryCatch(
        {
          kpr_sigma1 <- as.numeric(param_mat1[k, ][1])
          kpr_lambda1 <- as.numeric(param_mat1[k, ][2])
          
          results_mat1 <-  matrix(nrow = 0, ncol = 5)
          colnames(results_mat1) <- c("LLhood", "Diff", "para_diff", "Sigma", "lambda")
          
          #### Using GLM
          # X * beta1 for alpha
          kpr_alpha1 <- train_mat1[, "x"]  * kpr_b11
          kpr_b1 <- kpr_b01
          
          # likelihood
          kpr_llhood_old1 <- kpr_llhood(x = train_mat1[, "x"],
                                        y = train_mat1[, "y"],
                                        alpha = kpr_alpha1,
                                        b = kpr_b1,
                                        lambda = kpr_lambda1,
                                        sigma = kpr_sigma1)
          llhood_diff1 <- 1
          para_dif1 <- 1
          i <- 0
          
          ########### while loop ###########
          while(llhood_diff1 > 0.001 && para_dif1 > 0.01) {
            
            i <- i + 1
            
            # kernel matrix
            kpr_k_mat1 <- rbf_kernel(train_mat1[, "x"] , train_mat1[, "x"] , sigma = kpr_sigma1)
            kpr_alpha_mat1 <- matrix(kpr_alpha1, nrow = dim(kpr_k_mat1)[1])
            kpr_b_mat1 <- matrix(kpr_b1, nrow = dim(kpr_k_mat1)[1])
            
            # predicted mean for the ith fold data
            kpr_mu_hat1 <- exp( kpr_eta_funct(train_mat1[, "x"] , kpr_alpha1, kpr_b1, sigma = kpr_sigma1) )
            
            # first derivative
            kpr_frst_deri_mat11 <- kpr_k_mat1 %*% 
              ( kpr_mu_hat1 - train_mat1[, "y"] + 
                  kpr_lambda1*kpr_alpha_mat1 )
            
            # for scaled y
            # kpr_frst_deri_mat11 <- kpr_k_mat1 %*% 
            #                       ( kpr_mu_hat1 - as.matrix(scale(train_mat1[, "y"])) + 
            #                         kpr_lambda1*kpr_alpha_mat1 )
            
            kpr_frst_deri_mat21 <- sum(kpr_mu_hat1 - train_mat1[, "y"])
            
            # for scaled y
            # kpr_frst_deri_mat21 <- sum(kpr_mu_hat1 - as.matrix(scale(train_mat1[, "y"])) )
            
            # score matrix
            c11 <- dim(kpr_k_mat1)[1]
            c21 <- c11 + 1
            kpr_score_mat1 <- matrix(NA, nrow = c21, ncol = 1)
            kpr_score_mat1[1:c11, 1] <- kpr_frst_deri_mat11
            kpr_score_mat1[c21, 1] <- kpr_frst_deri_mat21
            
            # second derivative
            # first term (1,1)
            kpr_scnd_deri_mat11 <- kpr_k_mat1 %*% diag(c(kpr_mu_hat1)) %*% kpr_k_mat1 + kpr_lambda1*kpr_k_mat1
            # second term (1,2)
            kpr_scnd_deri_mat21 <- kpr_k_mat1 %*% kpr_mu_hat1
            # (2,1)
            kpr_scnd_deri_mat31 <- t(kpr_mu_hat1) %*% kpr_k_mat1
            # (2,2)
            kpr_scnd_deri_mat41 <- sum(kpr_mu_hat1)
            
            # hessian matrix
            kpr_hess_mat1 <- matrix(NA, nrow = c21, ncol = c21)
            kpr_hess_mat1[1:c11, 1:c11] <- kpr_scnd_deri_mat11
            kpr_hess_mat1[1:c11, c21] <- kpr_scnd_deri_mat21
            kpr_hess_mat1[c21, 1:c11] <- kpr_scnd_deri_mat31
            kpr_hess_mat1[c21, c21] <- c(kpr_scnd_deri_mat41)
            
            # predicted parameters
            kpr_alpha_b_mat1 <- matrix(c(kpr_alpha_mat1, kpr_b1), ncol = 1)
            kpr_alpha_b_mat_new1 <- kpr_alpha_b_mat1 - solve(kpr_hess_mat1) %*% kpr_score_mat1
            
            # para diff
            para_dif1 <- sum(abs(kpr_alpha_b_mat_new1 - kpr_alpha_b_mat1))
            
            kpr_alpha1 <- kpr_alpha_b_mat_new1[1:c11]
            kpr_b1 <- kpr_alpha_b_mat_new1[c21]
            
            # new likelihood
            kpr_llhood_new1 <- kpr_llhood(x = train_mat1[, "x"] ,
                                          y = train_mat1[, "y"],
                                          alpha = kpr_alpha1,
                                          b = kpr_b1,
                                          lambda = kpr_lambda1,
                                          sigma = kpr_sigma1)
            
            llhood_diff1 <- abs(kpr_llhood_new1 - kpr_llhood_old1)
            
            # update likelihood
            kpr_llhood_old1 <- kpr_llhood_new1
            
            # update results
            new_param1 = c(kpr_llhood_new1 ,llhood_diff1, para_dif1, kpr_sigma1, kpr_lambda1)
            results_mat1 = rbind(results_mat1, c(new_param1) )
            
            # results list
            results_list1[[k]] <- round(results_mat1, 5)
            
            # parameter list
            parameters_list1[[k]] <- kpr_alpha_b_mat_new1
            #print(round(results_mat1, 5) )
          }
          
          #
        },
        error = function(e){ return(NA) })
    }
    # results_list1
    
    # prediction on training set
    final_table1 <- matrix(nrow = 0, ncol = 4)
    # Assign column names
    colnames(final_table1) <- c("Sigma", "lambda", "LLhod", "Train_MSE")
    
    
    for (m in 1:length(results_list1) ){
      
      true_mean1 <- eta_true(x = train_mat1[, "x"] , z = train_mat1[, "z"])  
      
      best_kpr_alpha1 <- parameters_list1[[m]][1:c11]
      best_kpr_b1 <- parameters_list1[[m]][c21]
      best_kpr_sigma1 <- results_list1[[m]][,"Sigma"][1]
      best_kpr_lambda1 <- results_list1[[m]][,"lambda"][1]
      
      # predicted mean for the ith fold data
      kpr_mu_hat1 <- exp(kpr_eta_funct(x = train_mat1[, "x"],
                                       alpha = best_kpr_alpha1,
                                       b = best_kpr_b1,
                                       sigma = best_kpr_sigma1)
      )
      
      kpr_mse1 <- mean((true_mean1 - kpr_mu_hat1)^2)
      
      final_llhood1 <- kpr_llhood(x = train_mat1[, "x"],
                                  y = train_mat1[, "y"],
                                  alpha = best_kpr_alpha1,
                                  b = best_kpr_b1,
                                  lambda = best_kpr_lambda1,
                                  sigma = best_kpr_sigma1)
      
      # Add some rows to the dataframe
      new_rows1 <- c(best_kpr_sigma1, best_kpr_lambda1, final_llhood1, kpr_mse1)
      
      final_table1 <- rbind(final_table1, c(new_rows1))
    }
    # final_table1
    
    # based on maximum llhood
    batch_sigma1 <- as.numeric(final_table1[which.max(final_table1[, "LLhod"]), ][1])
    # batch lambda
    batch_lambda1 <- as.numeric(final_table1[which.max(final_table1[, "LLhod"]), ][2]) 
    # batch llhood
    batch_llhood1 <- as.numeric(final_table1[which.max(final_table1[, "LLhod"]), ][3])
    # mse for data batch
    batch_mse1 <- as.numeric(final_table1[which.max(final_table1[, "LLhod"]), ][4])
    # Add some rows to the dataframe
    new_max1 <- c(batch_mse1, batch_llhood1, batch_sigma1, batch_lambda1)
    replicates_max1 <- rbind(replicates_max1, c(new_max1))
    
    # based on minimum MSE
    batch_sigma4min <- as.numeric(final_table1[which.min(final_table1[, "Train_MSE"]), ][1])
    # batch lambda
    batch_lambda4min <- as.numeric(final_table1[which.min(final_table1[, "Train_MSE"]), ][2]) 
    # batch llhood
    batch_llhood4min <- as.numeric(final_table1[which.min(final_table1[, "Train_MSE"]), ][3])
    # mse for data batch
    batch_mse4min <- as.numeric(final_table1[which.min(final_table1[, "Train_MSE"]), ][4])
    # Add some rows to the dataframe
    new_min1 <- c(batch_mse4min, batch_llhood4min, batch_sigma4min, batch_lambda4min)
    replicates_min1 <- rbind(replicates_min1, c(new_min1))
    
  }
  else {
    # If there are no missing values, return an empty matrix for this iteration
    new_mse1 <- c(NA, NA, NA, NA)
    replicates_min1 <- rbind(replicates_min1, c(new_mse1))
    replicates_max1 <- rbind(replicates_max1, c(new_mse1))
  }
  
  # final1 <- list()
  # final1[[1]] <- replicates_min1
  # final1[[2]] <- replicates_max1
  
  df1 <- as.data.frame(replicates_min1)
  df2 <- as.data.frame(replicates_max1)
  df_list <- c(df_1 = df1, df_2 = df2)
  
  # results1 <- list(max_llhood = replicates_max1, min_mse = replicates_min1)
  return(df_list)
}


heart_df1 <- heart_data[c("target", "age", "sex", "exang")]
heart_df1.cols <- c("y", "x", "v1", "v2")
colnames(heart_df1) <- heart_df1.cols
ht.mat1 <- as.matrix(heart_df1)
ht.mat1

kpr_simulation4both_1(ht.mat1)

{
  # Fit a Poisson regression model
  poi_mod <- glm(y ~ x, data = heart_df1, family = "poisson")
  
  # using intercept as b
  kpr_b0 <- as.numeric(poi_mod$coefficients[1])
  
  # using intercept as b
  kpr_b1 <- as.numeric(poi_mod$coefficients[2])
  
  # X * beta1 for alpha
  kpr_alpha <- ht.mat1[, "x"] * kpr_b1
  
  kpr_b <- kpr_b0
}

predicted <- predict(poi_mod, newdata = heart_df1, type = "response")
# Calculate the Mean Squared Error (MSE)
mse0 <- mean((heart_df1$y - predicted)^2)
mse0



###########################################################
############# Using LBFSG. ################################
###########################################################


# kpr gradient
kpr_gradient <- function(kpr.params, data, sigma.val = 1, lambda.val = 1){
  cc1 <- as.numeric(nrow(data))
  cc2 <- cc1 + 1
  
  # alpha.mat <- matrix(kpr.params[1:cc1], ncol = 1)
  alpha.mat <- kpr.params[1:cc1]
  b.not <- kpr.params[cc2]
  
  x <- data[, "x"]
  y <- data[, "y"]
  
  # kernel matrix
  kernel.mat <- rbf_kernel(x, x, sigma = sigma.val)
  
  b.not.mat <- matrix(b.not, nrow = cc1)
  
  # predicted mean 
  pred_mu <- exp( kpr_eta_funct(x, alpha.mat, b.not.mat, sigma = sigma.val) )
  
  # first derivative
  alpha.deri <- kernel.mat %*% ( pred_mu - y + lambda.val*alpha.mat )
  
  b.deri <- sum(pred_mu - y)
  
  # score matrix
  # score_mat <- matrix(NA, nrow = cc2, ncol = 1)
  # score_mat[1:cc1, 1] <- alpha.deri
  # score_mat[cc2, 1] <- b.deri
  # score_mat
  c(alpha.deri, b.deri)
  
}

# kpr llhood
kpr_llhood.new <- function(kpr.params, data, sigma.val, lambda.val){
  nn2 <- as.numeric(nrow(data))
  nn3 <- nn2 + 1
  
  alpha.mat <- kpr.params[1:nn2]
  b.not <- kpr.params[nn3]
  
  x <- data[, "x"]
  y <- data[, "y"]
  
  kpr.eta <- kpr_eta_funct(x, alpha.mat, b.not, sigma = sigma.val)
  
  kernel.mat <- rbf_kernel(x, x, sigma = sigma.val)
  
  mean( exp(kpr.eta) - y*kpr.eta ) + (lambda.val/2) * t(alpha.mat) %*% kernel.mat %*% alpha.mat
}


kpr_gradient(c(kpr_alpha, kpr_b0), ht.mat1, 1, 1)

kpr_llhood.new(c(kpr_alpha, kpr_b0), ht.mat1, 1, 1)


initial_params <- c(kpr_alpha, kpr_b0)
initial_params <- c( rep(0, nrow(ht.mat1) + 1 ) )

# Optimization using lbfgs
lbfgs(kpr_llhood.new,
      kpr_gradient,
      data = ht.mat1,
      sigma.val = 2, 
      lambda.val = 1,
      initial_params
)


rezult <- optim(initial_params,
                kpr_llhood.new,
                kpr_gradient,
                method = "L-BFGS-B",
                ht.mat1,
                sigma.val = 0.1, 
                lambda.val = 0.00005
)

# predicted mean for the ith fold data
pred <- exp(kpr_eta_funct(x = ht.mat1[, "x"],
                          alpha = rezult$par[1:303],
                          b = rezult$par[304],
                          sigma = 0.1)
)

mean((ht.mat1[, "y"] - pred)^2)


######################################################################
#####################      PLK      ##################################
######################################################################
# plk gradient
plk_gradient <- function(plk.params, data, sigma.val, lambda.val){
  x <- data[, "x"]
  y <- data[, "y"]
  v_s <- data[, c("v1", "v2")]
  
  cat.cols <- ncol(v_s)
  cc1 <- as.numeric(nrow(data))
  cc2 <- cc1 + 1
  cc3 <- cc1 + cat.cols
  cc4 <- cc3 + 1
  
  
  # alpha.mat <- matrix(kpr.params[1:cc1], ncol = 1)
  alpha.mat <- plk.params[1:cc1]
  beta.mat <- plk.params[cc2:cc3]
  b.not <- plk.params[cc4]
  
  # kernel matrix
  kernel.mat <- rbf_kernel(x, x, sigma = sigma.val)
  b.not.mat <- matrix(b.not, nrow = cc1)
  
  # predicted mean 
  pred_mu <- exp(plk_eta_funct(x = x,
                               z = v_s,
                               alpha = alpha.mat,
                               beta = beta.mat,
                               b = b.not.mat,
                               sigma = sigma.val)
  )
  
  # first derivative
  alpha.deri <- kernel.mat %*% (pred_mu - y + lambda.val*alpha.mat)
  beta.deri <- t(v_s) %*% (pred_mu - y) 
  b.deri <- sum(pred_mu - y) 
  c(alpha.deri, beta.deri, b.deri)
  
}

# pl llhood
plk_llhood.new <- function(plk.params, data, sigma.val, lambda.val){
  x <- data[, "x"]
  y <- data[, "y"]
  v_s <- data[, c("v1", "v2")]
  
  cat.cols <- ncol(v_s)
  cc1 <- as.numeric(nrow(data))
  cc2 <- cc1 + 1
  cc3 <- cc1 + cat.cols
  cc4 <- cc3 + 1
  
  alpha.mat <- plk.params[1:cc1]
  beta.mat <- plk.params[cc2:cc3]
  b.not <- plk.params[cc4]
  b.not.mat <- matrix(b.not, nrow = cc1)
  
  plk.eta <- plk_eta_funct(x = x,
                           z = v_s,
                           alpha = alpha.mat,
                           beta = beta.mat,
                           b = b.not.mat,
                           sigma = sigma.val)
  
  kernel.mat <- rbf_kernel(x, x, sigma = sigma.val)
  mean( exp(plk.eta) - y*plk.eta ) + (lambda.val/2) * t(alpha.mat) %*% kernel.mat %*% alpha.mat
  
}

plk_initial_params <- c( rep(0.1, nrow(ht.mat1) + 3 ) )

# Optimization using lbfgs
lbfgs(plk_llhood.new,   
      plk_gradient,
      data = ht.mat1,
      sigma.val = 2, 
      lambda.val = 1,
      plk_initial_params
)


rezult1 <- optim(plk_initial_params,
                 plk_llhood.new,
                 plk_gradient,
                 method = "L-BFGS-B",
                 ht.mat1,
                 sigma.val = 1,
                 lambda.val = 1)

# predicted mean for the ith fold data
plk_pred <- exp(plk_eta_funct(x = ht.mat1[, "x"],
                              z = ht.mat1[, c("v1", "v2")],
                              alpha = rezult1$par[1:303],
                              beta = rezult1$par[304:305],
                              b = rezult1$par[306],
                              sigma = 1)
)

mean((ht.mat1[, "y"] - plk_pred)^2)




######################################################################
#####################      PKP      ##################################
######################################################################

pkp_gradient <- function(pkp.params, data, sigma.val = 1, beta.val = 1, lambda.val = 1){
  cc1 <- as.numeric(nrow(data))
  cc2 <- cc1 + 1
  
  alpha.mat <- pkp.params[1:cc1]
  b.not <- pkp.params[cc2]
  
  x <- data[, "x"]
  y <- data[, "y"]
  z <- data[, "v1"]
  
  # kernel matrix
  kernel.mat <- hamin_dist(x = x,
                           y = y,
                           z1 = z,
                           z2 = z,
                           sigma = sigma.val,
                           beta = beta.val
  )
  
  b.not.mat <- matrix(b.not, nrow = cc1)
  
  # predicted mean 
  pred_mu <- exp(pkp_eta_funct.new(x, y, z, alpha.mat, b.not, sigma = sigma.val, beta = beta.val))
  
  # first derivative
  alpha.deri <- kernel.mat %*% ( pred_mu - y + lambda.val*alpha.mat )
  b.deri <- sum(pred_mu - y)
  c(alpha.deri, b.deri)
  
}

pkp_llhood.new <- function(pkp.params, data, sigma.val = 1, beta.val = 1, lambda.val = 1){
  nn2 <- as.numeric(nrow(data))
  nn3 <- nn2 + 1
  
  alpha.mat <- pkp.params[1:nn2]
  b.not <- pkp.params[nn3]
  
  x <- data[, "x"]
  y <- data[, "y"]
  z <- data[, "v1"]
  
  pkp.eta <- pkp_eta_funct.new(x, y, z, alpha.mat, b.not, sigma = sigma.val, beta = beta.val)
  
  kernel.mat <- hamin_dist(x = x, y = y, z1 = z, z2 = z,
                           sigma = sigma.val, beta = beta.val)
  
  mean( exp(pkp.eta) - y*pkp.eta ) + (lambda.val/2) * t(alpha.mat) %*% kernel.mat %*% alpha.mat
}

#########################################################
pkp_initial_params <- c( rep(0.1, nrow(ht.mat1) + 1 ) )

{
  pkp.sig <- 10
  pkp.bet <- 1
  pkp.lab <- 0.1
}

# Optimization using lbfgs
pkp.lbf <- lbfgs(pkp_llhood.new,   
                 pkp_gradient,
                 data = ht.mat1,
                 sigma.val = pkp.sig, 
                 beta.val = pkp.bet,
                 lambda.val = pkp.lab,
                 pkp_initial_params
)

optim(par = pkp_initial_params,
      fn = pkp_llhood.new,
      gr = pkp_gradient,
      method = "L-BFGS-B",
      ht.mat1,
      sigma.val = pkp.sig,
      beta.val = pkp.bet,
      lambda.val = pkp.lab)

rezult2 <- optim(pkp_initial_params,
                 pkp_llhood.new,
                 pkp_gradient,
                 method = "L-BFGS-B",
                 ht.mat1,
                 sigma.val = pkp.sig,
                 beta.val = pkp.bet,
                 lambda.val = pkp.lab)

# predicted mean for the ith fold data
pkp_pred <- exp(pkp_eta_funct.new(x =  ht.mat1[, "x"],
                                  y =  ht.mat1[, "x"],
                                  z =  ht.mat1[, "v1"],
                                  alpha = rezult2$par[1:303],
                                  b = rezult2$par[304],
                                  sigma = pkp.sig,
                                  beta = pkp.bet) )

pkp_pred1 <- exp(pkp_eta_funct.new(x =  ht.mat1[, "x"],
                                   y =  ht.mat1[, "x"],
                                   z =  ht.mat1[, "v1"],
                                   alpha = pkp.lbf$par[1:303],
                                   b = pkp.lbf$par[304],
                                   sigma = pkp.sig,
                                   beta = pkp.bet) )

mean((ht.mat1[, "y"] - pkp_pred1)^2)
