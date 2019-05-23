library(MASS)

pseudo_inverse <- function(X) {
	# Single point for changing how we invert matrices.
	#output <- ginv(X)
	output <- solve(X)
	return (output)
}

generate_gamma <- function(d, c_num) {
	# Generate a random cholesky factor for correlation amongh variables. 
	d_elems <- runif(d, min = c_num, max = 0.1)
	d_matrix <- diag(d_elems)
	u_matrix <- matrix(rnorm(d*d), nrow=d, ncol=d)

	cov_matrix <- u_matrix %*% d_matrix %*% t(u_matrix)		
	
	gamma <- chol(cov_matrix)
	return (gamma)
}
	
generate_x <- function(n, gamma) {
	# Generate correlated values of X, d dimensions by n observations. The first column is a constant intercept value 1.0.
	d <- nrow(gamma)
	X <- matrix(rnorm((d-1)*n), nrow=n, ncol=(d-1))
	X <- cbind(rep(1, n), X)
	
	X <- X %*% gamma
	return (X)
}

generate_theta <- function(d, gamma) {
	X <- generate_x(10000, gamma)
	
	for(i in 1:100) {
		theta_val <- rnorm(d)		
		Y <- generate_y(X, theta_val, 10000)	
		
		# Make sure these thetas don't result in wildly unbalanced labels.	
		mu <- mean(Y)
		print(paste("Generating Theta[", i, "]: ", mu))
			
		if(mu > 0.35 && mu < 0.65) {
			return (theta_val)
		}
	}
	
	stop("Unable to generate good theta.")
}

# The next 6 functions are specific to logistic regression. 
# To try this with a different model, change these. 

logistic <- function(power_score) {
	output <- 1.0 / (1.0 + exp(-power_score))
	return (output)
}


y_means <- function(X, theta) {
	# computes E[y] for logistic model. 
	power_scores <- X %*% c(theta)
	prob <- logistic(power_scores)
	return (prob)
}

calc_entropy <- function(X, Y, theta) {
	# Computes E[entropy] for logistic model
	mean_y <- y_means(X, theta)
	v <- -log((Y * mean_y) + ((1 - Y) * (1 - mean_y)))
	return (mean(v))
}

gradient_matrix <- function(X, Y, theta) {
	# Computes pointwise gradients for logistic model
	mean_y <- y_means(X, theta)
	scale <- Y*(1.0 - mean_y) - (1.0 - Y)*mean_y
	sm <- matrix(rep(scale, ncol(X)), ncol=ncol(X), nrow=nrow(X))
	gradient_matrix <- X * -sm
	return (gradient_matrix)
}


hessian <- function(X, theta) {
	# Computes hessian matrix for logistic model
	mean_y <- y_means(X, theta)
	scale <- sqrt(mean_y * (1.0 - mean_y));
	sm <- matrix(rep(scale, ncol(X)), ncol=ncol(X), nrow=nrow(X))
	X1 <- sm * X
	
	h_matrix <- (t(X1) %*% X1) / nrow(X)
	return (h_matrix)
}


generate_y <- function(X, theta, n) {
	# Generates Y through random sampling for logistic model
	y_prob <- runif(n)
	y_mean <- y_means(X, theta)
	y <- y_prob < y_mean
	return (y)
}


# Everything after this should be unchanged from one model to the next

error_oos <- function(theta, true_theta, gamma, n) {
	X <- generate_x(n, gamma)
	Y <- generate_y(X, true_theta, n)
	entropy <- calc_entropy(X, Y, theta)
	return (entropy)
}

lambda_xval_en <- function(X, Y, lambda1, lambda2, target, mod) {
	# An elastic net cross-val objective function.
	test <- seq(target, nrow(X), mod)
	train <- seq(1, nrow(X), 1)
	train <- train[-test]
	
	X_train <- X[train,]
	X_test <- X[test,]
	Y_train <- matrix(Y[train], length(train), 1)
	Y_test <- matrix(Y[test], length(test), 1)

	fit_theta <- theta_en_solve(X_train, Y_train, lambda1, lambda2)
	oos1 <- sample_error(X_test, Y_test, fit_theta)
	return (oos1)
}

lambda_objective_ridge <- function(X, Y,lambda) {
	n_orig <- nrow(X)
	
	oos_total <- 0.0
	n_splits <- 10
	
	for(k in 1:n_splits) {
		oos_next <- lambda_xval_en(X, Y, 0.0, lambda, k, n_splits)		
		oos_total <- oos_total + oos_next
	}
	
	oos_total <- oos_total / n_splits
	return (oos_total)
}

lambda_objective_lasso <- function(X, Y,lambda) {
	n_orig <- nrow(X)
	
	oos_total <- 0.0
	n_splits <- 10
	
	for(k in 1:n_splits) {
		oos_next <- lambda_xval_en(X, Y, lambda, 0.0, k, n_splits)		
		oos_total <- oos_total + oos_next
	}
	
	oos_total <- oos_total / n_splits
	return (oos_total)
}

lambda_opt_approx_ridge <- function(X, Y) {	
	lambda_0 <- optimize(lambda_objective_ridge, X = X, Y = Y, c(0.0, 0.5))
	return (lambda_0$minimum)
    }

lambda_opt_approx_lasso <- function(X, Y) {
	lambda_0 <- optimize(lambda_objective_lasso, X = X, Y = Y, c(0.0, 0.5))
	return (lambda_0$minimum)
	}

sample_error <- function(X, Y, theta) {
	return (calc_entropy(X, Y, theta))
}

error_oos <- function(theta, true_theta, gamma, n) {
	X <- generate_x(n, gamma)
	Y <- generate_y(X, true_theta, n)
	return (sample_error(X, Y, theta))	
}

reg_objective <- function(theta, x, y, l1, l2) {
	se <- 2 * calc_entropy(x, y, theta)
	reg_term2 <- l2 * (t(theta) %*% theta)
	reg_term <- l1 * sum(abs(theta))
	
	output <- t(se) + reg_term + reg_term2
	return (output)
}

ice_objective <- function(theta, X, Y) {
	n <- nrow(X)
	entropy <- calc_entropy(X, Y, theta)

	h_matrix <- hessian(X, theta)
	grad_matrix <- gradient_matrix(X, Y, theta)
	gradient <- colMeans(grad_matrix) / n
	V <- cov(grad_matrix)
	
	fisher_information = V + (gradient %*% t(gradient))
	
	inv_hessian <- pseudo_inverse(h_matrix)
	
	ice_term <- sum(diag(fisher_information %*% inv_hessian)) / n
	
	# We demand that the correction term is positive, and not larger than the entropy term
	# itself. This basically reduces our vulnerability to ill-conditioned problems and numerical issues. 
	ice_term <- max(0.0, min(entropy, ice_term));
	
	output <- (entropy + ice_term)
	return (output)
}

theta_calc <- function(X1, Y1, lambda1, lambda2, start) {
	reg_val <- optim(start, reg_objective, gr = NULL, x = X1, y = Y1, l1 = lambda1, l2 = lambda2, control=list(maxit=1000))
	return (reg_val$par)
}

# Everything else is allowed to start from the MLE solution, so let the MLE calc
# also have two passes. 
theta_mle <- function(X1, Y1) {
	start <- rep(0.0, ncol(X1))
	t1 <- theta_calc(X1, Y1, 0.0, 0.0, start)
	t2 <- theta_calc(X1, Y1, 0.0, 0.0, t1)
	return (t2)
}

theta_en_solve <- function(X1, Y1, lambda1, lambda2) {
	# Elastic net theta solving. Allows us to do L1 or L2. 
	start <- rep(0.0, ncol(X1))
	
	# Allow the L2 and TIK approaches to start from the MLE solution.
	t1 <- theta_calc(X1, Y1, 0.0, 0.0, start)
	t2 <- theta_calc(X1, Y1, lambda1, lambda2, t1)
	return (t2)
}

theta_ice <- function(X, Y) {
	start <- rep(0.0, ncol(X))
	
	# Allow the ICE approach to start from the MLE solution.
	t1 <- theta_calc(X, Y, 0.0, 0.0, start)
	
	ice_val <- optim(t1, ice_objective, gr = NULL, X = X, Y = Y, control=list(maxit=1000))
	return (ice_val$par)
}

do_test_suite <- function(d, n, gamma, theta) {
	X <- generate_x(n, gamma)
	Y <- generate_y(X, theta, n)

	lambda_ridge <- lambda_opt_approx_ridge(X, Y)
	lambda_lasso <- lambda_opt_approx_lasso(X, Y)

	mle <- theta_mle(X, Y)
	ice <- theta_ice(X, Y)
	l2_aprx <- theta_en_solve(X, Y, 0.0, lambda_ridge)
	l1_aprx <- theta_en_solve(X, Y, lambda_lasso, 0.0)
	
	n_test <- (10000 + 2*n)
	
	X1 <- generate_x(n_test, gamma)
	Y1 <- generate_y(X1, theta, n_test)
	
	
	err_true <- sample_error(X1, Y1, theta)
	err_mle <- sample_error(X1, Y1, mle)
	err_l2 <- sample_error(X1, Y1, l2_aprx)
	err_l1 <- sample_error(X1, Y1, l1_aprx)
	err_ice <- sample_error(X1, Y1, ice)

	mse_mle_aprx <- err_mle - err_true
	mse_l2_aprx <- err_l2 - err_mle
	mse_l1_aprx <- err_l1 - err_mle
	mse_ice_aprx <- err_ice - err_mle

	print(paste("N*lambda_ridge: ", n * lambda_ridge, ", N*Lambda_lasso: ", n* lambda_lasso))

	output <- c(mse_mle_aprx, mse_l2_aprx, mse_l1_aprx, mse_ice_aprx)
	
	print(paste("Computed[", n, "][", d, "]: ", (n* lambda_ridge), " ", (n* lambda_lasso), " ", paste(output, collapse=" ")))
	
	return  (output)
}




space_study <- function(dims, missing, sizes, reps, c_num, seed) {
	all_results = c()
	
	set.seed(seed)
	
	for(w in 1:length(dims)) {
		d <- dims[w]
		ms <- missing[w]
	
		for(k in 1:length(sizes)) {
			n <- sizes[k]
			results <- c()
	
			print(paste("Starting[", n, "][", d, "]"))
			
			for(z in 1:reps) {
				
				if((z %% 10) == 1) {
					gamma <- generate_gamma(d, c_num)
					theta <- generate_theta(d, gamma)
	
					for(q in 2:(2 + ms)) {
						theta[q] <- 0.0				
					}
				}
				
				print(paste("Starting[", n, "][", d, "][", z, "]:", rcond(gamma)))
				
				next_val <- try(do_test_suite(d, n, gamma, theta))
				
				if(inherits(next_val, "try-error")) {
					print("Error in test suite, skipping.")
				}
				else if (!all(is.finite(next_val))) {
					print("Infinite or NaN values, skipping.")
				}
				else {
					results <- cbind(c(next_val), results)
				
					if(z > 2) {
						res_mu <- colMeans(t(results))
						res_sigma <- cov(t(results))
						
						sample_count <- ncol(results)
						res_var <- diag(res_sigma) / sample_count
					
						t_stat <- res_mu / sqrt(res_var)
						print(paste("T-Statistics[", n, "][", d, "][", ms, "][", z, "]: ", paste(t_stat, collapse=" ")))
					}					
				}
							
	
			}
			
			res_mu <- colMeans(t(results))
			res_sigma <- cov(t(results))		
			sample_count <- ncol(results)
			res_var <- diag(res_sigma) / sample_count
			t_stat <- res_mu / sqrt(res_var)
			
			print(paste("Done[", n, "][", d, "][", ms, "]"))
			print(paste("T-Statistics[", n, "][", d, "][", ms, "]: ", paste(t_stat, collapse=" ")))
			print("mle, l2, l1, ice")
			print(res_mu)
			print(res_sigma)
			print(cor(t(results)))
			
			combined_res = c(n, d, ms, res_mu, t_stat)
			all_results <- cbind(c(combined_res), all_results)
		}	
	}
	
	print("n, d, ms, mle_err, l2_err, l1_err, ice_err, mle_t, l2_t, ice_t")
	paste(t(all_results))

	return (t(all_results))
}

# Used to generate Table 4 in the paper
convergence_study <- function(d, ms, sizes, c_num, seed) {
	set.seed(seed)
	gamma <- generate_gamma(d, c_num)
	theta <- generate_theta(d, gamma)
	
	for(q in 2:(2 + ms)) {
		theta[q] <- 0.0				
	}
	
	
	X_all <- generate_x(max(sizes), gamma)
	Y_all <- generate_y(X_all, theta, max(sizes))
	
	n_test <- (10000 + 2*max(sizes))
		
	X1 <- generate_x(n_test, gamma)
	Y1 <- generate_y(X1, theta, n_test)
	
	results <- c()
	
	for(k in 1:length(sizes)) {
		n <- sizes[k]
	
		print(paste("Starting[", n, "][", d, "]"))
			
		X <- head(X_all, n)
		Y <- head(Y_all, n)
		
		lambda_ridge <- lambda_opt_approx_ridge(X, Y)
		lambda_lasso <- lambda_opt_approx_lasso(X, Y)
	
		mle <- theta_mle(X, Y)
		ice <- theta_ice(X, Y)
		l2_aprx <- theta_en_solve(X, Y, 0.0, lambda_ridge)
		l1_aprx <- theta_en_solve(X, Y, lambda_lasso, 0.0)
		
		err_true <- sample_error(X1, Y1, theta)
		err_mle <- sample_error(X1, Y1, mle)
		err_l2 <- sample_error(X1, Y1, l2_aprx)
		err_l1 <- sample_error(X1, Y1, l1_aprx)
		err_ice <- sample_error(X1, Y1, ice)
	
		output <- c(seed, d, ms, n, err_true, err_mle - err_true, err_l2 - err_true, err_l1 - err_true, err_ice - err_true)	
		
		print(paste("Done[", n, "][", d, "]: ", paste(output, collapse=" ")))
		results <- cbind(c(output), results)
		
		print(t(results))
	}	
	
	print("seed, d, ms, n, h_true, err_mle, err_l2, err_l1, err_ice")
	print(t(results))
}


# parameter space study.
dims <- c(5, 10, 20)
missing <- c(2, 4, 8)
sizes <- c(5e2, 1e3, 2e3, 5e3)
reps <- 300
c_num <- 0.0001

# used to generate Table 3 in the paper

space_study(dims, missing, sizes, reps, c_num, 123456)


# Convergence study
p <- 10 # dimensionality of input space
m <- 4 # model misspecification - missing variables
sizes <- c(5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5)
c_num <- 0.0001

convergence_study(p, m, sizes, c_num, 123456)

	

