

#install.packages("matrixStats")

library(matrixStats)
library(MASS)

## Consider a gaussian with mean = 1, variance = 0.2. 
## Fit (mu, sigma) using maximum likelihood, L2, and ICE. Also consider the version where we fit 1/sigma instad of sigma. 

buckets <- 100.0



calcYmean <- function(x, theta) {
	y_mean = theta[1]*sin(pi * x[,1]*x[,2]) + theta[2]*(x[,3] - theta[3])**2 + theta[4]*x[,4] + theta[5]*x[,5]
	return (y_mean)
}

probCalculationGaussian <- function(X, mu, sigma) {
	objective_var <- sigma * sigma
	centered <- X - mu
	exponent <- (centered * centered) / (2 * objective_var)
	rawProb <- exp(-exponent) / sqrt(2 * pi * objective_var)
	rawProb <- rawProb / buckets
	rawProb
}

probCalculation <- function(x, y, theta) {
	y_diff = calcYmean(x, theta) - y; # std normal if theta = theta_0
	sigma2 = theta[6] * theta[6];
	prob_ydiff = exp(-0.5*y_diff * y_diff / sigma2) / sqrt(2 * pi * sigma2)
	prob_ydiff = prob_ydiff / buckets;
	return (prob_ydiff);
}

computeEntropy <- function(x, y, theta) {
	entropy = -log(probCalculation(x, y, theta));
	return (entropy);
}

computeMuGradient <- function(x, y, theta) {
	partial_theta = matrix(c(sin(pi * x[,1]*x[,2]), (x[,3] - theta[3])**2, -2*theta[2]*(x[,3] - theta[3]), x[,4],x[,5], 0.0*x[,5]), nrow=nrow(x), ncol=6)
	return (partial_theta)
}

computeGradient <- function(x, y, theta) {
	partial_theta = computeMuGradient(x, y, theta);
	centered = (calcYmean(x, theta) - y);
	
	sigma = theta[6];
	
	grad = centered * partial_theta / (sigma * sigma);
	grad[,6] = -(centered * centered / (sigma * sigma * sigma)) + (1 / sigma);
	return (grad)
}

computeI <- function(x, y, theta) {
	grad = computeGradient(x, y, theta);
	
	avg_I = grad[1,] %*% t(grad[1,]);
	
	for(i in 2 : nrow(x)) {
		avg_I = avg_I + (grad[i,] %*% t(grad[i,]));
	}
	
	avg_I = avg_I / nrow(x);
	
	return (avg_I);
}

computeJ <- function(x, y, theta) {
	mu_gradient = computeMuGradient(x, y, theta);

	sigma = theta[6];

	base_J = 0.0 * (mu_gradient[1,] %*% t(mu_gradient[1,]));
	adjustment_A = 0.0;
	adjustment_B = 0.0;
	coeff_vector = (calcYmean(x, theta) - y);
		
	mu_adjustment = 0.0 * mu_gradient[1,];
	sigma_adjustment = 0.0;
		
	for(i in 1 : nrow(x)) {
		base_J = base_J + (mu_gradient[i,] %*% t(mu_gradient[i,]));
		coeff_i = coeff_vector[i];
		term_a = -2*(x[i,3] - theta[3]) * coeff_i;
		term_b = 2 * (x[i,2]) * coeff_i
		adjustment_A = adjustment_A + term_a;
		adjustment_B = adjustment_B + term_b;
		
		coeff2 = coeff_i / sigma;
		
		mu_adjustment = mu_adjustment + (-2 * coeff2 * mu_gradient[i,]);
		
		sigma_adjustment = sigma_adjustment + ((3 * coeff2 * coeff2) - 1);
	}
	
	for(i in 1 : 5) {
		base_J[6, i] = mu_adjustment[i];
		base_J[i, 6] = mu_adjustment[i];
	}
	
	base_J[6, 6] = sigma_adjustment;
	
	base_J = base_J / nrow(x);
	base_J = base_J / (sigma * sigma);
	adjustment_A = adjustment_A / nrow(x);
	adjustment_B = adjustment_B / nrow(x);

	# Add the adjustments
	base_J[3,2] = base_J[3,2] + adjustment_A;
	base_J[2,3] = base_J[2,3] + adjustment_A;
	base_J[3,3] = base_J[3,3] + adjustment_B;
	
	return (base_J);
}


objective_mle <- function(x, y, theta, l1, l2) {
	entropy <- computeEntropy(x, y, theta)	
	mean_entropy = mean(entropy)
	
	l1Norm <- sum(abs(theta))
	l2Norm <- sum(theta * theta)
	
	result <- mean_entropy + (l1 * l1Norm) + (l2 * l2Norm)
	
	return (result)
}

check_derivatives <- function(x, y, theta, scale) {
	check_theta = theta + scale * rnorm(6);
	
	grad_theta = colMeans(computeGradient(x, y, theta));
	grad_check = colMeans(computeGradient(x, y, check_theta));
	matrix_J = computeJ(x, y, theta);
	
	grad_delta = (check_theta - theta) %*% matrix_J;
	grad_delta2 = grad_check - grad_theta;
	
	theta_entropy = mean(computeEntropy(x, y, theta));
	check_entropy = mean(computeEntropy(x, y, check_theta));
	
	entropy_delta = check_entropy - theta_entropy;
	entropy_delta2 = sum((check_theta - theta) * grad_theta);
	
}

objective_ice <- function(x, y, theta) {
	entropy <- computeEntropy(x, y, theta)	
	mean_entropy = mean(entropy)
	
	matrix_I = computeI(x, y, theta);
	matrix_J = computeJ(x, y, theta);
	
	# for(i in 1 : nrow(matrix_I)) {
		# for(j in 1 : nrow(matrix_I)) {
			# if(i != j) {
				# matrix_J[i, j] = 0.0;
			# }
			
			# # if(i != 2 && i != 3) {
			# # # if(i> 3) {
				# # matrix_I[i, j] = (i == j);
				# # matrix_J[i, j] = (i == j);
			# # }
			# # # if(j > 3) {
			# # if(j != 2 && j != 3) {
				# # matrix_I[i, j] = (i == j);
				# # matrix_J[i, j] = (i == j);
			# # }
		# }
	# }
	
	
	# J_inverse = solve(matrix_J);
	J_inverse = ginv(matrix_J);
	matrix_IJ = matrix_I %*% J_inverse;
	
	# We demand that the correction term is positive, and not larger than the entropy term
	# itself. This basically reduces our vulnerability to ill-conditioned problems and numerical issues. 
	iceCorrection = max(0.0, sum(diag(matrix_IJ))) / nrow(x);
	iceCorrection = min(iceCorrection, mean_entropy);
	
	result <- mean_entropy + iceCorrection
	return (result)
}



solve_mle <- function(X, Y, startTheta, L1, L2) {
	if(L2 != 0 && L1 != 0) {
		startTheta = solve_mle(X, Y, startTheta, 0.0, 0.0)
	}
	
	mle_val <- optim(startTheta, objective_mle, gr = NULL, x = X, y=Y, l1 = L1, l2 = L2, control=list(maxit=10000))
	return (mle_val$par)
}

solve_ice <- function(X, Y, startTheta) {
	# We must always start from the MLE estimate.
	starting_point = solve_mle(X, Y, startTheta, 0.0, 0.0)

	mle_val <- optim(starting_point, objective_ice, gr = NULL, x = X, y=Y, control=list(maxit=10000))
	return (mle_val$par)
}

snapToGrid <- function(inputVal, bucketCount) {
	return (floor(inputVal * bucketCount) / bucketCount)
}

lambda_crossval_objective <- function(X, Y, start, lambda1, lambda2) {
	oos_total <- 0.0
	n_splits <- 2
	
	for(k in 1:n_splits) {
		test <- seq(k, nrow(X), n_splits)
		train <- seq(1, nrow(X), 1)
		train <- train[-test]
		
		X_train <- X[train,]
		X_test <- X[test,]
		Y_train <- Y[train];
		Y_test <- Y[test];
		
		nextTheta <- solve_mle(X_train, Y_train, start, lambda1, lambda2)
		oos_next <- objective_mle(X_test, Y_test, nextTheta, lambda1, lambda2)
		oos_total <- oos_total + oos_next
	}
	
	return (oos_total)
}

lambda_crossval_ridge <- function(X, Y, start, lambda) {
	return (lambda_crossval_objective(X, Y, start, lambda1 = 0.0, lambda))
}

lambda_opt_ridge <- function(X, Y, start) {
	lambda_0 <- optimize(lambda_crossval_ridge, X = X, Y=Y, start = start, c(0.0, 0.5))
	return (lambda_0$minimum)
}

compute_bucketset <- function(mu, sigma) {	
	lower <- mu - (10 * sigma);
	upper <- mu + (10 * sigma);
	intervalSize <- upper - lower
	bucketSet <- lower + (0:floor(intervalSize * buckets) / buckets)
	
	return (bucketSet)
}


rho <- function(X, theta0, theta) {
	trueY <- calcYmean(X, theta0);
	calcY <- calcYmean(X, theta);
	
	trueSigma2 = theta0[6] * theta0[6];
	sigma2 = theta[6] * theta[6];
	
	
	kl = 0.0;
	
	for(i in 1:nrow(X)) {
		centerPoint = trueY[i];
		calcCenter = calcY[i];
		centerDiff = abs(centerPoint - calcCenter);
		dev = centerDiff + 2.0;
		
		bucketSet <- compute_bucketset(centerPoint, dev)
		p <- probCalculationGaussian(bucketSet, centerPoint, sqrt(trueSigma2));
		q <- probCalculationGaussian(bucketSet, calcY[i], sqrt(sigma2));
		
		c1 <- bucketSet - centerPoint;
		c2 <- bucketSet - calcY[i];
		c1 <- (0.5 * c1 * c1 / trueSigma2) + 0.5 * log(trueSigma2);
		c2 <- (0.5 * c2 * c2 / sigma2) + 0.5 * log(sigma2);
		
		ratio2 <- c2 - c1;
		kl_term <- sum(p * ratio2);

		pSum = sum(p);
		qSum = sum(q);

		if(pSum < 0.9999 || qSum < 0.9999) {
			print(paste("Bucket sum: ", pSum, " : ", qSum, " (", kl_term, ")")); # Should be extremely close to 1.0
			print(paste("True Y: ", trueY[i]));
			print(paste("Calc Y: ", calcY[i]));
		}
		
		kl = kl + kl_term;
	}
	
	kl = kl / nrow(X);
	return (kl)
}




run_test <- function(trueTheta, startTheta, n, iteration) {
	#Generate test data.
	# LDDP this thing with 5 sigma cutoffs. 
	
	# X = matrix(rnorm(n * 5), nrow=n, ncol=5)
	X = matrix(runif(n * 5), nrow=n, ncol=5)
	Y = calcYmean(X, trueTheta) + rnorm(n);
	#X = snapToGrid(X, buckets)
	# Y = snapToGrid(Y, buckets)
	
	test_size = 100 + n;
	#XPrime = matrix(rnorm(n * 5), nrow=n, ncol=5)
	XPrime = matrix(runif(test_size * 5), nrow=test_size, ncol=5)
	# YPrime = calcYmean(XPrime, trueTheta) + rnorm(n);
	# XPrime = snapToGrid(XPrime, buckets)
	# YPrime = snapToGrid(YPrime, buckets)
	
	mle <- solve_mle(X, Y, startTheta, 0.0, 0.0)
	mle_entropy <- objective_mle(X, Y, mle, 0.0, 0.0)

	lambda2 = lambda_opt_ridge(X, Y, startTheta)
	l2_theta = solve_mle(X, Y, startTheta, 0.0, lambda2)
	l2_entropy <- objective_mle(X, Y, l2_theta, 0.0, 0.0)

	ice_theta = solve_ice(X, Y, startTheta)

	mle_kl = rho(XPrime, trueTheta, mle)
	l2_kl = rho(XPrime, trueTheta, l2_theta)
	ice_kl = rho(XPrime, trueTheta, ice_theta)

	#print(paste("MLE: ", mle_entropy, " -> ", mle_kl))
	#print(paste("L2[", lambda2, "]: ", l2_entropy, " -> ", l2_kl))
	
	print(paste(n, ", ", iteration, ", ", mle_kl, ", ", l2_kl, ", ", ice_kl, ", ", (l2_kl - mle_kl), ", ", (ice_kl - mle_kl)));
	
	
	return (c(mle_kl, l2_kl, ice_kl, (l2_kl - mle_kl), (ice_kl - mle_kl)))
}

# Conducts a K=2 friedman test given the vector of r_value's. (these are binary indicators of which "side" won each trial.)
friedman_pvalue <- function(r_value) {
	count = length(r_value);
	

	
	r_normalized = (sum(r_value) / count);
	r_n2 = 1.0 - r_normalized;
	
	# Ranks should be 1 and 2, not 0 and 1. 
	r_normalized = 1.0 + r_normalized;
	r_n2 = 1.0 + r_n2;
	
	# Now subtract k+1/2 = 3/2
	r_normalized = r_normalized -  3/2;
	r_n2 = r_n2 - 3/2;
	
	friedman_base = (r_normalized*r_normalized) + (r_n2 * r_n2);
	friedman_chi = 2.0 * count * friedman_base;
	friedman_p = pchisq(friedman_chi, df=2, lower.tail=FALSE);
	return (friedman_p);
}


run_test_repeated <- function(trueTheta, startTheta, n, count) {
	results <- c()

	print("n, iteration, mleKl, l2Kl, iceKl, l2Diff, iceDiff")
	
	# See https://en.wikipedia.org/wiki/Friedman_test
	# This is the r_j for the ICE estimator. 
	friedman_ice_vs_MLE = 0.0;
	friedman_ice_vs_L2 = 0.0;
	
	for(k in 1:count) {
		output <- run_test(trueTheta, startTheta, n, k);
		results <- cbind(c(output), results);
		
		if(results[3] < results[1]) {
			# ICE outperformed MLE
			friedman_ice_vs_MLE = friedman_ice_vs_MLE + 1.0;
		}
		if(results[3] < results[2]) {
			# ICE outperformed L2. 
			friedman_ice_vs_L2 = friedman_ice_vs_L2 + 1.0;
		}
	}
	
	repeated_results <- colMeans(t(results))
	repeated_devs <- sqrt(colVars(t(results)) / n)
	l2_diffs = (t(results[2,]) - t(results[1,]))[1,]
	ice_diffs = (t(results[3,]) - t(results[1,]))[1,]
	
	l2_dev = sqrt(var(l2_diffs) / count)
	ice_dev = sqrt(var(ice_diffs) / count)
	
	# T-statistics for the ICE and L2 estimators. 
	l2_z = mean(l2_diffs) / l2_dev;
	ice_z = mean(ice_diffs) / ice_dev;
	
	print(paste("MLE KL: ", repeated_results[1]))
	print(paste("L2 KL: ", repeated_results[2], " (worse by ", (repeated_results[2] - repeated_results[1]), ") [", l2_z, "]"))
	print(paste("ICE KL: ", repeated_results[3], " (worse (or better if negative) by ", (repeated_results[3] - repeated_results[1]), ") [", ice_z , "]"))
	
	# Calculate the p-value that ICE is better than MLE. 
	ice_vs_mle_p = friedman_pvalue(1.0 * (results[1,] > results[3,]));
	ice_vs_l2_p = friedman_pvalue(1.0 * (results[2,] > results[3,]));
	l2_vs_mle_p = friedman_pvalue(1.0 * (results[1,] > results[2,]));
	
	print(paste("ICE vs. MLE Friedman-P: ", ice_vs_mle_p));
	print(paste("ICE vs. L2 Friedman-P: ", ice_vs_l2_p));
	print(paste("L2 vs. MLE Friedman-P: ", l2_vs_mle_p));
	
	return (repeated_results)
}

true_theta = c(10.0, 20.0, 0.5, 10.0, 5.0, 1.0);
start_theta = true_theta + 0.1*rnorm(6);

# X = matrix(rnorm(10 * 5), nrow=10, ncol=5)
# Y = calcYmean(X, theta) + rnorm(ncol(X));

# X = matrix(runif(10 * 5), nrow=10, ncol=5)
# Y = calcYmean(X, true_theta) + rnorm(ncol(X));

r# un_test(true_theta, start_theta, 50, 0)

#arbitrary seed value. 
set.seed(0xcafebabe);

# This should almost always show that L2 is worse than MLE. 
repeated_results <- run_test_repeated(true_theta, start_theta, 8, 200);
repeated_results <- run_test_repeated(true_theta, start_theta, 16, 200);
repeated_results <- run_test_repeated(true_theta, start_theta, 32, 200);
repeated_results <- run_test_repeated(true_theta, start_theta, 64, 200);
repeated_results <- run_test_repeated(true_theta, start_theta, 128, 200);
repeated_results <- run_test_repeated(true_theta, start_theta, 256, 200);
repeated_results <- run_test_repeated(true_theta, start_theta, 512, 200);
repeated_results <- run_test_repeated(true_theta, start_theta, 1024, 200);



