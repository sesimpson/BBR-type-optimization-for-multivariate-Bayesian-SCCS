# test BBR-like procedure for SCCS

#args = commandArgs(trailing=TRUE)

# name of text file (output from perl) that is ready to be fed into R
#datafile = args[1]

# avoid printing numbers in scientific notation
options(scipen=30)

# functions ------------------------------------------------------------------

	# multivariate form for sccs likelihood, takes sparse matrix X as input
	#	where X represents a matrix of size (total numpds x num drugs)
	l = function(par, X) {
	
		X_par = X %*% par	
		exp_X_par = exp(X_par)
		offs_exp_X_par = Matrix(offs * exp_X_par)
	
		denom_pid = t(pid_mat) %*% offs_exp_X_par

		#denom_pid = apply(PIDlist, 1, function(pid) { 
		#	sum(offs_exp_X_par[dat$PID == pid]) })
	
		loglik = as.numeric( (t(y) %*% X_par) - (t(n_pid) %*% log(denom_pid)) )

		return( loglik )
	}
	
	
	# calculate gradient of log likelihood
	#		par = vector of parameters
	#		X = SPARSE matrix of covariates (sum_i(r_i) x p dimensional)
	#			where r_i = num of risk periods for person i 
	grad = function(par, X) {
		
		eta = as.matrix(dat$EVT)
		offs = dat$OFFS
	
		X_par = X %*% par	
		exp_X_par = exp(X_par)
		offs_exp_X_par = Matrix(offs * exp_X_par)

		# denominator term for "pi" vector
		denom_pid = t(pid_mat) %*% offs_exp_X_par
			
		scale_pid = pid_mat %*% (n_pid/denom_pid)
			
#		scale_pid = apply(as.matrix(dat$PID), 1, 
#			function(pid) {
#				index = which(PIDlist == pid);
#				return( n_pid[index] / denom_pid[index] )
#			} )

		# holds the (n_i * pi_ik) terms
		scale_pi_pid = scale_pid * offs_exp_X_par
	
		gradient = as.numeric(t(y - scale_pi_pid) %*% X)
	
		return(gradient)
	}


	# fisher information where X is the sparse matrix of covariates (drugs), 
	# representing a matrix of dimension (total numpds x num drugs)
	fisherinfo <- function(par, X) {
		p = dim(X)[2]
	
		X_par = X %*% par	
		exp_X_par = exp(X_par)
		offs_exp_X_par = offs * exp_X_par

		# denominator term for "pi" vector
		denom_pid = t(pid_mat) %*% offs_exp_X_par
					
		x_offs_exp_X_par = X.sM * as.numeric(offs_exp_X_par)

		numer_pid = t(pid_mat) %*% x_offs_exp_X_par
				
		quot_pid = numer_pid / as.numeric(denom_pid)
		
		term1 = t(apply(quot_pid, 1, function(x) { x %o% x }))
				
		X_outer = t(apply(X, 1, function(x) { x %o% x }))
		
		term2 = t(apply(pid_mat, 2, function(pid) { 
			t(offs_exp_X_par[pid]) %*% X_outer[pid,]
		})) / as.numeric(denom_pid)		
#		# for first term
#		x_offs_exp_X_par = as.matrix(apply(pid_mat, 2, function(pid) {
#				return(
#					as.numeric(t(offs_exp_X_par[pid]) %*% X[pid,])
#				)
#			}))

#		# scale by denominator
#		x_offs_exp_X_par = t(x_offs_exp_X_par) / as.numeric(denom_pid) 
#		x_offs_exp_X_par_outer = t(apply(x_offs_exp_X_par, 1, function(x) { x %o% x }))

#		X_outer = t(apply(X, 1, function(x) { x %o% x }))
	
#		term2 = t(apply(pid_mat, 2, function(pid) { 
#			t(offs_exp_X_par[pid]) %*% X_outer[pid,]
#		})) / as.numeric(denom_pid)

#		total = matrix(t(n_pid) %*% (x_offs_exp_X_par_outer - term2), nrow=p, ncol=p)

		total = Matrix(as.numeric(t(n_pid) %*% (term1 - term2)), nrow=p, ncol=p)

		return(-total)
	}


	# BBR type update for parameter component beta_j
	update_beta_j <- function(j, beta, X, type) {
	
			X_beta = as.numeric(X %*% beta)
			exp_X_beta = exp(X_beta)
			offs_exp_X_beta = offs * exp_X_beta
			x_offs_exp_X_beta = X[,j] * offs_exp_X_beta
		
		# denominator term
		denom_pid = apply(PIDlist, 1, function(pid) { sum(offs_exp_X_beta[dat$PID == pid]) })
			
		# numerator term
		numer_pid = apply(PIDlist, 1, function(pid) { sum(x_offs_exp_X_beta[dat$PID == pid]) })
			
		# scale by denominator
		t1 = numer_pid / denom_pid
		
		if (type == "normal") {		# normal prior
			
			g_d1 = (n_pid %*% t1) - (X[,j] %*% eta) + (beta[j]/sigma2_beta)
			g_d2 = (n_pid %*% (t1 * (1-t1))) + (1/sigma2_beta)

			update = -as.numeric(g_d1/g_d2)

		} else {								# laplace prior
			
			g_d2 = n_pid %*% (t1 * (1-t1))
			g_d1_first = (n_pid %*% t1) - (X[,j] %*% eta)

			if ( sign(beta[j]) == 0 ) {
				neg_update = -as.numeric((g_d1_first - lambda)/g_d2)
				pos_update = -as.numeric((g_d1_first + lambda)/g_d2)
				
				# only one of these conditions will hold, from convexity
					update = ifelse( neg_update < 0, neg_update,
							ifelse( pos_update > 0, pos_update, 0) )
							
			} else {
				update = -as.numeric((g_d1_first + lambda*sign(beta[j]))/g_d2)
				
				# update would change sign, so return beta_j to 0
				if (sign(beta[j] + update) != sign(beta[j])) {
					update = -beta[j]
				}
			}
		}
		
		return( update )
	}
	
	
	
# RUN -------------------------------------------------------------------------
	# read in data
	datafile = "/Users/ses/Desktop/sccs_multivar_package/ge_bbrtest/cond_373474_ge_bbrtest_Rin.txt"
	datafile = "/Users/ses/Desktop/sccs_multivar_package/sm_ge_format_multivar_OUT_53drugs.txt"

	# data in format of flat tab-delimited text file, header: 
	#	PID	EVT	OFFS	D1	D2	D3	... etc.
	dat = read.table(datafile, fill=TRUE, row.names=NULL, 
				header=TRUE, sep="\t")        # will automatically pad with NAs
												        # for rows w/o max num of drugs
	nrows = dim(dat)[1]
	ncols = dim(dat)[2]
	drugcols = paste("D", 1:(ncols-3), sep="")

	dat$PID = as.numeric(dat$PID)
		
	# vector containing all (unique) drug id numbers
		druglist = unique(as.numeric(as.matrix(dat[drugcols])))
		druglist = druglist[!is.na(druglist)]
		druglist = druglist[druglist != 0]
		druglist = sort(druglist)
		ndrugs = length(druglist)

		# unique PIDS
		PIDlist = as.matrix(unique(dat$PID))

		# sparseMatrix: covariate matrix, in sparse format (in "Matrix" package)
			suppressPackageStartupMessages(
				require(Matrix, warn.conflicts=FALSE, quietly=TRUE)
			)
			
			indices = unlist(apply(dat[,drugcols], 2, function(col) {
				rbind(
					which(!is.na(col)),
					pmatch(col[!is.na(col)], druglist, duplicates.ok=TRUE)
					)}), use.names=FALSE)
			
			indices = array(indices, dim=c(2, length(indices)/2))
			
			indices = indices[,which(!is.na(indices[2,]))]
			
			X.sM = sparseMatrix(i=indices[1,], j=indices[2,], dims=c(nrows, ndrugs))

			X.sM = as(X.sM, "dgCMatrix")
	
	y = Matrix(dat$EVT)
	offs = dat$OFFS
		
		
	# blockwise matrix to keep track of PIDs
	pid_mat = sparseMatrix(i=1:nrows, j=pmatch(dat$PID, PIDlist, duplicates.ok=TRUE), 
				dims=c(nrows, length(PIDlist)))
				
	pid_mat = as(pid_mat, "dgCMatrix")
				
	# number of events for each person
	n_pid = t(pid_mat) %*% y
				

# loop to do BBR-type procedure ---------------------------------------------------

	# set prior parameters
		sigma2_beta = 100
		lambda = sqrt(2/20)
	
	# other initializations
		type = "normal"
		done = FALSE
		beta = numeric(ndrugs)
		delta = rep(2, ndrugs)
		iter = 1

		X_beta = as.numeric(X.sM %*% beta)
		e_Xbeta = exp(X_beta)
		offs_e_Xbeta = Matrix(offs * e_Xbeta)
		
system.time(
		while(!done) {
			beta_old = beta
		
			for (j in 1:ndrugs) {

#			for (j in sample(1:ndrugs)) {
			#j = sample(1:ndrugs, 1)
			
				# newton-raphson step for coordinate j
				x_j = Matrix(X.sM[,j])
				
				x_offs_e_Xbeta = x_j * offs_e_Xbeta

				numer_j = t(pid_mat) %*% x_offs_e_Xbeta
				denom_j = t(pid_mat) %*% offs_e_Xbeta
			
				quot_j = numer_j / denom_j
				
				grad_j = (t(n_pid) %*% quot_j) - (t(y) %*% x_j)
				hess_j = (t(n_pid) %*% (quot_j * (1-quot_j)))
				
				
				if (type == "normal") { # NORMAL
					
					grad_j = grad_j + (beta[j]/sigma2_beta)
					hess_j = hess_j + (1/sigma2_beta)
			
					del_beta = -as.numeric(grad_j / hess_j)
					
				} else {	 # LAPLACIAN	

					if (beta[j] == 0) {
						neg_update = -as.numeric((grad_j - lambda)/hess_j)
						pos_update = -as.numeric((grad_j + lambda)/hess_j)
					
						# only one of these conditions will hold, from convexity
						del_beta = ifelse( neg_update < 0, neg_update,
								ifelse( pos_update > 0, pos_update, 0) )
					} else {
						
						del_beta = -as.numeric((grad_j + lambda*sign(beta[j]))/hess_j)
				
						# update would change sign, so return beta_j to 0
						if (sign(beta[j] + del_beta) != sign(beta[j])) {
							del_beta = -beta[j]
						}
					}
				}
				
				# limit update to trust region
				del_beta = ifelse(del_beta < -delta[j], -delta[j], 
							ifelse(del_beta > delta[j], delta[j], del_beta))
			
				delta[j] = max(2*abs(del_beta), 0.5*delta[j])
				beta[j] = beta[j] + del_beta

				
				# update vector offs * exp( X * beta)
				del_beta_xj = del_beta * x_j
				X_beta = X_beta + del_beta_xj
				offs_e_Xbeta = offs_e_Xbeta * exp(del_beta_xj)
			}

			# update log posterior
			loglik = (t(y) %*% (log(offs) + X_beta)) - (t(n_pid) %*% denom_j)
				
			if (type=="normal") {
				logpost = loglik - 0.5*ndrugs*log(2*pi*sigma2_beta) - 0.5*(t(beta) %*% beta)/sigma2_beta
					
			} else {
				logpost = loglik + ndrugs*log(0.5*lambda) - lambda*sum(abs(beta))
			}
				
			logpost = as.numeric(logpost)
			
			#iter = iter + 1

			cat(beta, "\n")
			cat("logpost: ", logpost, "\n\n")
			
			#if (iter %% ndrugs == 0) {
				
#				if (type == "laplace") { # laplace
#					cat("log post:", l(beta, X.sM) + ndrugs*log(0.5*lambda) - lambda*sum(abs(beta)), "\n\n")
#				
#				} else {						# normal
#					cat("log post:", l(beta, X.sM) - 0.5*ndrugs*log(2*pi*sigma2_beta) 
#						- (0.5/sigma2_beta)*sum(beta^2), "\n\n")
#				}
		
				# calculate and check convergence criteria
					conv = as.numeric(t(abs(X.sM %*% (beta - beta_old))) %*% y)
					conv = conv / as.numeric(1 + (t(abs(X.sM %*% beta_old)) %*% y) )
		
					if (as.numeric(conv	) <= 0.00001) { done = TRUE }			#}
		}
)	# end system.time

	
# ---------------------- print out results ----------------------
	for (i in 1:ndrugs) {
		cat(druglist[i], "\t", beta[i], "\n", sep="")
	}


# ML procedures
	library(maxLik)
	ML = maxLik(logLik=function(par) { l(par, X.sM) }, grad=function(par) { grad(par, X.sM) },
		hess = function(par) { -fisherinfo(par, X.sM) }, start=numeric(ndrugs))
	

