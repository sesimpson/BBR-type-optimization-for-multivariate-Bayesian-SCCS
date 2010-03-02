# test BBR-like procedure for SCCS

args = commandArgs(trailing=TRUE)

# name of text file (output from perl) that is ready to be fed into R
datafile = args[1]

# avoid printing numbers in scientific notation
options(scipen=30)

# functions ------------------------------------------------------------------

	# multivariate form for sccs likelihood, takes sparse matrix X as input
	#	where X represents a matrix of size (total numpds x num drugs)
	l = function(par, X) {
	
		X_par = X %*% par	
		exp_X_par = exp(X_par)
		offs_exp_X_par = offs * exp_X_par
	
		denom_pid = apply(PIDlist, 1, function(pid) { 
			sum(offs_exp_X_par[dat$PID == pid]) })
	
		loglik = as.numeric( (eta %*% X_par) - (n_pid %*% log(denom_pid)) )

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
		offs_exp_X_par = offs * exp_X_par

		# number of events for each person
		n_pid = as.matrix(apply(PIDlist, 1, function(pid) { sum(eta[dat$PID == pid]) }))
	
		# denominator term for "pi" vector
		denom_pid = as.matrix(apply(PIDlist, 1, function(pid) { sum(offs_exp_X_par[dat$PID == pid]) }))

		scale_pid = apply(as.matrix(dat$PID), 1, 
			function(pid) {
				index = which(PIDlist == pid);
				return( n_pid[index] / denom_pid[index] )
			} )

		# holds the (n_i * pi_ik) terms
		scale_pi_pid = scale_pid * offs_exp_X_par
	
		gradient = as.numeric(t(eta - scale_pi_pid) %*% X)
	
		return(gradient)
	}


	# fisher information where X is the sparse matrix of covariates (drugs), 
	# representing a matrix of dimension (total numpds x num drugs)
	fisherinfo <- function(par, X) {
		p = dim(X)[2]
		eta = as.matrix(dat$EVT)
		offs = dat$OFFS
	
		X_par = X %*% par	
		exp_X_par = exp(X_par)
		offs_exp_X_par = offs * exp_X_par

		# number of events for each person
		n_pid = as.matrix(apply(PIDlist, 1, function(pid) { sum(eta[dat$PID == pid]) }))
	
		# denominator term for "pi" vector
		denom_pid = as.matrix(apply(PIDlist, 1, function(pid) { sum(offs_exp_X_par[dat$PID == pid]) }))
	
		# for first term
		x_offs_exp_X_par = as.matrix(apply(PIDlist, 1, function(pid) {
				indices = which(dat$PID == pid)
				return(
					as.numeric(t(offs_exp_X_par[indices]) %*% X[indices,])
				)
			}))

		# scale by denominator
		x_offs_exp_X_par = t(x_offs_exp_X_par) / as.numeric(denom_pid) 
		x_offs_exp_X_par_outer = t(apply(x_offs_exp_X_par, 1, function(x) { x %o% x }))

		X_outer = t(apply(X, 1, function(x) { x %o% x }))
	
		term2 = t(apply(PIDlist, 1, function(pid) { 
			t(offs_exp_X_par[dat$PID == pid]) %*% X_outer[dat$PID == pid,]
		})) / as.numeric(denom_pid)

		total = matrix(t(n_pid) %*% (x_offs_exp_X_par_outer - term2), nrow=p, ncol=p)

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

			i = NULL
			j = NULL
			for (d in 1:ndrugs) {
				this_row = which(dat[drugcols] == druglist[d], arr.ind=TRUE)[,1]
				i = c(i, this_row)
				j = c(j, rep(d, length(this_row)))
			}
			X.sM = sparseMatrix(i=i, j=j, dims=c(nrows, ndrugs))

			X.sM = as(X.sM, "dgCMatrix")
	
	eta = dat$EVT
	offs = dat$OFFS
		
	# number of events for each person
	n_pid = apply(PIDlist, 1, function(pid) { sum(eta[dat$PID == pid]) })



# loop to do BBR-type procedure ---------------------------------------------------

	# set prior parameters
		sigma2_beta = 1000
		lambda = sqrt(2/100)
	
	# other initializations
		type = "laplace"
		done = FALSE
		beta = numeric(ndrugs)
		delta = rep(5, ndrugs)
		iter = 1

#	system.time(
		while(!done) {
			beta_old = beta
		
			#for (j in sample(1:ndrugs, 1)) {
			j = sample(1:ndrugs, 1)
				update = update_beta_j(j, beta, X.sM, type)
			
				update = ifelse(update < -delta[j], -delta[j], 
							ifelse(update > delta[j], delta[j], update))
			
				delta[j] = max(2*abs(update), 0.5*delta[j])
			
				beta[j] = beta[j] + update
			#}
			
			iter = iter + 1
		
			if (iter %% ndrugs == 0) {
#				cat(beta, "\n")
		
#				if (type == "laplace") { # laplace
#					cat("log post:", l(beta, X.sM) + ndrugs*log(0.5*lambda) - lambda*sum(abs(beta)), "\n\n")
#				
#				} else {						# normal
#					cat("log post:", l(beta, X.sM) - 0.5*ndrugs*log(2*pi*sigma2_beta) 
#						- (0.5/sigma2_beta)*sum(beta^2), "\n\n")
#				}
		
				# calculate and check convergence criteria
					conv = t(abs(X.sM %*% (beta - beta_old))) %*% eta
					conv = conv / (1 + (t(abs(X.sM %*% beta_old)) %*% eta) )
		
					if (as.numeric(conv	) <= 0.00001) { done = TRUE }			}
		}
#	)	# end system.time

	
# ---------------------- print out results ----------------------
	for (i in 1:ndrugs) {
		cat(druglist[i], "\t", beta[i], "\n", sep="")
	}

	
	

