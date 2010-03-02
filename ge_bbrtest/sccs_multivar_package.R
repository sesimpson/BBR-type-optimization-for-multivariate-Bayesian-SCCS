
args = commandArgs(trailing=TRUE)

# name of text file (output from perl) that is ready to be fed into R
datafile = args[1]


# multivariate form for sccs likelihood, takes sparse matrix X as input
#	where X represents a matrix of size (total numpds x num drugs)
	l = function(par, X) {
	
		X_par = X %*% par	
		exp_X_par = exp(X_par)
	
		eta = as.matrix(dat$EVT)
		offs = dat$OFFS

		offs_exp_X_par = offs * exp_X_par
	
		n_pid = as.matrix(apply(PIDlist, 1, function(pid) { sum(eta[dat$PID == pid]) }))
		denom_pid = as.matrix(apply(PIDlist, 1, function(pid) { sum(offs_exp_X_par[dat$PID == pid]) }))
	
		loglik = as.numeric((t(eta) %*% X_par) - (t(n_pid) %*% log(denom_pid)))

		return( loglik )
	}


# calculate gradient of log likelihood
#		par = vector of parameters
#		X = SPARSE matrix of covariates (sum_i(r_i) x p dimensional)
#			where r_i = # of risk periods for person i
# 
	grad <- function(par, X) {
		
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



# data from new format of flat tab-delimited text file, with header
#	PID	EVT	OFFS	D1	D2	D3	... etc.

	dat = read.table(datafile, fill=TRUE, row.names=NULL, 
				header=TRUE, sep="\t")        # will automatically pad with NAs
												        # for rows w/o max num of drugs
	
	nrows = dim(dat)[1]
	ncols = dim(dat)[2]
	drugcols = paste("D", 1:(ncols-3), sep="")

	dat$PID = as.numeric(dat$PID)


	# create covariate matrix, stored in sparse format (in "Matrix" package)
	
		# vector containing all (unique) drug id numbers
		druglist = unique(as.numeric(as.matrix(dat[drugcols])))
		druglist = druglist[!is.na(druglist)]
		druglist = druglist[druglist != 0]
		druglist = sort(druglist)
		ndrugs = length(druglist)

		# unique PIDS
		PIDlist = as.matrix(unique(dat$PID))

		# covariate matrix, in sparse format ("sparseMatrix" fxn in "Matrix" package)
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


# -------------------- RUN -----------------------------------------------------

#	cat("N = ", length(PIDlist), "\n", sep="") # N

	suppressPackageStartupMessages(
		require(maxLik, warn.conflicts=FALSE, quietly=TRUE)
	)
	
#	cat(date(), "\n")
	
	# with gradient, hessian (faster)
	ML = maxLik(function(beta) { l(beta, X.sM) }, 
				grad = function(beta) { grad(beta, X.sM) }, 
				hess = function(beta) { -fisherinfo(beta, X.sM) },
				start = numeric(ndrugs))

	#opt = optim(par=numeric(ndrugs), fn=function(beta) { l(beta,X.sM) },
	#			gr = function(beta) { grad(beta, X.sM) },
	#			hessian = TRUE)

#	print(summary(ML))
#	cat("\n")
	
#	cat(date(), "\n")

	beta_hat = ML$estimate

	info = fisherinfo(beta_hat, X.sM)
	SE_betahat = sqrt(diag(solve(info)))

	for (i in 1:ndrugs) {
		cat(druglist[i], "\t", beta_hat[i], "\t", SE_betahat[i], "\n", sep="")
	}

	# find 95% CIs (delta method)
#	z_ci = qnorm(0.975)
#	lo95CI = exp(beta_hat)-z_ci*exp(beta_hat)*SE_betahat	
#	hi95CI = exp(beta_hat)+z_ci*exp(beta_hat)*SE_betahat







