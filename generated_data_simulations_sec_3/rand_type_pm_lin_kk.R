source("common_pm.R")
source("common_pm_kk.R")

###########
#now we want to control for delta x


#sometimes the reservoir just isn't large enough
if (nRT <= 2 || nRC <= 2){
	XyD = matrix(NA, nrow = m, ncol = ncol(Xymatched))
	colnames(XyD) = colnames(Xymatched)
	
	for (im in 1 : m){
		Xyim = Xymatched[Xymatched$match_indic == im, ]
		XyimT = Xyim[Xyim$indic_T == 1, ]
		XyimC = Xyim[Xyim$indic_T == 0, ]
		XyD[im, ] = as.numeric(XyimT - XyimC)
	}
	linear_mod_matched = lm(y ~ . - match_indic - indic_T, data = as.data.frame(XyD))
	coefs = coef(summary(linear_mod_matched))
	
	beta_Ts[nsim] = coefs[1, 1]
	T_stats[nsim] = coefs[1, 3]
	pval = coefs[1, 4]
	pct_only_matchess[nsim] = 1
	
} else {
	#compute estimator from matched pairs by regression
	XyD = matrix(NA, nrow = m, ncol = ncol(Xymatched))
	colnames(XyD) = colnames(Xymatched)
	
	for (im in 1 : m){
		Xyim = Xymatched[Xymatched$match_indic == im, ]
		XyimT = Xyim[Xyim$indic_T == 1, ]
		XyimC = Xyim[Xyim$indic_T == 0, ]
		XyD[im, ] = as.numeric(XyimT - XyimC)
	}
	
	linear_mod_matched = lm(y ~ . - match_indic - indic_T, data = as.data.frame(XyD))
	coefs_matched = coef(summary(linear_mod_matched))
	beta_match_regression = coefs_matched[1, 1]	
	ssqd_match_regression = coefs_matched[1, 2]^2 #lin mod returns SE not VAR, so square it
	
	
	#compute estimator reservoir sample std error
	linear_mod_reservoir = lm(y ~ . - match_indic, data = Xyleft)
	
	coefs_reservoir = coef(summary(linear_mod_reservoir))
	beta_reservoir_regression = coefs_reservoir[p + 2, 1]
	ssqd_reservoir_regression = coefs_reservoir[p + 2, 2]^2 #lin mod returns SE not VAR, so square it
	
	gamma_star_hat = ssqd_reservoir_regression / (ssqd_reservoir_regression + ssqd_match_regression) #just a convenience for faster runtime	
	ssqr_over_sums[nsim] = gamma_star_hat	
	
	beta_Ts[nsim] = gamma_star_hat * beta_match_regression + (1 - gamma_star_hat) * beta_reservoir_regression #proper weighting	
	
	
	b_T_est_sd = sqrt(ssqd_match_regression * ssqd_reservoir_regression / (ssqd_match_regression + ssqd_reservoir_regression)) #analagous eq's
		
	T_stats[nsim] = beta_Ts[nsim] / b_T_est_sd #we calculate the T-stat against the null of zero effect	
	pval = 2 * (1 - pnorm(abs(T_stats[nsim]))) #approximate by using real Z
	
	ssqr_over_sums[nsim] = ssqd_reservoir_regression / (ssqd_reservoir_regression + ssqd_match_regression)	
	ssqr_eq_ssqdbar_pvals[nsim] = pf(ssqd_reservoir_regression * 2 / ssqd_match_regression, nRT + nRC - 2, m - 1, lower.tail = F)
	pct_only_matchess[nsim] = 0
}

