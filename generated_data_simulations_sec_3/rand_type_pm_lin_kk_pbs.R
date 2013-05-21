source("common_pm.R")
source("common_pm_kk.R")

#now we want to control for delta x and do parametric bootstrap

##First need the observed estimate

XyD = as.data.frame(matrix(NA, nrow = m, ncol = ncol(Xymatched)))
colnames(XyD) = colnames(Xymatched)

for (im in 1 : m){
	Xyim = Xymatched[Xymatched$match_indic == im, ]
	XyimT = Xyim[Xyim$indic_T == 1, ]
	XyimC = Xyim[Xyim$indic_T == 0, ]
	XyD[im, ] = as.numeric(XyimT - XyimC)
}

linear_mod_matched = lm(y ~ . - match_indic - indic_T, data = XyD)
coefs_matched = coef(summary(linear_mod_matched))

if (nRT <= 2 || nRC <= 2){
	b_T_obs = coefs_matched[1, 1]
	pct_only_matchess[nsim] = 1
} else {
	linear_mod_reservoir = lm(y ~ . - match_indic, data = Xyleft)
	
	coefs_reservoir = coef(summary(linear_mod_reservoir))
	beta_match_regression = coefs_reservoir[p + 2, 1]
	ssqd_reservoir_regression = coefs_reservoir[p + 2, 2]^2 #lin mod returns SE not VAR, so square it
	ssqd_match_regression = coefs_matched[1, 2]^2
	gamma_star_hat_samp = ssqd_reservoir_regression / (ssqd_reservoir_regression + ssqd_match_regression) #just a convenience for faster runtime	
	b_T_obs = gamma_star_hat_samp * beta_match_regression + (1 - gamma_star_hat_samp) * beta_match_regression #proper weighting	
	pct_only_matchess[nsim] = 0	
}

####get stats for parametric bootstrap and treat them as parameters
sD_plain = sqrt(ssqD_bar * m)
sR_plain = sqrt(ssqR / (1 / nRT + 1 / nRC))

## now we have to simulate from the parametric setup
b_T_sims = array(NA, Nsim_bootstrap)
for (nsim_bootstrap in 1 : Nsim_bootstrap){
	
	#sample from matched pairs
	ydiffs_sample = rnorm(m, 0, sD_plain)
	XyD_samp = XyD
	XyD_samp$y = ydiffs_sample	
	linear_mod_matched_samp = lm(y ~ . - match_indic - indic_T, data = XyD_samp)
	coefs_matched_samp = coef(summary(linear_mod_matched_samp))
	beta_match_regression_samp = coefs_matched_samp[1, 1]
	ssqd_match_regression_samp = coefs_matched_samp[1, 2]^2
	
	if (nRT <= 2 || nRC <= 2){
		b_T_sims[nsim_bootstrap] = beta_match_regression_samp
		
	} else {
		#sample from reservoir
		YleftT_samp = rnorm(nRT, 0, sR_plain)
		YleftC_samp = rnorm(nRC, 0, sR_plain)
		
		Xyleft_samp = Xyleft
		Xyleft_samp$y[1 : nRC] = YleftC_samp
		Xyleft_samp$y[(nRC + 1) : (nRC + nRT)] = YleftT_samp
		
		#compute estimator reservoir sample std error
		linear_mod_reservoir_samp = lm(y ~ . - match_indic, data = Xyleft_samp)		
		coefs_reservoir_samp = coef(summary(linear_mod_reservoir_samp))
		beta_reservoir_regression_samp = coefs_reservoir_samp[p + 2, 1]
		ssqd_reservoir_regression_samp = coefs_reservoir_samp[p + 2, 2]^2 #lin mod returns SE not VAR, so square it
		
		gamma_star_hat_samp = ssqd_reservoir_regression_samp / (ssqd_reservoir_regression_samp + ssqd_match_regression_samp)	
		
		b_T_sims[nsim_bootstrap] = gamma_star_hat_samp * beta_match_regression_samp + (1 - gamma_star_hat_samp) * beta_reservoir_regression_samp
	}
}

#hist(b_T_sims, br = 100)
#mean(b_T_sims)
#sum(b_T_sims_unc > 1) / Nsim_exact_test
#sum(b_T_sims_cond > 1) / Nsim_exact_test
#hist(b_T_sims_unc, br = 100)
#hist(b_T_sims_cond, br = 100)
#ks.test(b_T_sims_unc, b_T_sims_cond)

#this is the empirical two-sided p-value based on simulation
pval = sum(abs(b_T_obs) < abs(b_T_sims)) / Nsim_bootstrap
