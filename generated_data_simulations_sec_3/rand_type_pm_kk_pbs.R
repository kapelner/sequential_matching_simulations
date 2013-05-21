source("common_pm.R")
source("common_pm_kk.R")

b_T_obs = gamma_star_hat * d_bar + (1 - gamma_star_hat) * YbarRTMinusYbarRC

####get stats for parametric bootstrap and treat them as parameters
sD_plain = sqrt(ssqD_bar * m)
sR_plain = sqrt(ssqR / (1 / nRT + 1 / nRC))

## now we have to simulate from the parametric setup
b_T_sims = array(NA, Nsim_bootstrap)
for (nsim_bootstrap in 1 : Nsim_bootstrap){
	
	#sample from matched pairs
	ydiffs_sample = rnorm(m, 0, sD_plain)	
	d_bar_samp = mean(ydiffs_sample)
	ssqD_bar_samp = var(ydiffs_sample) / m
	
	if (nRC <= 1 || nRT <= 1){
		b_T_sims[nsim_bootstrap] = d_bar_samp
	} else {
		#sample from reservoir
		YleftT_samp = rnorm(nRT, 0, sR_plain)
		YleftC_samp = rnorm(nRC, 0, sR_plain)
		YbarRTMinusYbarRC_samp = mean(YleftT_samp) - mean(YleftC_samp)
		ssqR_samp = (var(YleftT_samp) * (nRT - 1) + var(YleftC_samp) * (nRC - 1)) / (nR - 2) * (1 / nRT + 1 / nRC)	
		
		#now do the estimation
		gamma_star_hat_samp = ssqR_samp / (ssqR_samp + ssqD_bar_samp)
		b_T_sims[nsim_bootstrap] = gamma_star_hat_samp * d_bar_samp + (1 - gamma_star_hat_samp) * YbarRTMinusYbarRC_samp
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

beta_Ts[nsim] = NaN
