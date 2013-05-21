source("common_pm.R")
source("common_pm_kk.R")

b_T_obs = gamma_star_hat * d_bar + (1 - gamma_star_hat) * YbarRTMinusYbarRC

## now we have to monte-carlo the exact test
b_T_sims = array(NA, Nsim_exact_test)
for (nsim_exact_test in 1 : Nsim_exact_test){
	
	#### unconditional permuting - flip coins for each (equivalent to conditional permuting)
	trt_vec_multiple = (rbinom(max(match_indic), 1, prob_trt) - 0.5) * 2
	ydiffs_copy = ydiffs * trt_vec_multiple
	d_bar_samp = mean(ydiffs_copy)
	ssqD_bar_samp = var(ydiffs_copy) / length(ydiffs_copy)
	
	### conditional permuting for the reservoir
	Xyleft$indic_T = sample(c(rep(1, nRT), rep(0, nRC)))	
	YleftT_samp = Xyleft[Xyleft$indic_T == 1, ]$y
	YleftC_samp = Xyleft[Xyleft$indic_T == 0, ]$y
	YbarRTMinusYbarRC_samp = mean(YleftT_samp) - mean(YleftC_samp)
	ssqR_samp = (var(YleftT_samp) * (nRT - 1) + var(YleftC_samp) * (nRC - 1)) / (nR - 2) * (1 / nRT + 1 / nRC)	
	
	#now compute b_T_sim from the permuted stats
	gamma_star_hat_samp = ssqR_samp / (ssqR_samp + ssqD_bar_samp)
	
	b_T_sims[nsim_exact_test] = gamma_star_hat_samp * d_bar_samp + (1 - gamma_star_hat_samp) * YbarRTMinusYbarRC_samp
	
}

#hist(b_T_sims, br = 100)
#mean(b_T_sims)
#sum(b_T_sims_unc > 1) / Nsim_exact_test
#sum(b_T_sims_cond > 1) / Nsim_exact_test
#hist(b_T_sims_unc, br = 100)
#hist(b_T_sims_cond, br = 100)
#ks.test(b_T_sims_unc, b_T_sims_cond)

#this is the empirical two-sided p-value based on simulation
pval = sum(abs(b_T_obs) < abs(b_T_sims)) / Nsim_exact_test

beta_Ts[nsim] = NaN
