source("common_pm.R")
source("common_pm_kk.R")

b_T_obs = mean(yTs) - mean(yCs)

## now we have to monte-carlo the exact test
b_T_sims = array(NA, Nsim_exact_test)
for (nsim_exact_test in 1 : Nsim_exact_test){
	#permute matched pairs
	permute_odds = rbinom(m, 1, prob_trt)
	permute_evens = 1 - permute_odds
	Xymatched$indic_T = as.vector(t(cbind(permute_odds, permute_evens)))
	
	#permute reservoir
	Xyleft$indic_T = sample(c(rep(1, nRT), rep(0, nRC)))
	
	#get ybarT - ybarC and record it
	yTs = c(Xymatched[Xymatched$indic_T == 1, ]$y, Xyleft[Xyleft$indic_T == 1, ]$y)
	yCs = c(Xymatched[Xymatched$indic_T == 0, ]$y, Xyleft[Xyleft$indic_T == 0, ]$y)
	b_T_sims[nsim_exact_test] = mean(yTs) - mean(yCs)
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
