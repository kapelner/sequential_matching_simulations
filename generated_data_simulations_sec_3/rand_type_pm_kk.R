source("common_pm.R")
source("common_pm_kk.R")


#sometimes the reservoir just isn't large enough
if (nRT <= 1 || nRC <= 1){
	beta_Ts[nsim] = d_bar
	b_T_est_sd = sqrt(ssqD_bar)
	pct_only_matchess[nsim] = 1
	
} else {
	ssqR = (var(YleftT) * (nRT - 1) + var(YleftC) * (nRC - 1)) / (nR - 2) * (1 / nRT + 1 / nRC) #pooled variance in reservoir	
	gamma_star_hat = ssqR / (ssqR + ssqD_bar) #just a convenience for faster runtime	
	
	
	beta_Ts[nsim] = gamma_star_hat * d_bar + (1 - gamma_star_hat) * YbarRTMinusYbarRC #proper weighting	
	b_T_est_sd = sqrt(ssqR * ssqD_bar / (ssqR + ssqD_bar)) #see eq's

	#diagnostic stuf
	ssq_from_ssqr = ssqR / (1 / nRT + 1 / nRC)
	ssq_from_ssq_Dbar = ssqD_bar * m
	ssqr_over_sums[nsim] = ssq_from_ssqr / (ssq_from_ssqr + ssq_from_ssq_Dbar)
	ssqr_eq_ssqdbar_pvals[nsim] = pf(ssq_from_ssqr * 2 / ssq_from_ssq_Dbar, nRT + nRC - 2, m - 1, lower.tail = F)
	
	yMT = Xymatched[Xymatched$indic_T == 1, ]$y
	yCT = Xymatched[Xymatched$indic_T == 0, ]$y
	matched_correlations[nsim] = cor(yMT, yCT)
#	corr_eq_zero_pvals[nsim] = cor.test(yMT, yCT, alternative = "greater", method = "pearson")$p.value
	pct_only_matchess[nsim] = 0
}

T_stats[nsim] = beta_Ts[nsim] / b_T_est_sd #we calculate the T-stat against the null of zero effect	
pval = 2 * (1 - pnorm(abs(T_stats[nsim]))) #approximate by using real Z