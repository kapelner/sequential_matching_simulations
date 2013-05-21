max_std_diff_balances = array(0, Nsim_per_block)
max_ks_stats = array(0, Nsim_per_block)
beta_Ts = array(NA, Nsim_per_block)
T_stats = array(NA, Nsim_per_block)
final_reservoir_size = array(NA, Nsim_per_block)
Rsqs = array(NA, Nsim_per_block)
Ha_acceptances = array(NA, Nsim_per_block)
pct_Ts = array(NA, Nsim_per_block)
ssqr_over_sums = array(NA, Nsim_per_block)
ssqr_eq_ssqdbar_pvals = array(NA, Nsim_per_block)
matched_correlations = array(NA, Nsim_per_block)
corr_eq_zero_pvals = array(NA, Nsim_per_block)
pct_only_matchess = array(NA, Nsim_per_block)

for (nsim in 1 : Nsim_per_block){
	if (nsim %% 50 == 0){
		cat(".")
	}
	
	#generate data
	source(paste("data_and_err_gen_", data_and_err_gen, ".R", sep = ""))
	
	#run one run of whatever simulation type
	source(paste("rand_type_", randomization_type, ".R", sep = ""))
	Ha_acceptances[nsim] = ifelse(pval < 0.05, 1, 0)

	#calculate balance metrics
	xTs = Xy[Xy$indic_T == 1, 1 : p]
	xCs = Xy[Xy$indic_T == 0, 1 : p]
	
	max_std_diff_balance = -999999
	for (j in 1 : p){
		std_diff_balance = abs(mean(xTs[, j]) - mean(xCs[, j])) / sqrt(var(xTs[, j]) / length(xTs[, j]) + var(xCs[, j]) / length(xCs[, j]))
		if (std_diff_balance > max_std_diff_balance){
			max_std_diff_balances[nsim] = std_diff_balance
		}
	}
	max_ks_stat = -99999
	for (j in 1 : p){
		ks_stat = ks.test(xTs[, j], xCs[, j])$statistic
		if (ks_stat > max_ks_stat){
			max_ks_stats[nsim] = ks_stat
		}
	}	
	
	#cross match test???? Don't really care since we don't relaly care about balance to begin with...
	
	
	#plot response models
	if (nsim < NUM_PLOTS_PER_RESP){
		tryCatch({plot_response_model(nsim, beta_Ts[nsim], final_reservoir_size[nsim])})
	}
	pct_Ts[nsim] = nrow(xTs) / n
}
#tryCatch({plot_balance_and_betas(randomization_type, n)}, error = function(){}) 

cat("\n")
#add to results
results["avg_max_std_diff_bal", randomization_type] = mean(max_std_diff_balances, na.rm = TRUE)
results["avg_beta_T", randomization_type] = mean(beta_Ts, na.rm = TRUE)	
results["avg_abs_bias", randomization_type] = mean(abs(beta_Ts - beta_T), na.rm = TRUE)
results["avg_max_ks_stat", randomization_type] = mean(max_ks_stats, na.rm = TRUE)
results["std_err_beta_T", randomization_type] = sd(beta_Ts, na.rm = TRUE)
results["power", randomization_type] = mean(Ha_acceptances, na.rm = TRUE)
results["pct_trt_diff", randomization_type] = mean(abs(pct_Ts - prob_trt), na.rm = TRUE)
#results["Rsq_avg", randomization_type] = mean(Rsqs, na.rm = TRUE)

nR_PMF = table(final_reservoir_size * n)
write.csv(nR_PMF, file = paste("nR_PMF_n_", n, "_lambda_", prob_match_cutoff_alpha, ".csv", sep  = ""), row.names = FALSE)
pdf(file = paste("nR_PMF_n_", n, "_lambda_", prob_match_cutoff_alpha, ".pdf", sep  = ""))
barplot(nR_PMF, main = paste("nR_PMF_n_", n, "_lambda_", prob_match_cutoff_alpha))
dev.off()
windows()
barplot(nR_PMF, main = paste("nR_PMF_n_", n, "_lambda_", prob_match_cutoff_alpha))


#only save reservoir data for matching algorithm
if (length(grep("pm", randomization_type)) > 0){
	results["res_end_prop_avg", randomization_type] = mean(final_reservoir_size)
}

if (randomization_type == "pm_kk" || randomization_type == "pm_linear_kk"){
	results["ssqr_over_sum", randomization_type] = mean(ssqr_over_sums, na.rm = TRUE)
	results["ssqr_eq_ssqd_pval", randomization_type] = mean(ssqr_eq_ssqdbar_pvals, na.rm = TRUE)
	results["pct_only_matches", randomization_type] = mean(pct_only_matchess, na.rm = TRUE)
}
if (randomization_type == "pm_kk"){
	results["match_corr", randomization_type] = mean(matched_correlations, na.rm = TRUE)
	results["corr_eq_zero_pval", randomization_type] = mean(corr_eq_zero_pvals, na.rm = TRUE)
}