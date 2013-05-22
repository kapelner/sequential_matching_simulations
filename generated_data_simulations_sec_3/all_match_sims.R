NOT_GRID = length(grep("wharton.upenn.edu", Sys.getenv(c("HOSTNAME")))) == 0

if (NOT_GRID){
	setwd("c:/users/kapelner/workspace/StatTurk/matching_sims")
}

#tryCatch(library(optmatch), error = function(e){install.packages("optmatch")}, finally = library(optmatch))
#tryCatch(library(lme4), error = function(e){install.packages("lme4")}, finally = library(lme4))

###set simulation parameters here
args = commandArgs(TRUE)
if (length(args) > 0){
	for (i in 1 : length(args)){
		eval(parse(text = args[[i]]))
	}
}
if (NOT_GRID){
	iter_num = 6
}
print(paste("iter_num:", iter_num))

#how precise do we want our simulated results?
Nsim_per_block = 1000
Nsim_exact_test = 1000
Nsim_bootstrap = Nsim_exact_test

#the treatment effect is fixed so we can check unbiasedness, etc
beta_T = 1

#file name header
RESULTS_FILE_HEAD = paste(ifelse(beta_T == 0, "null_", ""), "results_quads_and_ints_betas_", sep = "")

#balanced assignment of treatment and control
prob_trt = 0.5 

#how should the X's and the errors be generated?
data_and_err_gen = "just_multnorm"

#what kind of models are we running?
sim_type = "series_of_quadratics_and_interactions"


#how many subjects enter the sequential experiment?
ns_to_test = c( 
	50,
	100,
	200,
	500
)

#cutoff parameter for matching
prob_match_cutoff_alphas = c(
	0.01,
	0.05,
	0.075,
	0.10,
	0.20,
	0.3,
	0.4,
	0.5,
	0.75,
	0.95
)

#how do we randomize subjects into treatment or control?
randomization_types = c(
	"van_bau",
	"efron_bau",
	"strat_bau",
	"ps_min_bau",
#	"pm_bau",
	"pm_kk",
#	"pm_kk_pbs",
	"van_linear",
	"efron_lin",
	"strat_lin",
	"ps_min_lin",
#	"pm_lin",
	"pm_lin_kk",
#	"pm_lin_kk_pbs",
	"van_exact",
	"efron_exact",	
	"strat_exact",
	"ps_min_exact",
#	"pm_exact",
	"pm_kk_exact"
)


#what do we measure for each set of simulations?
metrics_for_each_run = c(
	"avg_max_std_diff_bal", 
	"avg_max_ks_stat",
	"avg_beta_T", 
	"avg_abs_bias",
	"std_err_beta_T",
	"power",
	"pct_trt_diff",
	"res_end_prop_avg",
	"ssqr_over_sum",
	"ssqr_eq_ssqd_pval",
	"match_corr",
	"corr_eq_zero_pval",
	"pct_only_matches"
)

source("plot_response_models_balances_and_betas.R")
#run sims here
source(paste("sim_type_", sim_type, ".R", sep = ""))

