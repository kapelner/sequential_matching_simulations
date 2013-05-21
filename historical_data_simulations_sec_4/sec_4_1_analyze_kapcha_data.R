setwd("c:/users/kapelner/workspace/StatTurk/data_analysis")
library(MASS)
X_orig = as.data.frame(read.table("dumps/kapcha_dump_coded_09_06_10__03_04_clean.txt", sep = "|", fill = TRUE, header = TRUE))
dim(X_orig)
table(X_orig$treatment_run)



#we only need control and treatments
#Xtc = X[X$treatment_run %in% c("control", "kapcha"), ]
#Xtc = X[X$treatment_run %in% c("kapcha"), ]
#
#dim(Xtc)
subjects = unique(X_orig$mturk_worker_id)
#kill black subjects
subjects = setdiff(subjects, "")

X_orig[X_orig$mturk_worker_id == "A1WFZTP2MKTHH7", ]

table(X_orig$qu_signature)
names(table(X_orig$qu_signature))
table(X_orig$num_windows_switches)
table(X_orig$treatment_run)




#covariates
covariates_from_matrix_embedded_in_qu_signature = c(
	"beer_question_thaler_1985",			#this is our response variable y	
	"dem_01_birth_year",
	"dem_02_gender",
	"dem_03_education",
#	"dem_04_why_work_on_mturk",				#too much of a pain to convert
	"dem_05_earnings",
	"dem_06_hours",
	"dem_07_percent_surveys",
	"dem_08_multitasking",
#		"feedback_question",				#can't use feedback because it really is after the experiment
#		"football_question_thaler_1985",	#we're only focusing on the beer question
	"imc_1_oppenheimer_2009",
	"motivation_gauge_oppenheimer_2009",
	"nfc_01_problems",
	"nfc_02_situation",
	"nfc_03_thinking_fun",
	"nfc_04_requires_thought",
	"nfc_05_think_in_depth"
#	"nfc_06_deliberation",					#the NFC question are pretty colinear so let's just use the first 5
#	"nfc_07_think_as_hard",					#the NFC question are pretty colinear so let's just use the first 5
#	"nfc_08_short_long_term",				#the NFC question are pretty colinear so let's just use the first 5
#	"nfc_09_require_little_thought",		#the NFC question are pretty colinear so let's just use the first 5
#	"nfc_10_relying_on_thought",			#the NFC question are pretty colinear so let's just use the first 5
#	"nfc_11_new_solutions",					#the NFC question are pretty colinear so let's just use the first 5
#	"nfc_12_new_ways",						#the NFC question are pretty colinear so let's just use the first 5
#	"nfc_13_puzzles",						#the NFC question are pretty colinear so let's just use the first 5
#	"nfc_14_thinking_abstractly",			#the NFC question are pretty colinear so let's just use the first 5
#	"nfc_15_intellectual",					#the NFC question are pretty colinear so let's just use the first 5
#	"nfc_16_mental_effort",					#the NFC question are pretty colinear so let's just use the first 5
#	"nfc_17_utilitarian",					#the NFC question are pretty colinear so let's just use the first 5
#	"nfc_18_deliberating_issues"			#the NFC question are pretty colinear so let's just use the first 5
)


#now design the new matrix with the columns where they should be.
other_covariates = c("num_windows_switches", "treatment_q1_place", "time_reading_directions", "treatment_run")
Xy = as.data.frame(matrix(NA, nrow = length(subjects), ncol = length(covariates_from_matrix_embedded_in_qu_signature) + length(other_covariates)))
colnames(Xy) = c(covariates_from_matrix_embedded_in_qu_signature, other_covariates)


#X[X$mturk_worker_id == subjects[i], c("qu_signature", "qu_responses")]
for (i in 1 : length(subjects)){
#	X[X$mturk_worker_id == subjects[i], ]
	for (cov in covariates_from_matrix_embedded_in_qu_signature){
		Xijraw = X_orig[X_orig$mturk_worker_id == subjects[i] & X_orig$qu_signature == cov, "qu_responses"]
		if (length(Xijraw) == 0){
			next
		}		
		if (cov == "imc_1_oppenheimer_2009"){ #the IMC pass has to be handled separately to convert to binary
			if (as.character(Xijraw) == "pass"){
#				cat(cov, "Xijraw:", Xijraw, "Xij:", 1, "\n")
				Xy[i, cov] = 1
			} else {
#				cat(cov, "Xijraw:", Xijraw, "Xij:", 0, "\n")
				Xy[i, cov] = 0
			}			
		} else {
			Xij = as.numeric(as.character(Xijraw))
#			cat(cov, "Xijraw:", Xijraw, "Xij:", Xij, "\n")

			Xy[i, cov] = Xij			
		}

	}
	Xy[i, "num_windows_switches"] = X_orig[X_orig$mturk_worker_id == subjects[i], "num_windows_switches"][1]
	trt = X_orig[X_orig$mturk_worker_id == subjects[i], "treatment_q1_place"][1]
	if (trt == "fancy"){
		Xy[i, "treatment_q1_place"] = 1
	} else if (trt == "run_down") {
		Xy[i, "treatment_q1_place"] = 0
	}	
	
	started_at_date = X_orig[X_orig$mturk_worker_id == subjects[i], "started_at_date"][1]
	started_at_time = X_orig[X_orig$mturk_worker_id == subjects[i], "started_at_time"][1] 
	read_directions_at_date = X_orig[X_orig$mturk_worker_id == subjects[i], "read_directions_at_date"][1]
	read_directions_at_time = X_orig[X_orig$mturk_worker_id == subjects[i], "read_directions_at_time"][1]
	
	t_started = as.POSIXct(strptime(paste(started_at_date, started_at_time), "%y/%m/%d %H:%M:%S"))
	t_read_directions = as.POSIXct(strptime(paste(read_directions_at_date, read_directions_at_time), "%y/%m/%d %H:%M:%S"))
	Xy[i, "time_reading_directions"] = as.numeric(t_read_directions - t_started)
	
	Xy[i, "treatment_run"] = as.character(X_orig[X_orig$mturk_worker_id == subjects[i], "treatment_run"][1])
}
head(Xy)

#kill missing data
Xy = na.omit(Xy)

#kill all where people paid more than $10 (2 rows)
Xy = Xy[Xy$beer_question_thaler_1985 <= 10, ]
#hist(Xy$beer_question_thaler_1985, br =50)
#sort(Xy$beer_question_thaler_1985)
#table(Xy$dem_08_multitasking)

#adjust birth year
#table(Xy$dem_01_birth_year)
Xy$age = 110 - Xy$dem_01_birth_year #Why 110? because 110 <=> 2010, the year the study was ran on MTurk
Xy$dem_01_birth_year = NULL

#table(Xy$treatment_run)

#make demographics terms factors
Xy$dem_02_gender = as.factor(Xy$dem_02_gender)
Xy$dem_03_education = as.factor(Xy$dem_03_education)
Xy$dem_05_earnings = as.factor(Xy$dem_05_earnings)
Xy$dem_06_hours = as.factor(Xy$dem_06_hours)
Xy$dem_07_percent_surveys = as.factor(Xy$dem_07_percent_surveys)
Xy$dem_08_multitasking = as.factor(Xy$dem_08_multitasking)

Xy$treatment_run = as.factor(Xy$treatment_run)

#Xy$dem_02_gender = as.numeric(Xy$dem_02_gender)
#Xy$dem_03_education = as.numeric(Xy$dem_03_education)
#Xy$dem_05_earnings = as.numeric(Xy$dem_05_earnings)
#Xy$dem_06_hours = as.numeric(Xy$dem_06_hours)
#Xy$dem_07_percent_surveys = as.numeric(Xy$dem_07_percent_surveys)
#Xy$dem_08_multitasking = as.numeric(Xy$dem_08_multitasking)


Xy_kapcha = Xy[Xy$treatment_run == "kapcha", ]
Xy_kapcha$treatment_run = NULL
Xy_exhortation = Xy[Xy$treatment_run == "exhortation", ]
Xy_exhortation$treatment_run = NULL
Xy_timing = Xy[Xy$treatment_run == "timing", ]
Xy_timing$treatment_run = NULL
Xy_control = Xy[Xy$treatment_run == "control", ]
Xy_control$treatment_run = NULL

#run some regressions
mod = lm(beer_question_thaler_1985 ~ ., Xy)
summary(mod)
mod = lm(beer_question_thaler_1985 ~ . - age + cut(age, breaks = seq(from = 10, to = 90, by = 10)), Xy)
summary(mod)
mod = lm(beer_question_thaler_1985 ~ 0 + ., Xy)
summary(mod)
mod = lm(beer_question_thaler_1985 ~ treatment_q1_place, Xy)
summary(mod)

mod = lm(beer_question_thaler_1985 ~ ., Xy_kapcha)
summary(mod)
mod = lm(beer_question_thaler_1985 ~ . - age + cut(age, breaks = seq(from = 10, to = 90, by = 10)), Xy_kapcha)
summary(mod)
mod = lm(beer_question_thaler_1985 ~ treatment_q1_place, Xy_kapcha)
summary(mod)

mod = lm(beer_question_thaler_1985 ~ ., Xy_exhortation)
summary(mod)
mod = lm(beer_question_thaler_1985 ~ . - age + cut(age, breaks = seq(from = 10, to = 90, by = 10)), Xy_exhortation)
summary(mod)
mod = lm(beer_question_thaler_1985 ~ treatment_q1_place, Xy_exhortation)
summary(mod)

mod = lm(beer_question_thaler_1985 ~ ., Xy_timing)
summary(mod)
mod = lm(beer_question_thaler_1985 ~ . - age + cut(age, breaks = seq(from = 10, to = 90, by = 10)), Xy_timing)
summary(mod)
mod = lm(beer_question_thaler_1985 ~ treatment_q1_place, Xy_timing)
summary(mod)

mod = lm(beer_question_thaler_1985 ~ ., Xy_control)
summary(mod)
mod = lm(beer_question_thaler_1985 ~ . - age + cut(age, breaks = seq(from = 10, to = 90, by = 10)), Xy_control)
summary(mod)
mod = lm(beer_question_thaler_1985 ~ treatment_q1_place, Xy_control)
summary(mod)
mod = lm(beer_question_thaler_1985 ~ treatment_q1_place + age + I(dem_05_earnings == 9) + I(dem_08_multitasking == 2) + nfc_01_problems, Xy_control)
summary(mod)

mod = lm(beer_question_thaler_1985 ~ age + motivation_gauge_oppenheimer_2009 + imc_1_oppenheimer_2009 + 
		as.numeric(dem_05_earnings) + as.numeric(dem_06_hours) + as.numeric(dem_02_gender) + 
		as.numeric(dem_08_multitasking) + 
		nfc_01_problems + nfc_05_think_in_depth + treatment_q1_place, Xy_control)

summary(mod)

#How much better can we do with BART and RF?
Xy_control2 = as.data.frame(model.matrix(lm(beer_question_thaler_1985 ~ age + motivation_gauge_oppenheimer_2009 + imc_1_oppenheimer_2009 + 
						as.numeric(dem_05_earnings) + as.numeric(dem_06_hours) + as.numeric(dem_02_gender) + 
						as.numeric(dem_08_multitasking) + 
						nfc_01_problems + nfc_05_think_in_depth + treatment_q1_place, Xy_control)))
summary(mod)
Xy_control2 = Xy_control2[, 2 : ncol(Xy_control2)]
Xy_control2$y = Xy_control$beer_question_thaler_1985
set_bart_machine_num_cores(4)
bart_machine = build_bart_machine(Xy = Xy_control2,
		num_trees = 200,
		num_burn_in = 1000,
		num_iterations_after_burn_in = 1000)
bart_machine
plot_y_vs_yhat(bart_machine)
plot_convergence_diagnostics(bart_machine)

library(randomForest)
rf = randomForest(x = Xy_control2[, 1 : (ncol(Xy_control2) - 1)], y = Xy_control2$y)
#find R^2
yhat = predict(rf, Xy_control2[, 1 : (ncol(Xy_control2) - 1)])
plot(Xy_control2$y, yhat)
abline(a = 0, b = 1)
sse = sum((Xy_control2$y - yhat)^2)
sst = sum((Xy_control2$y - mean(Xy_control2$y))^2)
Pseudo_Rsq_rf = 1 - sse / sst
Pseudo_Rsq_rf

Xy_control2 = as.data.frame(model.matrix(lm(beer_question_thaler_1985 ~ treatment_q1_place + age + I(dem_05_earnings == 9) + I(dem_08_multitasking == 2) + nfc_01_problems, Xy_control)))
Xy_control2 = Xy_control2[, 2 : ncol(Xy_control2)]
Xy_control2$y = Xy_control$beer_question_thaler_1985
set_bart_machine_num_cores(4)
bart_machine = build_bart_machine(Xy = Xy_control2,
		num_trees = 200,
		num_burn_in = 1000,
		num_iterations_after_burn_in = 1000)
bart_machine
plot_y_vs_yhat(bart_machine)
plot_convergence_diagnostics(bart_machine)

library(randomForest)
rf = randomForest(x = Xy_control2[, 1 : (ncol(Xy_control2) - 1)], y = Xy_control2$y)
#find R^2
yhat = predict(rf, Xy_control2[, 1 : (ncol(Xy_control2) - 1)])
plot(Xy_control2$y, yhat)
abline(a = 0, b = 1)
sse = sum((Xy_control2$y - yhat)^2)
sst = sum((Xy_control2$y - mean(Xy_control2$y))^2)
Pseudo_Rsq_rf = 1 - sse / sst
Pseudo_Rsq_rf

#How much better can we do with BART?
#Xy_control2 = Xy_control
#Xy_control2$y = Xy_control$beer_question_thaler_1985
#Xy_control2$beer_question_thaler_1985 = NULL
#Xy_control2$dem_02_gender = as.numeric(Xy_control2$dem_02_gender)
#Xy_control2$dem_03_education = as.numeric(Xy_control2$dem_03_education)
#Xy_control2$dem_05_earnings = as.numeric(Xy_control2$dem_05_earnings)
#Xy_control2$dem_06_hours = as.numeric(Xy_control2$dem_06_hours)
#Xy_control2$dem_07_percent_surveys = as.numeric(Xy_control2$dem_07_percent_surveys)
#Xy_control2$dem_08_multitasking = as.numeric(Xy_control2$dem_08_multitasking)
#head(Xy_control2)
#set_bart_machine_num_cores(4)
#bart_machine = build_bart_machine(Xy = Xy_control2,
#		num_trees = 200,
#		num_burn_in = 1000,
#		num_iterations_after_burn_in = 1000)
#bart_machine
#plot_y_vs_yhat(bart_machine)
#plot_convergence_diagnostics(bart_machine)
#
#library(randomForest)
#rf = randomForest(x = Xy_control2[, 1 : (ncol(Xy_control2) - 1)], y = Xy_control2$y)
##find R^2
#yhat = predict(rf, Xy_control2[, 1 : (ncol(Xy_control2) - 1)])
#plot(Xy_control2$y, yhat)
#abline(a = 0, b = 1)
#sse = sum((Xy_control2$y - yhat)^2)
#sst = sum((Xy_control2$y - mean(Xy_control2$y))^2)
#Pseudo_Rsq_rf = 1 - sse / sst
#Pseudo_Rsq_rf


#Xy_control2 = Xy_control
#Xy_control2$dem_03_education = NULL
#Xy_control2$dem_05_earnings_above_4 = ifelse(as.numeric(Xy_control$dem_05_earnings) > 4, 1, 0)
#Xy_control2$dem_05_earnings = NULL
#Xy_control2$dem_06_hours_above_4 = ifelse(as.numeric(Xy_control$dem_06_hours) > 4, 1, 0)
#Xy_control2$dem_06_hours = NULL
#Xy_control2$dem_07_percent_surveys = NULL
#Xy_control2$dem_08_multitasking_above_2 = ifelse(as.numeric(Xy_control$dem_08_multitasking) > 2, 1, 0)
#Xy_control2$dem_08_multitasking = NULL
#Xy_control2$nfc_02_situation = NULL    
#Xy_control2$nfc_03_thinking_fun = NULL    
#Xy_control2$nfc_04_requires_thought = NULL
#Xy_control2$num_windows_switches = NULL
#Xy_control2$age_less_than_30 = ifelse(as.numeric(Xy_control$age) < 30, 1, 0)
#Xy_control2$age_between_30_40 = ifelse(as.numeric(Xy_control$age) >= 30 & as.numeric(Xy_control$age) < 40, 1, 0)
#Xy_control2$age_between_40_50 = ifelse(as.numeric(Xy_control$age) >= 40 & as.numeric(Xy_control$age) < 50, 1, 0)
#Xy_control2$age = NULL
#head(Xy_control2)

#mod = lm(beer_question_thaler_1985 ~ ., Xy_control2)
#summary(mod)

#need sample var-cov matrix of design

#mod = lm(beer_question_thaler_1985 ~ . - dem_03_education - dem_07_percent_surveys - dem_08_multitasking - dem_06_hours - nfc_01_problems - nfc_02_situation  -nfc_04_requires_thought, Xy_control)
#summary(mod)
#mod = lm(beer_question_thaler_1985 ~ dem_05_earnings, Xy_control)
#table(Xy_control$dem_05_earnings)
#
#Xy_control_mod = Xy_control[Xy_control$dem_05_earnings %in% seq(1,5), ]
#table(Xy_control_mod$dem_05_earnings)
#
#Xy_control_mod$dem_05_earnings = as.numeric(Xy_control_mod$dem_05_earnings)
#Xy_control_mod$dem_05_earnings = as.factor(Xy_control_mod$dem_05_earnings)
#
#
#mod = lm(beer_question_thaler_1985 ~ dem_05_earnings, Xy_control_mod)
#table(Xy_control$dem_05_earnings)
#
#
#mm = model.matrix(beer_question_thaler_1985 ~ dem_05_earnings, Xy_control_mod)
#head(mm)
#mm = mm[, 2 : ncol(mm)]
#
#t(mm) %*% mm
#det(t(mm) %*% mm)
#var(mm)
#det(var(mm))
#
#
#
#names(mod)
#mm = model.matrix(beer_question_thaler_1985 ~ . - dem_03_education - dem_07_percent_surveys - dem_08_multitasking - dem_06_hours - nfc_01_problems - nfc_02_situation - nfc_04_requires_thought, Xy_control)
#mm = model.matrix(beer_question_thaler_1985 ~ ., Xy_control)
#mm = mm[, 2 : ncol(mm)]
#head(mm)
#det(var(mm))




Nsim = 200
n_study = 168 ### the pretend amount of people in this study
prob_match_cutoff_alpha = 0.10 ###found to be the best in the basic simulations

beta_matches = array(NA, Nsim)
beta_se_matches = array(NA, Nsim)
beta_regression = array(NA, Nsim)
beta_se_regression = array(NA, Nsim)
sub_ns = array(NA, Nsim)

Xy_controlT = Xy_control[Xy_control$treatment_q1_place == 1, ]
Xy_controlC = Xy_control[Xy_control$treatment_q1_place == 0, ]

for (nsim in 1 : Nsim){
	cat("nsim:", nsim," ")
	#scramble the order of the participants for each run
	X = Xy_control[sample(1 : nrow(Xy_control), n_study), ]
	#pull out true indicator of treatment and true y
	true_trt = X$treatment_q1_place
	y_true = X$beer_question_thaler_1985
	X$beer_question_thaler_1985 = NULL
	X$treatment_q1_place = NULL
	
	#make model matrix (turn factor cols into dummy cols)
#	Xmodel = as.data.frame(model.matrix(lm(y_true ~ 
#		age + motivation_gauge_oppenheimer_2009 + imc_1_oppenheimer_2009 + 
#		as.numeric(dem_05_earnings) + as.numeric(dem_06_hours) + as.numeric(dem_02_gender) + 
#		as.numeric(dem_08_multitasking) + 
#		nfc_01_problems + nfc_05_think_in_depth, X)))
	Xmodel = as.data.frame(model.matrix(lm(y_true ~ age + I(dem_05_earnings == 9) + I(dem_08_multitasking == 2) + nfc_01_problems, X)))


	
	#kill intercept
	Xmodel = Xmodel[, 2 : ncol(Xmodel)]

	indic_T = array(NA, n_study) 
	#initialize the reservoir
	match_indic = array(-1, n_study) #0 indicates reservoir
	
	for (i_match in 1 : n_study){
		
		xs_to_date = Xmodel[1 : i_match, ]
		p = ncol(xs_to_date)
		
		#if there is nothing in the reservoir, randomize and add it to reservoir
#		if (i_match > p && abs(det(var(xs_to_date))) < EFF_ZERO){
#			#we're going to add some noise to each column now
#			variances = diag(var(xs_to_date))
#			for (j in 1 : p){
#				xs_to_date[, j] = xs_to_date[, j] + rnorm(i_match, 0, ifelse(variances[j] == 0, 1, variances[j]) / 100)
#			}
#		}
		if (length(match_indic[match_indic == 0]) == 0 || i_match <= p + 1){ # || abs(det(var(xs_to_date))) < EFF_ZERO%
			indic_T[i_match] = true_trt[i_match]
			match_indic[i_match] = 0
			cat(".")
		} else {
			#first calculate the threshold we're operating at
			
			S_xs_inv = ginv(var(xs_to_date))
			F_crit =  qf(prob_match_cutoff_alpha, p, i_match - p)
			T_cutoff_sq = p * (n_study - 1) / (n_study - p) * F_crit
			#now iterate over all items in reservoir and take the minimum distance x
			reservoir_indices = which(match_indic == 0)
			x_star = Xmodel[i_match, ]
				
			sqd_distances = array(NA, length(reservoir_indices))
			for (r in 1 : length(reservoir_indices)){
				x_star_min_res = as.matrix(x_star - Xmodel[reservoir_indices[r], ])								
				sqd_distances[r] = 1 / 2 * x_star_min_res %*% S_xs_inv %*% t(x_star_min_res)	
			}
#			cat(paste("i_match", i_match, "sqd_distances", paste(sqd_distances, collapse = ", "), "T_cutoff_sq", T_cutoff_sq, "\n"))
			
			#find minimum distance index
			min_sqd_dist_index = which(sqd_distances == min(sqd_distances))
#			if (length(sqd_distances[min_sqd_dist_index]) > 1 || length(T_cutoff_sq) > 1){
#				cat(paste("i_match", i_match, "sqd_distances[min_sqd_dist_index]", paste(sqd_distances[min_sqd_dist_index], collapse = "!!!!!!!!!!!!!!!!"), "T_cutoff_sq", paste(T_cutoff_sq, collapse = "!!!!!!!!!!!!!!!!"), "\n"))
#				print(xs_to_date)
#			}
			#if it's smaller than the threshold, we're in business: match it
			#but ONLY if it matches the true randomization
#			cat("sqd_distances[min_sqd_dist_index]", sqd_distances[min_sqd_dist_index], "T_cutoff_sq", T_cutoff_sq, "\n")
			min_distance = sqd_distances[min_sqd_dist_index][1]
			min_distance_ind = reservoir_indices[min_sqd_dist_index][1]
			if (min_distance < T_cutoff_sq){
				if (1 - indic_T[min_distance_ind] == true_trt[i_match]){
					match_num = max(match_indic) + 1
					match_indic[min_distance_ind] = match_num
					match_indic[i_match] = match_num
					indic_T[i_match] = 1 - indic_T[min_distance_ind]
					cat("|")
#					cat("made match", match_num, "!\n")
				} else {
					cat("x")
#					cat("ditched a row!\n")
				}
			} else { #otherwise, randomize and add it to the reservoir ie keep the original randomization
				indic_T[i_match] = true_trt[i_match]
				match_indic[i_match] = 0
				cat(".")
			}
		}
	}
	
	indic_T_sub = indic_T[!is.na(indic_T)]
	y_true_sub = y_true[!is.na(indic_T)]
	match_indic_sub = match_indic[!is.na(indic_T)]
	Xy = Xmodel[!is.na(indic_T), ]
	sub_ns[nsim] = length(indic_T_sub)
	
	cat("  sub_n:", sub_ns[nsim])
	#now use our estimator
	#get reservoir data
	Xyleft = Xy[match_indic_sub == 0, ]
	YleftT = y_true_sub[indic_T_sub == 1]
	YleftC = y_true_sub[indic_T_sub == 0]
	YbarRTMinusYbarRC = mean(YleftT) - mean(YleftC)
	
	#get reservoir sample sizes
	nR = nrow(Xyleft)
	nRT = length(YleftT)
	nRC = length(YleftC)
	
	#get d_i's and its sample avg
	match_diffs = array(NA, max(match_indic_sub))
	matched_ts = array(NA, max(match_indic_sub))
	matched_cs = array(NA, max(match_indic_sub))
	
	for (match_id in 1 : max(match_indic_sub)){
		yT = y_true_sub[indic_T_sub == 1 & match_indic_sub == match_id]
		yC = y_true_sub[indic_T_sub == 0 & match_indic_sub == match_id]
		match_diffs[match_id] = yT - yC
		matched_ts[match_id] = yT
		matched_cs[match_id] = yC
	}
	Dbar = mean(match_diffs)
	cor_test = cor.test(matched_ts, matched_cs)
	cat("  r: ", cor_test$estimate, "pval:", cor_test$p.value, "\n")
	
	#compute reservoir sample variance
	ssqR = (var(YleftT) * (nRT - 1) + var(YleftC) * (nRC - 1)) / (nR - 2) * (1 / nRT + 1 / nRC)
	ssqM = var(match_diffs) / length(match_diffs)
	
	gamma_star_hat = ssqR / (ssqR + ssqM)
	
	b_T_est = gamma_star_hat * Dbar + (1 - gamma_star_hat) * YbarRTMinusYbarRC
	b_T_est_var = ssqR * ssqM / (ssqR + ssqM)
	b_T_est_sd = sqrt(b_T_est_var)
	
	beta_matches[nsim] = b_T_est
	beta_se_matches[nsim] = b_T_est_sd
	
	#we now do a two sample pooled t-test using lm on the same number of subjects keeping n_T and n_C fixed
	n_T = sum(indic_T_sub == 1)
	n_C = sum(indic_T_sub == 0)
	
	#do a linear regression for comparison purposes with same sub_n
	X = rbind(Xy_controlT[sample(1 : nrow(Xy_controlT), n_T), ], Xy_controlC[sample(1 : nrow(Xy_controlC), n_C), ])
#	Xmodel = as.data.frame(model.matrix(lm(beer_question_thaler_1985 ~ ., X)))
#	#kill intercept
#	Xmodel = Xmodel[, 2 : ncol(Xmodel)]
	
	linear_mod = lm(beer_question_thaler_1985 ~ treatment_q1_place, X)
	beta_regression[nsim] = coef(summary(linear_mod))[2, 1]
	beta_se_regression[nsim] = coef(summary(linear_mod))[2, 2]
}

#xs_to_date_match = cbind(xs_to_date, match_indic)
#xs_to_date_match[order(xs_to_date_match$match_indic), ]

mean(beta_se_matches^2, na.rm=T)
mean(beta_se_regression^2, na.rm=T)

#efficiency via matching - paired test due to sub_n being different run-run
mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
###approx sample size reduction due to method
1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1

t.test(beta_se_matches, beta_se_regression)
mean(sub_ns)

############## models for y_true ~
#		age + motivation_gauge_oppenheimer_2009 + imc_1_oppenheimer_2009 + 
#		as.numeric(dem_05_earnings) + as.numeric(dem_06_hours) + as.numeric(dem_02_gender) + 
#		as.numeric(dem_08_multitasking) + 
#		nfc_01_problems + nfc_05_think_in_depth


#> mean(beta_se_matches^2, na.rm=T)
#[1] 0.3006662
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.3367747
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.843083
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.4574307

#> mean(sub_ns)
#[1] 37.825

#> mean(beta_se_matches^2, na.rm=T)
#[1] 0.1682897
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.1687765
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.203385
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.169011
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = -0.4104, df = 385.951, p-value = 0.6818
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.01805057  0.01181653 
#sample estimates:
#		mean of x mean of y 
#0.4019146 0.4050316 
#
#> mean(sub_ns)
#[1] 72.035

#> #xs_to_date_match = cbind(xs_to_date, match_indic)
#		> #xs_to_date_match[order(xs_to_date_match$match_indic), ]
#		> 
#		> mean(beta_se_matches^2, na.rm=T)
#[1] 0.1058364
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.1019394
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.056688
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.0536468
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = 1.2322, df = 370.754, p-value = 0.2186
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.002863977  0.012478008 
#sample estimates:
#		mean of x mean of y 
#0.3223534 0.3175464 
#
#> mean(sub_ns)
#[1] 116.055






############## models for y_true ~ age + I(dem_05_earnings == 9) + I(dem_08_multitasking == 2) + nfc_01_problems


#[1] 0.2456825
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.333376
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 2.005802
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.5014463
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = -5.0213, df = 388.702, p-value = 7.834e-07
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.11401085 -0.04985087 
#sample estimates:
#		mean of x mean of y 
#0.4708249 0.5527558 
#
#> mean(sub_ns)
#[1] 34.875



#> mean(beta_se_matches^2, na.rm=T)
#[1] 0.1356702
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.1723224
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.596562
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.3736541
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = -6.0675, df = 388.359, p-value = 3.088e-09
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.06473161 -0.03304753 
#sample estimates:
#		mean of x mean of y 
#0.3590183 0.4079079 
#
#> mean(sub_ns)
#[1] 67.77



#> 
#		> mean(beta_se_matches^2, na.rm=T)
#[1] 0.07605734
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.1115186
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.570489
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.3632555
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = -17.4246, df = 393.965, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.0655511 -0.0522587 
#sample estimates:
#		mean of x mean of y 
#0.2735053 0.3324102 
#
#> mean(sub_ns)
#[1] 112.14


#nsim: 1  .............|.xx....|....xxx.||xx.|.xx||x||xx.x|.  sub_n: 37  r:  -0.2383489 pval: 0.507236 
#nsim: 2  ................xxx..x.xx||.xxx.|.|.||..|x.xxxxx.|  sub_n: 35  r:  0.5262348 pval: 0.1803348 
#nsim: 3  ................x.x..||.xxx|.x.xxxxxxx|.x.|x.||x.|  sub_n: 34  r:  -0.4245226 pval: 0.2944836 
#nsim: 4  ...........|...|...|x.|.|.x|.|x.|....xx|||x.||x.x|  sub_n: 42  r:  0.1485558 pval: 0.6122606 
#nsim: 5  ................x.x||.|.|..xx|xx.|.x.x.||.|.x|xx..  sub_n: 39  r:  0.303559 pval: 0.3938454 
#nsim: 6  ................x..xx.x..xx.||.x|x|xxx.||x|.|.xxx|  sub_n: 35  r:  0.2873965 pval: 0.4533366 
#nsim: 7  .............x.|..||..|xx.x..|.|x|.xx|x.|.x|xxx|x|  sub_n: 37  r:  -0.3785392 pval: 0.2249887 
#nsim: 8  ............xx.|x|xx|..x.xxx.x..|x||.|.xx.x..x|x.|  sub_n: 34  r:  0.356562 pval: 0.3462404 
#nsim: 9  ..............|...xx||.x|.|||.....|.||x|x.|xx..|.x  sub_n: 42  r:  -0.09723812 pval: 0.7519872 
