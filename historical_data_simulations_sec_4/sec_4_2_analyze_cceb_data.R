setwd("c:/users/kapelner/workspace/StatTurk/data_analysis")

library(sas7bdat)

Xraw = read.sas7bdat("dumps/final_ic01_02oct08.sas7bdat")
dim(Xraw)
head(Xraw)
colnames(Xraw)

#[1] "subjid"                       "apf2"                         "apf9"                        
#[4] "date_collected2"              "date_collected9"              "egf2"                        
#[7] "egf9"                         "hb_egf2"                      "hb_egf9"                     
#[10] "initials2"                    "initials9"                    "perc_inhib2"                 
#[13] "perc_inhib9"                  "se_egf2"                      "se_egf9"                     
#[16] "se_hb_egf2"                   "se_hb_egf9"                   "se_inhib2"                   
#[19] "se_inhib9"                    "small2"                       "small9"                      
#[22] "Utr2a2"                       "Utr2a9"                       "Utr7a2"                      
#[25] "Utr7a9"                       "Utr11a2"                      "Utr11a9"                     
#[28] "Utr122"                       "Utr129"                       "ccid"                        
#[31] "sk"                           "responder"                    "SUBJ"                        
#[34] "ARM"                          "sex"                          "age"                         
#[37] "racial_type"                  "b_pain1"                      "b_urgency1"                  
#[40] "b_freq1"                      "b_freq_cat1"                  "b_pain2"                     
#[43] "b_urgency2"                   "b_freq2"                      "b_freq_cat2"                 
#[46] "b_pain"                       "b_urgency"                    "b_frequency"                 
#[49] "b_freq_cat"                   "withdrew"                     "compgrp"                     
#[52] "group"                        "comp_percent"                 "gra_orig"                    
#[55] "max6wk"                       "max12wk"                      "cciddose"                    
#[58] "dose"                         "GRA_bin_exclmiss"             "GRA_scale8"                  
#[61] "GRA_bin_scale8"               "GRA_carryforward"             "GRA_bin_carryforward"        
#[64] "SST1A"                        "SST1B"                        "SST1C"                       
#[67] "SST1D"                        "SST1E"                        "SST1F"                       
#[70] "SST1G"                        "SST1H"                        "SST1I"                       
#[73] "b_HAD15"                      "wk6_HAD15"                    "wk12_HAD15"                  
#[76] "chg_base6_HAD15"              "chg_base12_HAD15"             "b_IIEF1"                     
#[79] "wk6_IIEF1"                    "wk12_IIEF1"                   "chg_base6_IIEF1"             
#[82] "chg_base12_IIEF1"             "b_IIEF2"                      "wk6_IIEF2"                   
#[85] "wk12_IIEF2"                   "chg_base6_IIEF2"              "chg_base12_IIEF2"            
#[88] "b_WISIC"                      "wk6_WISIC"                    "wk12_WISIC"                  
#[91] "chg_base6_WISIC"              "chg_base12_WISIC"             "b_WISIC_IC7"                 
#[94] "wk6_WISIC_IC7"                "wk12_WISIC_IC7"               "chg_base6_WISIC_IC7"         
#[97] "chg_base12_WISIC_IC7"         "b_PCS"                        "wk6_PCS"                     
#[100] "wk12_PCS"                     "chg_base6_PCS"                "chg_base12_PCS"              
#[103] "b_MCS"                        "wk6_MCS"                      "wk12_MCS"                    
#[106] "chg_base6_MCS"                "chg_base12_MCS"               "b_FSFI_total"                
#[109] "wk6_FSFI_total"               "wk12_FSFI_total"              "chg_base6_FSFI_total"        
#[112] "chg_base12_FSFI_total"        "b_FSFI_desire"                "wk6_FSFI_desire"             
#[115] "wk12_FSFI_desire"             "chg_base6_FSFI_desire"        "chg_base12_FSFI_desire"      
#[118] "b_FSFI_arousal"               "wk6_FSFI_arousal"             "wk12_FSFI_arousal"           
#[121] "chg_base6_FSFI_arousal"       "chg_base12_FSFI_arousal"      "b_FSFI_lubrication"          
#[124] "wk6_FSFI_lubrication"         "wk12_FSFI_lubrication"        "chg_base6_FSFI_lubrication"  
#[127] "chg_base12_FSFI_lubrication"  "b_FSFI_orgasm"                "wk6_FSFI_orgasm"             
#[130] "wk12_FSFI_orgasm"             "chg_base6_FSFI_orgasm"        "chg_base12_FSFI_orgasm"      
#[133] "b_FSFI_satisfaction"          "wk6_FSFI_satisfaction"        "wk12_FSFI_satisfaction"      
#[136] "chg_base6_FSFI_satisfaction"  "chg_base12_FSFI_satisfaction" "b_FSFI_pain"                 
#[139] "wk6_FSFI_pain"                "wk12_FSFI_pain"               "chg_base6_FSFI_pain"         
#[142] "chg_base12_FSFI_pain"         "b_SIQTOT"                     "wk6_SIQTOT"                  
#[145] "wk12_SIQTOT"                  "chg_base6_SIQTOT"             "chg_base12_SIQTOT"           
#[148] "b_PIQTOT"                     "wk6_PIQTOT"                   "wk12_PIQTOT"                 
#[151] "chg_base6_PIQTOT"             "chg_base12_PIQTOT"            "b_VoidFreq"                  
#[154] "wk6_VoidFreq"                 "wk12_VoidFreq"                "chg_base6_VoidFreq"          
#[157] "chg_base12_VoidFreq"          "b_VoidFreqVLGSleep"           "wk6_VoidFreqVLGSleep"        
#[160] "wk12_VoidFreqVLGSleep"        "chg_base6_VoidFreqVLGSleep"   "chg_base12_VoidFreqVLGSleep" 
#[163] "wk6_pain"                     "wk12_pain"                    "chg_base6_pain"              
#[166] "chg_base12_pain"              "wk6_urgency"                  "wk12_urgency"                
#[169] "chg_base6_urgency"            "chg_base12_urgency"           "wk6_freq"                    
#[172] "wk12_freq"                    "chg_base6_freq"               "chg_base12_freq"             
#[175] "ct"                           "b_pain_cat"                   "b_urge_cat"                  
#[178] "b_frequency_cat"              "age_group"                    "white"                       
#[181] "missing_bm"                   "racenum"                      "armbin"                      
#[184] "b_pain_cat_collapsed"         "b_urge_cat_collapsed"         "b_freq_cat_collapsed"        
#[187] "hisp"                         "wk1_stc4a"                    "wk1_stc4b"                   
#[190] "wk1_stc4c"                    "wk1_stc4d"                    "wk1_stc5a"                   
#[193] "wk1_stc5b"                    "wk1_stc5c"                    "wk1_stc5d"                   
#[196] "wk2_stc4a"                    "wk2_stc4b"                    "wk2_stc4c"                   
#[199] "wk2_stc4d"                    "wk2_stc5a"                    "wk2_stc5b"                   
#[202] "wk2_stc5c"                    "wk2_stc5d"                    "wk3_stc4a"                   
#[205] "wk3_stc4b"                    "wk3_stc4c"                    "wk3_stc4d"                   
#[208] "wk3_stc5a"                    "wk3_stc5b"                    "wk3_stc5c"                   
#[211] "wk3_stc5d"                    "wk4_stc4a"                    "wk4_stc4b"                   
#[214] "wk4_stc4c"                    "wk4_stc4d"                    "wk4_stc5a"                   
#[217] "wk4_stc5b"                    "wk4_stc5c"                    "wk4_stc5d"                   
#[220] "wk5_stc4a"                    "wk5_stc4b"                    "wk5_stc4c"                   
#[223] "wk5_stc4d"                    "wk5_stc5a"                    "wk5_stc5b"                   
#[226] "wk5_stc5c"                    "wk5_stc5d"                    "wk6_stc4a"                   
#[229] "wk6_stc4b"                    "wk6_stc4c"                    "wk6_stc4d"                   
#[232] "wk6_stc5a"                    "wk6_stc5b"                    "wk6_stc5c"                   
#[235] "wk6_stc5d"                    "wk12_stc4a"                   "wk12_stc4b"                  
#[238] "wk12_stc4c"                   "wk12_stc4d"                   "wk12_stc5a"                  
#[241] "wk12_stc5b"                   "wk12_stc5c"                   "wk12_stc5d"                  
#[244] "dose_bin"                     "wk1_gra"                      "wk2_gra"                     
#[247] "wk3_gra"                      "wk4_gra"                      "wk5_gra"                     
#[250] "wk6_gra"                      "wk12_gra"                     "wk1_STC1"                    
#[253] "wk2_STC1"                     "wk3_STC1"                     "wk4_STC1"                    
#[256] "wk5_STC1"                     "wk6_STC1"                     "wk12_STC1"                   
#[259] "wk1_STC2"                     "wk2_STC2"                     "wk3_STC2"                    
#[262] "wk4_STC2"                     "wk5_STC2"                     "wk6_STC2"                    
#[265] "wk12_STC2"                    "wk1_STC2A"                    "wk2_STC2A"                   
#[268] "wk3_STC2A"                    "wk4_STC2A"                    "wk5_STC2A"                   
#[271] "wk6_STC2A"                    "wk12_STC2A"                   "wk1_STC2B"                   
#[274] "wk2_STC2B"                    "wk3_STC2B"                    "wk4_STC2B"                   
#[277] "wk5_STC2B"                    "wk6_STC2B"                    "wk12_STC2B"                  
#[280] "wk1_STC2D"                    "wk2_STC2D"                    "wk3_STC2D"                   
#[283] "wk4_STC2D"                    "wk5_STC2D"                    "wk6_STC2D"                   
#[286] "wk12_STC2D"                   "wk1_responder"                "wk2_responder"               
#[289] "wk3_responder"                "wk4_responder"                "wk5_responder"               
#[292] "wk6_responder"                "wk12_responder"               "lastvisit_vnum"              
#[295] "fuptime_cat"                  "fuptime_cont"                 "RC_ID"                       
#[298] "MED3"                         "uti"                          "pid"                         
#[301] "end"                          "vul"                          "MED8"                        
#[304] "MED9"                         "prost"                        "MED11"                       
#[307] "MED12"                        "asthma"                       "drugaller"                   
#[310] "foodaller"                    "skinaller"                    "sin"                         
#[313] "rhi"                          "latex"                        "MED12H"                      
#[316] "MED13"                        "bowel"                        "MED13B"                      
#[319] "MED13C"                       "MED13D"                       "MED14"                       
#[322] "diab"                         "MED14B"                       "MED14C"                      
#[325] "MED14D"                       "MED15"                        "fatigue"                     
#[328] "MED15B"                       "MED15C"                       "MED15D"                      
#[331] "MED16"                        "MED16A"                       "MED16B"                      
#[334] "MED16C"                       "MED16D"                       "MED16E"                      
#[337] "sextrans"                     "MED17A"                       "MED17B"                      
#[340] "MED17C"                       "MED17D"                       "MED17E"                      
#[343] "MED17F"                       "MED17G"                       "MED17H"                      
#[346] "MED17I"                       "MED18"                        "MED18A"                      
#[349] "MED18B"                       "MED18C"                       "MED18D"                      
#[352] "MED18E"                       "MED19"                        "lum"                         
#[355] "MED19B"                       "MED19C"                       "mig"                         
#[358] "MED19E"                       "MED19F"                       "MED20"                       
#[361] "fib"                          "auto"                         "MED20C"                      
#[364] "MED21"                        "MED21A"                       "hyst"                        
#[367] "MED21C"                       "surg"                         "MED22A"                      
#[370] "MED22ACK"                     "MED22B"                       "MED22BCK"                    
#[373] "MED23"                        "everdiag"                     "age1diag"                    
#[376] "eversym"                      "age1sym"                      "income"                      
#[379] "education"                    "employment"                   "partner"                     
#[382] "dob"                          "agecat"                       "fatifue"                     
#[385] "med1adif"                     "med2adif"  

hist(Xraw$age)
table(Xraw$age, useNA = "always")
table(Xraw$sex, useNA = "always")
table(Xraw$sextrans, useNA = "always")
#table(Xraw$sin, useNA = "always")
#table(Xraw$skinaller, useNA = "always")
table(Xraw$uti, useNA = "always")
table(Xraw$white, useNA = "always")
table(Xraw$hisp, useNA = "always")
#hist(Xraw$age1sym)
#table(Xraw$age1sym, useNA = "always")
table(Xraw$education, useNA = "always")
#table(Xraw$surg, useNA = "always")
#table(Xraw$asthma, useNA = "always")
table(Xraw$employment, useNA = "always")
table(Xraw$b_frequency, useNA = "always")
table(Xraw$b_urgency1, useNA = "always")
#table(Xraw$diab, useNA = "always")
#table(Xraw$fatigue, useNA = "always") #kill 88's
table(Xraw$partner, useNA = "always")
#table(Xraw$pid, useNA = "always")
#table(Xraw$fib, useNA = "always")
#table(Xraw$rhi, useNA = "always")
#table(Xraw$foodaller, useNA = "always")

#table(Xraw$hyst, useNA = "always")
#table(Xraw$income, useNA = "always")
#table(Xraw$lum, useNA = "always")
#table(Xraw$mig, useNA = "always")
#######TRT
table(Xraw$arm, useNA = "always")
table(Xraw$b_SIQTOT, useNA = "always")
table(Xraw$b_HAD15, useNA = "always")
#table(Xraw$b_IIEF1, useNA = "always")
table(Xraw$b_MCS, useNA = "always")
table(Xraw$b_PCS, useNA = "always")
table(Xraw$b_PIQTOT, useNA = "always")
#table(Xraw$b_VoidFreq, useNA = "always")
#table(Xraw$b_VoidFreqVLGSleep, useNA = "always")
table(Xraw$b_WISIC, useNA = "always")
table(Xraw$b_freq1, useNA = "always")
table(Xraw$b_pain, useNA = "always")
table(Xraw$b_pain1, useNA = "always")
table(Xraw$b_urgency, useNA = "always")

######RESP
hist(Xraw$chg_base12_pain, useNA = "always")

#make the analysis frame

Xy = Xraw[, c("age", "sex", "sextrans", "uti", "white", 
				"hisp", "education", "employment", 
				"b_SIQTOT", "b_HAD15", "b_MCS", "b_PCS", "b_PIQTOT", "b_WISIC",
				"b_frequency", "b_pain", "b_urgency", "partner", "armbin", 
				"chg_base12_pain")]
Xy$sex = Xy$sex - 1
Xy$sextrans = Xy$sextrans - 1
Xy$uti = Xy$uti - 1
Xy$education = Xy$education - 1
Xy$employment = Xy$employment - 1


Xy$y = Xy$chg_base12_pain
Xy$chg_base12_pain = NULL



Xy = na.omit(Xy)
dim(Xy)
mod = lm(y ~ ., Xy)
summary(mod)

#duplicates table 2 row 1 page 1856 of Foster et al (2010)
t.test(Xy[Xy$armbin == 0, ]$y, Xy[Xy$armbin == 1, ]$y)
#t = -1.1268, df = 225.86, p-value = 0.261
#95 percent confidence interval: -1.003345  0.273312 
#ybarT - ybarC = -0.353768
se = abs((-1.003345 - -0.353768) / 1.96)
-0.353768 + 1.96 * se
-0.353768 - 1.96 * se


variance = se^2
variance_with_our_method = variance / 1.2
se_with_our_method = sqrt(variance_with_our_method)
-0.353768 + 1.96 * se_with_our_method
-0.353768 - 1.96 * se_with_our_method





###BART AND RF
set_bart_machine_num_cores(4)
bart_machine = build_bart_machine(Xy = Xy,
		num_trees = 200,
		num_burn_in = 1000,
		num_iterations_after_burn_in = 1000)
bart_machine
plot_y_vs_yhat(bart_machine)
plot_convergence_diagnostics(bart_machine)

library(randomForest)
rf = randomForest(x = Xy[, 1 : (ncol(Xy) - 1)], y = Xy$y)
#find R^2
yhat = predict(rf, Xy[, 1 : (ncol(Xy) - 1)])
plot(Xy$y, yhat, ylim = c(-10, 4))
abline(a = 0, b = 1)
sse = sum((Xy$y - yhat)^2)
sst = sum((Xy$y - mean(Xy$y))^2)
Pseudo_Rsq_rf = 1 - sse / sst
Pseudo_Rsq_rf






Nsim = 200
n_study = 224 ### the pretend amount of people in this study (224 total)
prob_match_cutoff_alpha = 0.10 ###found to be the best in the basic simulations

beta_matches = array(NA, Nsim)
beta_se_matches = array(NA, Nsim)
beta_regression = array(NA, Nsim)
beta_se_regression = array(NA, Nsim)
sub_ns = array(NA, Nsim)

Xy_controlT = Xy[Xy$armbin == 1, ]
Xy_controlC = Xy[Xy$armbin == 0, ]

for (nsim in 1 : Nsim){
	cat("nsim:", nsim," ")
	#scramble the order of the participants for each run
	X = Xy[sample(1 : nrow(Xy), n_study), ]
	#pull out true indicator of treatment and true y
	true_trt = X$armbin
	y_true = X$y
	X$y = NULL
	X$armbin = NULL
	
	#make model matrix (turn factor cols into dummy cols)
	Xmodel = as.data.frame(model.matrix(lm(y_true ~ 
							age + sex + sextrans + uti + white + 
							hisp + education + employment + b_frequency +
							b_SIQTOT + b_HAD15 + b_MCS + b_PIQTOT + b_WISIC +
							b_pain + b_urgency +
							partner, X)))
#	Xmodel = as.data.frame(model.matrix(lm(y_true ~ partner + b_pain + b_WISIC + b_frequency, X)))
	
	
	
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
	Xyused = Xmodel[!is.na(indic_T), ]
	sub_ns[nsim] = length(indic_T_sub)
	
	cat("  sub_n:", sub_ns[nsim])
	#now use our estimator
	#get reservoir data
	Xyleft = Xyused[match_indic_sub == 0, ]
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
#	cor_test = cor.test(matched_ts, matched_cs)
#	cat("  r: ", cor_test$estimate, "pval:", cor_test$p.value, "\n")
	cat("\n")	

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
	
	linear_mod = lm(y ~ armbin, X) #same as 2-sample t-test with pooled variance assumption
	beta_regression[nsim] = coef(summary(linear_mod))[2, 1]
	beta_se_regression[nsim] = coef(summary(linear_mod))[2, 2]
}

#xs_to_date_match = cbind(xs_to_date, match_indic)
#xs_to_date_match[order(xs_to_date_match$match_indic), ]

mean(beta_se_matches^2, na.rm=T)
mean(beta_se_regression^2, na.rm=T)
mean(beta_se_matches, na.rm=T)
mean(beta_se_regression, na.rm=T)

#efficiency via matching - paired test due to sub_n being different run-run
mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
###approx sample size reduction due to method
1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1

t.test(beta_se_matches, beta_se_regression)
mean(sub_ns)


############## models for (y_true ~ 
#	age + sex + sextrans + uti + white + 
#	hisp + education + employment + b_frequency +
#	b_SIQTOT + b_HAD15 + b_MCS + b_PIQTOT + b_WISIC +
#	b_freq1 + b_pain + b_urgency +
#	b_pain1 + b_urgency1 + partner

#> mean(beta_se_matches^2, na.rm=T)
#[1] 0.5705806
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.6429635
#> mean(beta_se_matches, na.rm=T)
#[1] 0.7448645
#> mean(beta_se_regression, na.rm=T)
#[1] 0.7966292
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.298644
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.2299658
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = -4.7035, df = 363.614, p-value = 3.637e-06
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.07340702 -0.03012235 
#sample estimates:
#		mean of x mean of y 
#0.7448645 0.7966292 
#
#> mean(sub_ns)
#[1] 38.905


#> mean(beta_se_matches^2, na.rm=T)
#[1] 0.3113475
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.3263608
#> mean(beta_se_matches, na.rm=T)
#[1] 0.5542728
#> mean(beta_se_regression, na.rm=T)
#[1] 0.5698078
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.101373
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.0920428
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = -2.8753, df = 337.923, p-value = 0.004292
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.026162591 -0.004907465 
#sample estimates:
#		mean of x mean of y 
#0.5542728 0.5698078 
#
#> mean(sub_ns)
#[1] 75.2

#> mean(beta_se_matches^2, na.rm=T)
#[1] 0.2145067
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.2188113
#> mean(beta_se_matches, na.rm=T)
#[1] 0.4611856
#> mean(beta_se_regression, na.rm=T)
#[1] 0.4671154
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.051974
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.04940643
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = -1.6972, df = 319.915, p-value = 0.09062
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.0128035553  0.0009439116 
#sample estimates:
#		mean of x mean of y 
#0.4611856 0.4671154 
#
#> mean(sub_ns)
#[1] 111.255


#		> mean(beta_se_matches^2, na.rm=T)
#[1] 0.1405304
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.1483538
#> mean(beta_se_matches, na.rm=T)
#[1] 0.374056
#> mean(beta_se_regression, na.rm=T)
#[1] 0.3849709
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.07204
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.06719932
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = -5.5718, df = 291.616, p-value = 5.736e-08
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.014770305 -0.007059348 
#sample estimates:
#		mean of x mean of y 
#0.3740560 0.3849709 
#
#> mean(sub_ns)
#[1] 165.545




############## models for y_true ~ partner + b_pain + b_WISIC + b_frequency

#> mean(beta_se_matches^2, na.rm=T)
#[1] 0.6007976
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.6599561
#> mean(beta_se_matches, na.rm=T)
#[1] 0.7633119
#> mean(beta_se_regression, na.rm=T)
#[1] 0.8070775
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.2693
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.2121639
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = -3.7759, df = 352.79, p-value = 0.000187
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.06656111 -0.02097000 
#sample estimates:
#		mean of x mean of y 
#0.7633119 0.8070775 
#
#> mean(sub_ns)
#[1] 38.01


#> mean(beta_se_matches^2, na.rm=T)
#[1] 0.2846047
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.3357241
#> mean(beta_se_matches, na.rm=T)
#[1] 0.5304094
#> mean(beta_se_regression, na.rm=T)
#[1] 0.5778175
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.231276
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.1878346
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = -9.3449, df = 369.608, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.05738397 -0.03743217 
#sample estimates:
#		mean of x mean of y 
#0.5304094 0.5778175 
#
#> mean(sub_ns)
#[1] 72.86

#
#> #xs_to_date_match = cbind(xs_to_date, match_indic)
#		> #xs_to_date_match[order(xs_to_date_match$match_indic), ]
#		> 
#		> mean(beta_se_matches^2, na.rm=T)
#[1] 0.2019438
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.2258122
#> mean(beta_se_matches, na.rm=T)
#[1] 0.447516
#> mean(beta_se_regression, na.rm=T)
#[1] 0.4744269
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.152904
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.1326252
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = -7.7416, df = 345.059, p-value = 1.099e-13
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.03374802 -0.02007378 
#sample estimates:
#		mean of x mean of y 
#0.4475160 0.4744269 
#
#> mean(sub_ns)
#[1] 108.55

#
#[1] 0.1370412
#> mean(beta_se_regression^2, na.rm=T)
#[1] 0.1527293
#> mean(beta_se_matches, na.rm=T)
#[1] 0.3695245
#> mean(beta_se_regression, na.rm=T)
#[1] 0.3905596
#> 
#		> #efficiency via matching - paired test due to sub_n being different run-run
#		> mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T)
#[1] 1.128392
#> ###approx sample size reduction due to method
#		> 1 - (mean(beta_se_regression^2 / beta_se_matches^2, na.rm=T))^-1
#[1] 0.1137832
#> 
#		> t.test(beta_se_matches, beta_se_regression)
#
#Welch Two Sample t-test
#
#data:  beta_se_matches and beta_se_regression 
#t = -11.3352, df = 333.853, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#		-0.02468558 -0.01738473 
#sample estimates:
#		mean of x mean of y 
#0.3695245 0.3905596 
#
#> mean(sub_ns)
#[1] 160.33






#select 4 best variables
rf_mod = randomForest(Xy[, 1 : (ncol(Xy) - 1)], y = Xy$y, importance = TRUE)
varImpPlot(rf_mod, type = 1)


head(Xy)
Xyred = Xy[, c("b_pain", "partner", "b_WISIC", "b_frequency", "y")]

modred = lm(y ~ ., Xyred)
summary(modred)

###BART AND RF
set_bart_machine_num_cores(4)
bart_machine = build_bart_machine(Xy = Xyred,
		num_trees = 200,
		num_burn_in = 1000,
		num_iterations_after_burn_in = 1000)
bart_machine
plot_y_vs_yhat(bart_machine)
plot_convergence_diagnostics(bart_machine)

library(randomForest)
rf = randomForest(x = Xyred[, 1 : (ncol(Xyred) - 1)], y = Xyred$y)
#find R^2
yhat = predict(rf, Xyred[, 1 : (ncol(Xyred) - 1)])
plot(Xyred$y, yhat, ylim = c(-10, 4))
abline(a = 0, b = 1)
sse = sum((Xyred$y - yhat)^2)
sst = sum((Xyred$y - mean(Xyred$y))^2)
Pseudo_Rsq_rf = 1 - sse / sst
Pseudo_Rsq_rf
