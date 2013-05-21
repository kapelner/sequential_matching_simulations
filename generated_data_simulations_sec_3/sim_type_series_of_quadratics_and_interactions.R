response_model = "linear_quad_and_interact_medley"

#beta_1s = c(
#	0,
#	1,
#	2,
#	3
#)
#beta_2s = c(
#	0,
#	2
#)
#beta_3s = c(
#	0, 
#	0.5, 
#	1
#)
#beta_4s = c(
#	0, 
#	0.5, 
#	1
#)
#beta_5s = c(
#	0, 
#	0.5, 
#	1
#)
#
##create permutation matrix of betas so we can pull out one row and parallelize the script by iteration number
#beta_mat = matrix(NA, nrow = length(beta_1s) * length(beta_2s) * length(beta_3s) * length(beta_4s) * length(beta_5s), ncol = 5)
#i_n = 1
#
#for (beta_1 in beta_1s){
#	for (beta_2 in beta_2s){
#		for (beta_3 in beta_3s){
#			for (beta_4 in beta_4s){
#				for (beta_5 in beta_5s){
#					beta_mat[i_n, ] = c(beta_1, beta_2, beta_3, beta_4, beta_5)
#					i_n = i_n + 1
#				}	
#			}
#		}			
#	}
#}
#

beta_mat = matrix(NA, nrow = 6, ncol = 5)
beta_mat[1, ] = c(0,0,0,0,0)
beta_mat[2, ] = c(2,0,0,0,0)
beta_mat[3, ] = c(2,2,0,0,0)
beta_mat[4, ] = c(2,1,1,1,0.5)
beta_mat[5, ] = c(2,2,-1,-1,1)
beta_mat[6, ] = c(1,1,1,1,1)

beta_1 = beta_mat[iter_num, 1]
beta_2 = beta_mat[iter_num, 2]
beta_3 = beta_mat[iter_num, 3]
beta_4 = beta_mat[iter_num, 4]
beta_5 = beta_mat[iter_num, 5]


					
#each prob_match and response gets their own results matrix and it gets saved to a different CSV file
results = matrix(NA, nrow = length(metrics_for_each_run), ncol = length(randomization_types))
colnames(results) = randomization_types
rownames(results) = array(NA, length(metrics_for_each_run))
#add row names to the results matrix
for (m in 1 : length(metrics_for_each_run)){
	rownames(results) = metrics_for_each_run	
}	

full_results_across_alpha = matrix(NA, nrow = 0, ncol = ncol(results) + 2)
colnames(full_results_across_alpha) = c("alpha", "n", randomization_types)

#now do the simulations
for (i_n in 1 : length(ns_to_test)){	
	for (prob_match_cutoff_alpha in prob_match_cutoff_alphas){
		
		n = ns_to_test[i_n]	
	
		for (randomization_type in randomization_types){
			#don't do duplicate runs if you don't have to - the alphas are ONLY for pair match runs
			if (length(grep("pm", randomization_type)) == 0 & prob_match_cutoff_alpha != prob_match_cutoff_alphas[1]){
				next
			}
			
			cat(randomization_type, "n =", n, "alpha =", prob_match_cutoff_alpha, " betas ", beta_1, " ", beta_2, " ", beta_3, " ", beta_4, " ", beta_5, "\n")
			
			source("inner_simulation.R")
		}
		full_results_across_alpha = rbind(
			full_results_across_alpha, 
			cbind(
				rep(prob_match_cutoff_alpha, nrow(results)), 
				rep(n, nrow(results)),
				results
			)
		)
		#save the results to the hard drive
		write.csv(full_results_across_alpha, paste(RESULTS_FILE_HEAD, beta_1, "_", beta_2, "_", beta_3, "_", beta_4, "_", beta_5, ".csv", sep = ""))		
	}
	
}
				