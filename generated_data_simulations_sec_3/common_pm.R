#do the standard matching procedure

#initialize the indicator treatment vector
indic_T = array(NA, n) 
#initialize the reservoir
match_indic = array(-1, n) #0 indicates reservoir

#now we're going to go through and do the matching
for (i_match in 1 : n){
	#if there is nothing in the reservoir, randomize and add it to reservoir
	if (length(match_indic[match_indic == 0]) == 0 || i_match <= p + 1){
		indic_T[i_match] = rbinom(1, 1, prob_trt)
		match_indic[i_match] = 0
	} else {
		#first calculate the threshold we're operating at
		xs_to_date = x_s[1 : i_match, ]
		S_xs_inv = solve(var(xs_to_date))
		F_crit =  qf(prob_match_cutoff_alpha, p, i_match - p)
		T_cutoff_sq = p * (n - 1) / (n - p) * F_crit
		#now iterate over all items in reservoir and take the minimum distance x
		reservoir_indices = which(match_indic == 0)
		x_star = x_s[i_match, ]
		sqd_distances = array(NA, length(reservoir_indices))
		for (r in 1 : length(reservoir_indices)){
			sqd_distances[r] = 1 / 2 * 
					t(x_star - x_s[reservoir_indices[r], ]) %*%
					S_xs_inv %*%
					(x_star - x_s[reservoir_indices[r], ])			
		}
#		cat(paste("i_match", i_match, "sqd_distances", paste(sqd_distances, collapse = ", "), "T_cutoff_sq", T_cutoff_sq, "\n"))
		
		#find minimum distance index
		min_sqd_dist_index = which(sqd_distances == min(sqd_distances))
		if (length(sqd_distances[min_sqd_dist_index]) > 1 || length(T_cutoff_sq) > 1){
			min_sqd_dist_index = min_sqd_dist_index[1] #if there's a tie, just take the first one
		}
		#if it's smaller than the threshold, we're in business: match it
		if (sqd_distances[min_sqd_dist_index] < T_cutoff_sq){
			match_num = max(match_indic) + 1
			match_indic[reservoir_indices[min_sqd_dist_index]] = match_num
			match_indic[i_match] = match_num
			indic_T[i_match] = 1 - indic_T[reservoir_indices[min_sqd_dist_index]]
		} else { #otherwise, randomize and add it to the reservoir
			indic_T[i_match] = rbinom(1, 1, prob_trt)
			match_indic[i_match] = 0		
		}
	}
}

#print(table(match_indic))

#create response vector
source(paste("create_response_", response_model, ".R", sep = ""))


#create design matrix
Xy = as.data.frame(cbind(x_s, indic_T, y, match_indic))
#write.csv(Xy, "matched_matrix.csv", row.names = FALSE)
#pull out yT, yC, xT, xC
yTs = Xy[Xy$indic_T == 1, "y"]
yCs = Xy[Xy$indic_T == 0, "y"]
xTs = Xy[Xy$indic_T == 1, 1 : p]
xCs = Xy[Xy$indic_T == 0, 1 : p]

#and record other useful info
final_reservoir_size[nsim] = length(match_indic[match_indic == 0]) / n
m = max(match_indic)

#now we want to just calculate the diffs inside matches and ignore the reservoir
ydiffs = array(NA, m)

for (match_id in 1 : m){
	yT = Xy[Xy$indic_T == 1 & Xy$match_indic == match_id, "y"]
	yC = Xy[Xy$indic_T == 0 & Xy$match_indic == match_id, "y"]
	ydiffs[match_id] = yT - yC
}

#the estimate is the one sample t test from Stat 101
if (length(ydiffs) <= 1){
	d_bar = 0
	beta_Ts[nsim] = d_bar
	pval = 1
	Rsqs[nsim] = 0
} else {
	d_bar = mean(ydiffs)
	ttest = t.test(ydiffs)
	beta_Ts[nsim] = d_bar
	pval = ttest$p.value
	Rsqs[nsim] = 1 - sum((ydiffs - beta_Ts[nsim])^2) / sum((ydiffs - mean(ydiffs))^2)
}