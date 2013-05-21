#get reservoir data
Xyleft = Xy[Xy$match_indic == 0, ] #pull out the reservoir data
final_reservoir_size[nsim] = nrow(Xyleft) / n

YleftT = Xyleft[Xyleft$indic_T == 1, ]$y #get the reservoir responses from the treatment
YleftC = Xyleft[Xyleft$indic_T == 0, ]$y #get the reservoir responses from the control
YbarRTMinusYbarRC = mean(YleftT) - mean(YleftC) #compute the classic estimator from the reservoir: ybar_T - ybar_C

#get reservoir sample sizes
nR = nrow(Xyleft) #how many observations are there in the reservoir?
nRT = length(YleftT) #how many treatment observations are there in the reservoir?
nRC = length(YleftC) #how many control observations are there in the reservoir?

ssqR = (var(YleftT) * (nRT - 1) + var(YleftC) * (nRC - 1)) / (nR - 2) * (1 / nRT + 1 / nRC)
ssqD_bar = var(ydiffs) / m

gamma_star_hat = ssqR / (ssqR + ssqD_bar)

Xymatched = Xy[Xy$match_indic > 0, ]
Xymatched = Xymatched[order(Xymatched$match_indic), ]

Xyleft = Xyleft[order(Xyleft$indic_T), ]