setwd("c:/users/kapelner/workspace/StatTurk/matching_sims")
POW = read.csv("all_powers.csv", row.names = NULL)

POW

n_50_col = "red"
n_100_col = "green"
n_200_col = "blue"
offset_classic = 2
offset_linear = 7
offset_exact = 12
offset_mod_I = 1
offset_mod_II = 4
offset_mod_III = 7
PCH_symbol = 16
METHOD_NAMES = c("C", "E", "S", "M", "SM")
CEX_SIZE = 2
GEN_FOR_FILES = TRUE
MARGINS_INNER = c(2.5,2.5,0.1,0.5)
MARGINS = c(2.5,4.5,0.1,0.5)
DOT_SIZE = 2
LINE_THICKNESS = 4

graphics.off()
par(mar=MARGINS)

########MODEL I Classic
if (GEN_FOR_FILES){
	pdf("model_I_classic.PDF")
	par(mar=MARGINS)
} else {
	windows()
}
plot(1:5, rep(0, 5), ylim = c(0, 1), xlab = "", ylab = "Power", cex.lab = CEX_SIZE, type = "n", xaxt = "n", cex.axis = CEX_SIZE)
axis(1, at = 1:5, labels = METHOD_NAMES, cex.axis = CEX_SIZE)
points(1:5, POW[offset_mod_I, (offset_classic + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_50_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[offset_mod_I, offset_classic + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_50_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_I + 1), (offset_classic + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_100_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_I + 1), offset_classic + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_100_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_I + 2), (offset_classic + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_200_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_I + 2), offset_classic + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_200_col, lwd = LINE_THICKNESS)
}
if (GEN_FOR_FILES){
	dev.off()
}


########MODEL I Linear
if (GEN_FOR_FILES){
	pdf("model_I_linear.PDF")
	par(mar=MARGINS_INNER)
} else {
	windows()
}
plot(1:5, rep(0, 5), ylim = c(0, 1), xlab = "", ylab = "Power", type = "n", xaxt = "n", cex.axis = CEX_SIZE)
axis(1, at = 1:5, labels = METHOD_NAMES, cex.axis = CEX_SIZE)
points(1:5, POW[offset_mod_I, (offset_linear + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_50_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[offset_mod_I, offset_linear + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_50_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_I + 1), (offset_linear + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_100_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_I + 1), offset_linear + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_100_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_I + 2), (offset_linear + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_200_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_I + 2), offset_linear + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_200_col, lwd = LINE_THICKNESS)
}
if (GEN_FOR_FILES){
	dev.off()
}



########MODEL I Exact
if (GEN_FOR_FILES){
	pdf("model_I_exact.PDF")
	par(mar=MARGINS_INNER)
} else {
	windows()
}
plot(1:5, rep(0, 5), ylim = c(0, 1), xlab = "", ylab = "Power", type = "n", xaxt = "n", cex.axis = CEX_SIZE)
axis(1, at = 1:5, labels = METHOD_NAMES, cex.axis = CEX_SIZE)
points(1:5, POW[offset_mod_I, (offset_exact + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_50_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[offset_mod_I, offset_exact + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_50_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_I + 1), (offset_exact + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_100_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_I + 1), offset_exact + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_100_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_I + 2), (offset_exact + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_200_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_I + 2), offset_exact + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_200_col, lwd = LINE_THICKNESS)
}
if (GEN_FOR_FILES){
	dev.off()
}




########MODEL II Classic
if (GEN_FOR_FILES){
	pdf("model_II_classic.PDF")
	par(mar=MARGINS)
} else {
	windows()
}
plot(1:5, rep(0, 5), ylim = c(0, 1), xlab = "", ylab = "Power", cex.lab = CEX_SIZE, type = "n", xaxt = "n", cex.axis = CEX_SIZE)
axis(1, at = 1:5, labels = METHOD_NAMES, cex.axis = CEX_SIZE)
points(1:5, POW[offset_mod_II, (offset_classic + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_50_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[offset_mod_II, offset_classic + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_50_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_II + 1), (offset_classic + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_100_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_II + 1), offset_classic + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_100_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_II + 2), (offset_classic + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_200_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_II + 2), offset_classic + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_200_col, lwd = LINE_THICKNESS)
}
if (GEN_FOR_FILES){
	dev.off()
}



########MODEL II Linear
if (GEN_FOR_FILES){
	pdf("model_II_linear.PDF")
	par(mar=MARGINS_INNER)
} else {
	windows()
}
plot(1:5, rep(0, 5), ylim = c(0, 1), xlab = "", ylab = "Power", type = "n", xaxt = "n", cex.axis = CEX_SIZE)
axis(1, at = 1:5, labels = METHOD_NAMES, cex.axis = CEX_SIZE)
points(1:5, POW[offset_mod_II, (offset_linear + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_50_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[offset_mod_II, offset_linear + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_50_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_II + 1), (offset_linear + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_100_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_II + 1), offset_linear + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_100_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_II + 2), (offset_linear + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_200_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_II + 2), offset_linear + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_200_col, lwd = LINE_THICKNESS)
}
if (GEN_FOR_FILES){
	dev.off()
}



########MODEL II Exact
if (GEN_FOR_FILES){
	pdf("model_II_exact.PDF")
	par(mar=MARGINS_INNER)
} else {
	windows()
}
plot(1:5, rep(0, 5), ylim = c(0, 1), xlab = "", ylab = "Power", type = "n", xaxt = "n", cex.axis = CEX_SIZE)
axis(1, at = 1:5, labels = METHOD_NAMES, cex.axis = CEX_SIZE)
points(1:5, POW[offset_mod_II, (offset_exact + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_50_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[offset_mod_II, offset_exact + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_50_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_II + 1), (offset_exact + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_100_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_II + 1), offset_exact + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_100_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_II + 2), (offset_exact + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_200_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_II + 2), offset_exact + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_200_col, lwd = LINE_THICKNESS)
}
if (GEN_FOR_FILES){
	dev.off()
}



########MODEL III Classic
if (GEN_FOR_FILES){
	pdf("model_III_classic.PDF")
	par(mar=MARGINS)
} else {
	windows()
}
plot(1:5, rep(0, 5), ylim = c(0, 1), xlab = "", ylab = "Power", cex.lab = CEX_SIZE, type = "n", xaxt = "n", cex.axis = CEX_SIZE)
axis(1, at = 1:5, labels = METHOD_NAMES, cex.axis = CEX_SIZE)
points(1:5, POW[offset_mod_III, (offset_classic + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_50_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[offset_mod_III, offset_classic + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_50_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_III + 1), (offset_classic + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_100_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_III + 1), offset_classic + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_100_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_III + 2), (offset_classic + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_200_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_III + 2), offset_classic + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_200_col, lwd = LINE_THICKNESS)
}
if (GEN_FOR_FILES){
	dev.off()
}



########MODEL III Linear
if (GEN_FOR_FILES){
	pdf("model_III_linear.PDF")
	par(mar=MARGINS_INNER)
} else {
	windows()
}
plot(1:5, rep(0, 5), ylim = c(0, 1), xlab = "", ylab = "Power", type = "n", xaxt = "n", cex.axis = CEX_SIZE)
axis(1, at = 1:5, labels = METHOD_NAMES, cex.axis = CEX_SIZE)
points(1:5, POW[offset_mod_III, (offset_linear + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_50_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[offset_mod_III, offset_linear + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_50_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_III + 1), (offset_linear + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_100_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_III + 1), offset_linear + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_100_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_III + 2), (offset_linear + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_200_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_III + 2), offset_linear + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_200_col, lwd = LINE_THICKNESS)
}
if (GEN_FOR_FILES){
	dev.off()
}



########MODEL III Exact
if (GEN_FOR_FILES){
	pdf("model_III_exact.PDF")
	par(mar=MARGINS_INNER)
} else {
	windows()
}
plot(1:5, rep(0, 5), ylim = c(0, 1), xlab = "", ylab = "Power", type = "n", xaxt = "n", cex.axis = CEX_SIZE)
axis(1, at = 1:5, labels = METHOD_NAMES, cex.axis = CEX_SIZE)
points(1:5, POW[offset_mod_III, (offset_exact + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_50_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[offset_mod_III, offset_exact + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_50_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_III + 1), (offset_exact + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_100_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_III + 1), offset_exact + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_100_col, lwd = LINE_THICKNESS)
}
points(1:5, POW[(offset_mod_III + 2), (offset_exact + 1:5)], ylim = c(0, 1), cex = DOT_SIZE, xlab = "", ylab = "Power", col = n_200_col, pch = PCH_symbol)
for (i in 1 : 5){
	mid = POW[(offset_mod_III + 2), offset_exact + i]
	moe = 1.96 * sqrt(mid * (1 - mid) / 1000)
	segments(i, mid - moe, i, mid + moe, col = n_200_col, lwd = LINE_THICKNESS)
}
if (GEN_FOR_FILES){
	dev.off()
}