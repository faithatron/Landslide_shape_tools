#############################
# Name: Inverse Gamma Landslide Fit 
# Author: Faith Taylor
# Last updated : 02/02/2014
# Purpose: Fit inverse gamma pdf to landslide data and test gooness of fit
# Copyright: (c) 2014 Faith Taylor
# This program is free software under the GNU General Public Licence (>=v2). 
##############################



# Load modules
library(MCMCpack)
library(MASS)
library(PearsonDS)
library(calibrate)
library(sfmisc)
# Initial definitions 
# 1. The inverse gamma log likelihood function

ll = function(par){
  if(par[1]>0 & par[2]>0 & par[3]<min(data)) return( -sum(log(dinvgamma(data-    par[3],par[1],par[2]))) )
  else return(Inf)
}
# 2. Make a series of x values (Length to width ratios)
xs = seq(1,100, 0.01)

# 3. Table of values for bin sizes
#lower_upper_vals = read.csv("lower_upper_vals.csv")

rm(A)
#parmar x x top right-side
par(mar= c(8,5,1,1), family = "serif")
LWRrainbow_pal = c("#902525", "#ff0707", "#f18c37", "#fcfc02", "#91cd4e", "#1EBB66", "#07b1f0", "#0069bd", "#001055", "#68259b")


lower = lower_upper_vals[1,2]
upper = lower_upper_vals[1,3]
phrase1 = bquote(bold(.(lower)~ m^2 ~""<= "" ~bolditalic(A[L]) ~ "<" ~bold(.(upper) ~ m^2)))
lower = lower_upper_vals[2,2]
upper = lower_upper_vals[2,3]
phrase2 = bquote(bold(.(lower)~ m^2 ~""<= "" ~bolditalic(A[L]) ~ "<" ~bold(.(upper) ~ m^2)))
lower = lower_upper_vals[3,2]
upper = lower_upper_vals[3,3]
phrase3 = bquote(bold(.(lower)~ m^2 ~""<= "" ~bolditalic(A[L]) ~ "<" ~bold(.(upper) ~ m^2)))
lower = lower_upper_vals[4,2]
upper = lower_upper_vals[4,3]
phrase4 = bquote(bold(.(lower)~ m^2 ~""<= "" ~bolditalic(A[L]) ~ "<" ~bold(.(upper) ~ m^2)))
lower = lower_upper_vals[5,2]
upper = lower_upper_vals[5,3]
phrase5 = bquote(bold(.(lower)~ m^2 ~""<= "" ~bolditalic(A[L]) ~ "<" ~bold(.(upper) ~ m^2)))
lower = lower_upper_vals[6,2]
upper = lower_upper_vals[6,3]
phrase6 = bquote(bold(.(lower)~ m^2 ~""<= "" ~bolditalic(A[L]) ~ "<" ~bold(.(upper) ~ m^2)))
lower = lower_upper_vals[7,2]
upper = lower_upper_vals[7,3]
phrase7 = bquote(bold(.(lower)~ m^2 ~""<= "" ~bolditalic(A[L]) ~ "<" ~bold(.(upper) ~ m^2)))
lower = lower_upper_vals[8,2]
upper = lower_upper_vals[8,3]
phrase8 = bquote(bold(.(lower)~ m^2 ~""<= "" ~bolditalic(A[L]) ~ "<" ~bold(.(upper) ~ m^2)))
lower = lower_upper_vals[9,2]
upper = lower_upper_vals[9,3]
phrase9 = bquote(bold(.(lower)~ m^2 ~ ""<= "" ~bolditalic(A[L]) ~ "<" ~bold(.(upper) ~ m^2)))
lower = lower_upper_vals[10,2]
upper = lower_upper_vals[10,3]
phrase10 = bquote(bold( ~bolditalic(A[L]) ~ ""<= "" ~bold(.(lower) ~ m^2)))


roundUp <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

#par(mfrow=c(4,3))
#alldata = northridge_new_lwrs[!is.na(northridge_new_lwrs)]
alldata = guat_data[!is.na(guat_data)]
max_data = max(alldata)
round_max_data = roundUp(max_data)

#plot(1,1,   col = "White",  xaxs ="i", yaxs = "i", ylim = c(0,20), xlim = c(0,35))
#for (h in 1:10){
#  boxplot(c(Umbria_lwrs[h], guat_elip_binlwrs[h], northridge_new_lwrs[h]), range = 0, lwd = 1, lty = 1, col = LWRrainbow_pal[h],  xaxs ="i", yaxs = "i", add = TRUE, at = c((1+(h*3)),(2+(h*3)),(3+(h*3))))
#}

# Create boxplot of all data by category
par(mar= c(10,5,1,1), family = "serif")
boxplot(nridge_data, range = 0, lwd = 1, lty = 1, col = LWRrainbow_pal, axes = FALSE, xaxs ="i", yaxs = "i", ylim = c(0,round_max_data+5))
axis(1, tck = 0, at = seq(0,11,1), labels = FALSE, lwd = 2)
axis(2, lwd = 2, at = seq(0, round_max_data+5, 5), las = 1, cex.axis = 1.5, tck = 0.03)
axis(3, tck = 0, labels = FALSE, at = c(0,11), lwd = 2)
axis(4, lwd = 2, at = seq(0, round_max_data+5, 5), las = 1, cex.axis = 1.5, labels = FALSE, tck = 0.03)
title(ylab = expression(bold("Length to Width Ratio, ") ~bolditalic("LWR")), cex.lab = 1.5, line = 3)
title(xlab = expression(bold("Landslide Area Category, ") ~bolditalic(C[A[L] ]) ~bold((m^2))), cex.lab = 1.5, line = 9)
labels = list( phrase1, phrase2, phrase3, phrase4, phrase5, phrase6, phrase7, phrase8, phrase9, phrase10)
text(1:10, par("usr")[3] - 0.25, srt = 45, adj = 1,labels = do.call("expression", labels), xpd = TRUE)
for (g in 1:10){
  dataset = nridge_data[g]
  max_cat = max(dataset[!is.na(dataset)]) + 1.5
  num = length(dataset[!is.na(dataset)])
  faith = bquote(bolditalic("n = " ~bold(.(num))))
  text(g,max_cat, srt = 90, labels = faith, col = "grey")
  med = median(dataset[!is.na(dataset)])
  print(med)
}

# Read data
for (A in 1:10){
  rm(shape_f, location_f, scale_f, Proportion, all_d, quants, m, ests, estsa, b, bin, bin_boundaries, binsizes, bootstrap_location, bootstrap_location, bootstrap_scale, bootstrap_scale_recip, bootstrap_shape, cumulative_x, data, dif, dval, emp_bin_mids, emp_cum_y, emp_freq_dens, emp_freqs, emp_prob_dens, empirical_D, g, h, i, k, label, loc_max, location_s, logdat, loghist, lower, lower_pdf, max_pdf, max_x_d, maxx, median_pdf,  mle, n, numbins, params, parset, pdf, phrase, phrase2, rand_D, rand_cdf, rand_lwrs, rand_lwrs_cumulative, location_1, scale_a, shape_a, shape_rho, smooth_teo_cum, startlocation, startscale, startshape, teo_cdf, teo_cum, testresult, upper, upper_pdf, val, yupper, yupper2)
  
  # 1. Initialising steps
  # Name the bin category (e.g. LW1 is landslides between 0 and 100 m2) 
  bin = paste("Northridge","LW", A, sep = "")
  # Add in here the text for lanslide bin size
  lower = lower_upper_vals[A,2]
  upper = lower_upper_vals[A,3]
  phrase = bquote(bold(.(lower)~ m^2 ~"<" ~bolditalic(A[L]) ~ "<" ~bold(.(upper) ~ m^2)))
  variablej = paste("phrase_", A, sep = "")
  assign(variablej, phrase)
  # Read a column of the table
  data = northridge_new_lwrs[A]
  #Strip out any non numeric rows
  data = data[!is.na(data)]
  #Calculate the number of observations
  n = length(data)
  phrase2 = bquote(bolditalic( "n = "~.(n)))

  
  # 2. Fitting an inverse gamma pdf to the raw data 
  # Use mle fitting to find the parameters that best fit the lwrs in this area category. This uses the optim function
  try((mle = optim(c(1,1,1),ll)), silent = FALSE)

  # Obtain initial parameter values 
  params = mle$par
  shape_rho = params [1]
  scale_a = params[2]
  location_s = params[3]


  # Produce a plot of the observed data and the pdf/cdf fit to the data (non cumulative then cumulative)
  # Evaluate the inverse gamma pdf for each of the x values using the parameters from MLE fitting
  pdf = dpearsonV(xs, shape_rho, location_s, (1/scale_a))

  # Bin the raw data 

  logdat = log(data)
  loghist = hist(logdat, plot = FALSE)
  bin_boundaries = exp(loghist$breaks)
  emp_bin_mids = exp(loghist$mids)
  emp_freqs = loghist$counts
  #calculate bin sizes
  numbins = length(emp_bin_mids)
  binsizes = numeric(numbins)
  for (k in 1:numbins){
    dif = bin_boundaries[k +1] - bin_boundaries[k]
    binsizes[k] = dif}
  
  emp_freq_dens = emp_freqs/binsizes
  emp_prob_dens = emp_freqs/(binsizes * n)
  

  # Plot the pdf and raw data
  par(mar= c(5,6,2,2))
  plot(xs, pdf, type = "l", col = rgb(1,0,0), lwd = 3, ylim= c(10^-8, 10), log = "xy", ylab = "", xlab=expression(bold("Length to Width Ratio, "~bolditalic(LWR))), axes = FALSE, yaxs ="i", xaxs = "i", cex.lab = 1.7 )
  ticks = c(0.1,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)
  title(ylab = expression(bold("Probability Density, " ~bolditalic(p(LWR)))), cex.lab = 1.7, line = 4)
  axis(1, at = c(1,10, 100), labels =( expression(~bold(10^0),~bold(10^1), ~bold(10^2))), lty = 1, lwd = 0, tck = 0, cex.axis = 1.7)
  yticks = c(0.000000001, 0.000000002, 0.000000003, 0.000000004, 0.000000005,  0.000000006, 0.000000007, 0.000000008, 0.000000009, 0.00000001, 0.00000002, 0.00000003, 0.00000004, 0.00000005,  0.00000006, 0.00000007, 0.00000008, 0.00000009,0.0000001, 0.0000002, 0.0000003, 0.0000004, 0.0000005,  0.0000006, 0.0000007, 0.0000008, 0.0000009, 0.000001, 0.000002, 0.000003,  0.000004, 0.000005, 0.000006, 0.000007, 0.000008, 0.000009, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.00008, 0.00009,0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10)
  axis(2, at = c(0.000000001, 0.00000001,0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1,10), labels =(expression(~bold(10^-9),~bold(10^-8),~bold(10^-7),~bold(10^-6),~bold(10^-5), ~bold(10^-4), ~bold(10^-3), ~bold(10^-2), ~bold(10^-1), ~bold(10^0) , ~bold(10^1))), lty = 1, lwd = 0, tck = 0, cex.axis = 1.7, las = 1)
  axis(1, ticks, labels = NA, lty = 1, lwd = 2, tck = 0.01)
  axis(2, yticks, labels = NA, lty = 1, lwd = 2, tck = 0.01)
  axis(3, ticks, labels = NA, lty = 1, lwd = 2, tck = 0.01)
  axis(4, yticks, labels = NA, lty = 1, lwd = 2, tck = 0.01)
  points(emp_bin_mids, emp_prob_dens, pch = 20)
  val <- substitute("value" == phrase, list(phrase = phrase)) 
  legend("topright", do.call("expression", list(phrase, "Observed", "Inverse Gamma pdf", phrase2)), lty = c(NA, NA, 1), col = c(NA, "black", "red"), pch = c(NA, 20, NA), lwd = c(NA, NA, 2),bty = 'o', cex = 1.1, box.lwd = 2, inset = 0.01)

  ## Now plot the cumulative
  # Sort the observed data
  cumulative_x = sort(data)
  emp_cum_y = (1:n)/n
  emp_cum_y = sort(emp_cum_y, decreasing = TRUE)
  teo_cum = ppearsonV(cumulative_x, shape_rho, location_s, (1/scale_a), lower.tail = FALSE)
  empirical_D = max(abs(emp_cum_y - teo_cum))
  smooth_teo_cum = ppearsonV(xs, shape_rho, location_s, (1/scale_a), lower.tail = FALSE)
  dev.new()
  par(mar= c(5,6,2,2))
  #Plot both distributions and check manually
  plot(xs, smooth_teo_cum, col = "red", type ="l", lwd = 2, axes = FALSE, xaxs = "i", yaxs = "i", xlim = c(1, (max(data))+5), ann = FALSE, ylim = c(0,1)) 
  points(cumulative_x, emp_cum_y, col = "blue", pch = 20)
  title(ylab = expression(bold("Cumulative Probability, " ~bolditalic(F(LWR)))), cex.lab = 1.5, line = 3)
  title(xlab = expression(bold("Length to Width Ratio, " ~bolditalic(LWR))), cex.lab = 1.5)
  axis(1, lwd = 2, tck =0.03, cex.axis = 1.5, font =2, at = c(1,(seq(0, ((max(data))+10), 5))))
  axis(2, lwd = 2, tck =0.03, cex.axis = 1.5, font =2, at = seq(0,1, by = 0.2), las =1)
  axis(3, lwd = 2, tck =0.03, labels = FALSE, at = c(1,(seq(0, ((max(data))+10), 5))))
  axis(4, lwd = 2, tck =0.03, labels = FALSE, at = seq(0,1, by = 0.2), las =1)
  legend("topright", do.call("expression", list(phrase, "Observed", "Inverse Gamma cdf", phrase2)), lty = c(NA, NA, 1), col = c(NA, "blue", "red"), pch = c(NA, 20, NA), lwd = c(NA, NA, 2),bty = 'o', cex = 1.1, box.lwd = 2, inset = 0.01)
  
  all_d = emp_cum_y - teo_cum
  g = match((min(emp_cum_y - teo_cum)), all_d)
  h = match((max(emp_cum_y - teo_cum)), all_d)
  if (abs(all_d[g]) > abs(all_d[h])) maxx = g else maxx = h
  max_x_d = cumulative_x[maxx]
  
  dev.new()
  par(mar= c(5,6,2,2))
  plot(cumulative_x, all_d, xlab = "LWR",  col ="darkgreen", pch = 20, xlim = c(1, (max(data))+5), ylim =c(-0.2, 0.2), xaxs ="i", yaxs = "i", ann = FALSE, axes = FALSE )
  axis(1, lwd = 2, tck =0.03, cex.axis = 1.5, font =2, at = c(1,(seq(0, ((max(data))+10), 5))))
  axis(3, lwd = 2, tck =0.03, labels = FALSE, at = c(1,(seq(0, ((max(data))+10), 5))))
  axis(2, at = seq( -0.2, 0.2, 0.05 ), lwd = 2, cex.lab = 1.2, las = 1, tck  = 0.03, cex.axis = 1.5)
  axis(4, at = seq(-0.2, 0.2, 0.05) , lwd = 2, labels = FALSE, tck  = 0.03)
  title(ylab = expression(bold("Obs. - Theoret.,  ")~bolditalic(F(LWR)[obs]~ " -" ~F(LWR)[theo])),cex.lab = 1.5, line = 4)
  title(xlab = expression(bold("Length to Width Ratio, " ~bolditalic(LWR))), cex.lab = 1.5)
  points(max_x_d, all_d[maxx], col = "red", pch = 8)
  label = round(max_x_d, 2)
  label = paste("LWR =", label)
  textxy(max_x_d, all_d[maxx],label, cx = 1.3, dcol = "red")
  abline(0,0, lty = 2, lwd = 2, col ="grey" )
  legend("topright", do.call("expression", list(phrase,  phrase2)), lty = c(NA, NA), col = c(NA, NA), pch = c(NA, NA), lwd = c(NA, NA),bty = 'o', cex = 1.1, box.lwd = 2, inset = 0.01)

  # Use monte carlo modelling to test goodness of fit 
  testresult = numeric(2500)
  dval = numeric(2500)
  

  rm(mle)
  rm(params)
  #Set up for loop
  for (i in 1:2500){
    # Generate n random values of LWR using the parameters fit to the observed data. N should be the same as the number of observed LWRS. 
    rand_lwrs = (rpearsonV(n, shape_rho,  location_s,  (1/scale_a)))
    rand_lwrs_cumulative = sort(rand_lwrs)
    # Use MLE to fit an inverse gamma to these randomly generated LWRs
  
    # Optimise the parameter values
    mle = optim(c(1.4,0.001,-0.0001),ll)
    params = mle$par
    shape_f = params [1]
    #Scale is the reciprocal of the actual scale value 
    scale_f = 1/(params[2])
    location_f = params [3]
    parset = list(shape_f, location_f, scale_f)
    teo_cdf = ppearsonV(rand_lwrs_cumulative, params = parset, lower.tail = FALSE)
    rand_cdf = (1:n)/n
    rand_cdf = sort(rand_cdf, decreasing = TRUE)
    rand_D= max(abs(rand_cdf - teo_cdf))
    dval[i] = rand_D
    if (rand_D <= empirical_D){ testresult[i] = 0}
    if (rand_D > empirical_D){ testresult[i] = 1}}

  ### The number of times the Empirical D value is smaller than the Monte carlo D value ###
  sum(testresult)
  Proportion = sum(testresult) / 2500
  ### Proportion of times the D value is smaller than the Monte Carlo D value ###
  dev.new()
  par(mar= c(5,6,2,2))
  yupper = (round(max(dval), 2)) + 0.1
  yupper2 = round(yupper, 1)
  boxplot(dval, range = 0, lwd = 2, lty = 1, axes = FALSE, xaxs ="i", yaxs = "i", ylim = c(0, yupper2) )
  axis(1, tck = 0, labels = FALSE, at = c(0,2), lwd = 2)
  axis(2, lwd = 2, at = seq(0, yupper2, 0.05), las = 1, cex.axis = 1.5, tck = 0.03)
  axis(3, tck = 0, labels = FALSE, at = c(0,2), lwd = 2)
  axis(4, lwd = 2, at = seq(0, yupper2, 0.05), las = 1, cex.axis = 1.5, labels = FALSE, tck = 0.03)
  title(ylab = expression(bold(D - Value)), cex.lab = 1.5, line = 4)
  points(empirical_D, col = "red", pch = "+", cex = 2)
  abline((quantile(dval, 0.9)), 0, lty = 2, col = "grey", lwd = 2)
  abline((quantile(dval, 0.95)), 0, lty = 2, col = "lightgrey", lwd = 2)
  legend("topright", do.call("expression", list(phrase,  phrase2)), lty = c(NA, NA), col = c(NA, NA), pch = c(NA, NA), lwd = c(NA, NA),bty = 'o', cex = 1.1, box.lwd = 2, inset = 0.01)

  # Now try producing bootstapped samples of the real data and fitting to each sample
  startshape = shape_rho
  startscale = scale_a
  startlocation = location_s


  estsa <- sapply(1:1000, function(i) {

  #Sample the data with replacement (the bootstrapping part)
    xi <- sample(data, size = n, replace = TRUE)
    ll = function(par){
    if(par[1]>0 & par[2]>0 & par[3]<min(xi)) return( -sum(log(dinvgamma(xi-    par[3],par[1],par[2]))) )
    else return(Inf)}
    
    optim(c(startshape, startscale, startlocation),ll)})

  bootstrap_shape = numeric(1000)
  bootstrap_location = numeric(1000)
  bootstrap_scale = numeric(1000)
  bootstrap_rollover = numeric(1000)
  for (b in 1:1000){
  bootstrap_shape [b] = estsa[,b]$par[1]
  bootstrap_location[b] = estsa[,b]$par[3]
  bootstrap_scale[b]  = estsa[,b]$par[2]}
  bootstrap_rollover[b] = (bootstrap_location[b]/(bootstrap_shape [b] + 1))+  bootstrap_scale[b]
  bootstrap_scale_recip = 1/bootstrap_scale


  m = matrix(c(bootstrap_shape, bootstrap_location, bootstrap_scale_recip), nrow= 1000, ncol =3)
  ests = apply(m, 1, function(x) dpearsonV(xs, params = x))

  #Plot the first distritbution 
  par(mar= c(5,6,2,2))
  plot(xs, ests[,1], type = "l", col=rgb(.6, .6, .6, .1), ylim= c(10^-8, 10), log = "xy", ylab = "", xlab=expression(bold("Length to Width Ratio, "~bolditalic(LWR))), axes = FALSE, yaxs ="i", xaxs = "i", cex.lab = 1.7 )
  ticks = c(0.1,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)
  title(ylab = expression(bold("Probability Density, " ~bolditalic(p(LWR)))), cex.lab = 1.7, line = 4)
  axis(1, at = c(1,10, 100), labels =( expression(~bold(10^0),~bold(10^1), ~bold(10^2))), lty = 1, lwd = 0, tck = 0, cex.axis = 1.7)
  yticks = c(0.000000001, 0.000000002, 0.000000003, 0.000000004, 0.000000005,  0.000000006, 0.000000007, 0.000000008, 0.000000009, 0.00000001, 0.00000002, 0.00000003, 0.00000004, 0.00000005,  0.00000006, 0.00000007, 0.00000008, 0.00000009,0.0000001, 0.0000002, 0.0000003, 0.0000004, 0.0000005,  0.0000006, 0.0000007, 0.0000008, 0.0000009, 0.000001, 0.000002, 0.000003,  0.000004, 0.000005, 0.000006, 0.000007, 0.000008, 0.000009, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.00008, 0.00009,0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10)
  axis(2, at = c(0.000000001, 0.00000001,0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1,10), labels =(expression(~bold(10^-9),~bold(10^-8),~bold(10^-7),~bold(10^-6),~bold(10^-5), ~bold(10^-4), ~bold(10^-3), ~bold(10^-2), ~bold(10^-1), ~bold(10^0) , ~bold(10^1))), lty = 1, lwd = 0, tck = 0, cex.axis = 1.7, las = 1)

   for(i in 2:500)lines(xs, ests[, i], col=rgb(.6, .6, .6, .1))
  quants <- apply(ests, 1, quantile, c(0.025, 0.5, 0.975))
  median_pdf = quants[2,]
  lower_pdf = quants[3,]
  upper_pdf = quants[1,]
  lines(xs, quants[1, ], col="red", lwd=1.5, lty=2)
  lines(xs, quants[3, ], col="red", lwd=1.5, lty=2)
  lines(xs, quants[2, ], col="darkred", lwd=2)
  axis(1, ticks, labels = NA, lty = 1, lwd = 2, tck = 0.01)
  axis(2, yticks, labels = NA, lty = 1, lwd = 2, tck = 0.01)
  axis(3, ticks, labels = NA, lty = 1, lwd = 2, tck = 0.01)
  axis(4, yticks, labels = NA, lty = 1, lwd = 2, tck = 0.01)
  legend("topright", do.call("expression", list(phrase, "Fit to bootstrapped samples", "Median fit", "5th/95th Percentile Fit", phrase2)), lty = c(NA, 1,1, 2, NA), col = c(NA, "Grey", "red", "red", NA),  lwd = c(NA, 0.5,2,1, NA),bty = 'o', cex = 1.1, box.lwd = 2, inset = 0.01)

  bootstrap_rollover = numeric(1000)
  for(i in 1:1000){
    max_pdf = max(ests[,i])
    loc_max = match(max_pdf, ests[,i])
    location_1 = xs[loc_max]
    bootstrap_rollover[i] = location_1}

  variable1 = paste("shape", bin, sep = "")
  assign(variable1, bootstrap_shape)
  variable2 = paste("location", bin, sep = "")
  assign(variable2, bootstrap_location)
  variable3 = paste("scale", bin, sep = "")
  assign(variable3, bootstrap_scale)
  variable4 = paste("rollover", bin, sep = "")
  assign(variable4, bootstrap_rollover)
  variable5 = paste("d", bin, sep = "")
  assign(variable5, empirical_D)
  variable6 = paste("dprop", bin, sep = "")
  assign(variable6, Proportion)
  variable7 = paste("median_pdf", bin, sep = "")
  assign(variable7, median_pdf)
  variable8 = paste("upper_pdf", bin, sep = "")
  assign(variable8, upper_pdf)
  variable9 = paste("lower_pdf", bin, sep = "")
  assign(variable9, lower_pdf)
  
  print(bin)
}


##################
# dev.new()
# par(mar= c(5,6,2,2))
# plot(xs, median_pdfTaiwanLW3, type = "l", col = "darkred", lwd = 3, ylim= c(10^-8, 10), log = "xy", ylab = "", xlab=expression(bold("Length to Width Ratio, "~bolditalic(LWR))), axes = FALSE, yaxs ="i", xaxs = "i", cex.lab = 1.7 )
# ticks = c(0.1,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100)
# title(ylab = expression(bold("Probability Density, " ~bolditalic(p(LWR)))), cex.lab = 1.7, line = 4)
# axis(1, at = c(1,10, 100), labels =( expression(~bold(10^0),~bold(10^1), ~bold(10^2))), lty = 1, lwd = 0, tck = 0, cex.axis = 1.7)
# yticks = c(0.000000001, 0.000000002, 0.000000003, 0.000000004, 0.000000005,  0.000000006, 0.000000007, 0.000000008, 0.000000009, 0.00000001, 0.00000002, 0.00000003, 0.00000004, 0.00000005,  0.00000006, 0.00000007, 0.00000008, 0.00000009,0.0000001, 0.0000002, 0.0000003, 0.0000004, 0.0000005,  0.0000006, 0.0000007, 0.0000008, 0.0000009, 0.000001, 0.000002, 0.000003,  0.000004, 0.000005, 0.000006, 0.000007, 0.000008, 0.000009, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.00008, 0.00009,0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10)
# axis(2, at = c(0.000000001, 0.00000001,0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1,10), labels =(expression(~bold(10^-9),~bold(10^-8),~bold(10^-7),~bold(10^-6),~bold(10^-5), ~bold(10^-4), ~bold(10^-3), ~bold(10^-2), ~bold(10^-1), ~bold(10^0) , ~bold(10^1))), lty = 1, lwd = 0, tck = 0, cex.axis = 1.7, las = 1)
# #polygon(c(lwr, rev(lwr)), c(lower_pdfTaiwanLW1, rev(upper_pdfTaiwanLW1)), col = rgb(1,0,0,0.3), border = NA)
# #polygon(c(lwr, rev(lwr)), c(lower_pdfTaiwanLW1, rev(upper_pdfTaiwanLW2)), col = rgb(1,0.5,0,0.3), border = NA)
# #polygon(c(lwr, rev(lwr)), c(lower_pdfTaiwanLW1, rev(upper_pdfTaiwanLW3)), col = rgb(1,1,0.2,0.3), border = NA)
# #polygon(c(lwr, rev(lwr)), c(lower_pdfTaiwanLW1, rev(upper_pdfTaiwanLW4)), col = rgb(0.5,1,0.2,0.3), border = NA)
# #polygon(c(lwr, rev(lwr)), c(lower_pdfTaiwanLW1, rev(upper_pdfTaiwanLW5)), col = rgb(0.2,0.5,1,0.3), border = NA)
# #polygon(c(lwr, rev(lwr)), c(lower_pdfTaiwanLW1, rev(upper_pdfTaiwanLW6)), col = rgb(0.5,0.5,1,0.3), border = NA)
# #polygon(c(lwr, rev(lwr)), c(lower_pdfTaiwanLW1, rev(upper_pdfTaiwanLW7)), col = rgb(0.5,0  ,1  ,0.3), border = NA)
# #polygon(c(lwr, rev(lwr)), c(lower_pdfTaiwanLW1, rev(upper_pdfTaiwanLW8)), col = rgb(1,1,1 ,0.3), border = NA)
# #polygon(c(lwr, rev(lwr)), c(lower_pdfTaiwanLW1, rev(upper_pdfTaiwanLW9)), col = rgb(1,0.5,1 ,0.3), border = NA)
# #polygon(c(lwr, rev(lwr)), c(lower_pdfTaiwanLW1, rev(upper_pdfTaiwanLW10)), col = rgb(0.3,0,0.5,0.3), border = NA)
# 
# #points(xs, median_pdfTaiwanLW1, col = "darkred", type = "l", lwd = 3)
# #points(xs, median_pdfTaiwanLW2, col = "red", type = "l", lwd = 3)
# points(xs, median_pdfTaiwanLW3, col = "orange", type = "l", lwd = 3)
# points(xs, median_pdfTaiwanLW4, col = "gold", type = "l", lwd = 3)
# points(xs, median_pdfTaiwanLW5, col = "yellow", type = "l", lwd = 3)
# points(xs, median_pdfTaiwanLW6, col = "lawngreen", type = "l", lwd = 3)
# points(xs, median_pdfTaiwanLW7, col = "darkgreen", type = "l", lwd = 3)
# points(xs, median_pdfTaiwanLW8, col = "blue", type = "l", lwd = 3)
# points(xs, median_pdfTaiwanLW9, col = "darkblue", type = "l", lwd = 3)
# points(xs, median_pdfTaiwanLW10, col = "magenta", type = "l", lwd = 3)
# axis(1, ticks, labels = NA, lty = 1, lwd = 2, tck = 0.01)
# axis(2, yticks, labels = NA, lty = 1, lwd = 2, tck = 0.01)
# axis(3, ticks, labels = NA, lty = 1, lwd = 2, tck = 0.01)
# axis(4, yticks, labels = NA, lty = 1, lwd = 2, tck = 0.01)
# legend("topright", do.call("expression", list(phrase_1, phrase_2, phrase_3, phrase_4, phrase_5, phrase_6, phrase_7, phrase_8, phrase_9, phrase_10)), lty = 1, col = c("darkred", "red", "Orange", "gold", "yellow", "lawngreen", "darkgreen", "Blue", "darkblue", "Magenta"),  lwd = 2 ,bty = 'o', cex = 0.8, box.lwd = 2, inset = 0.01)
# 
# # 
# dev.new()
# par(mar= c(8,6,2,2))
# boxplot(rolloverTaiwanLW1, rolloverTaiwanLW2, rolloverTaiwanLW3, rolloverTaiwanLW4, rolloverTaiwanLW5, rolloverTaiwanLW6, rolloverTaiwanLW7, rolloverTaiwanLW8, rolloverTaiwanLW9, rolloverTaiwanLW10, lty =1, lwd =1.5, range = 0, pch =1, axes = FALSE, xaxs = "i", yaxs = "i", ylim = c(1, 5))
# #boxplot(rolloverTaiwanLW3, rolloverTaiwanLW4, rolloverTaiwanLW5, rolloverTaiwanLW6, rolloverTaiwanLW7, rolloverTaiwanLW8, rolloverTaiwanLW9, rolloverTaiwanLW10, lty =1, lwd =1.5, range = 0, pch =1, axes = FALSE, xlim = c(0, 11), xaxs = "i", yaxs = "i", ylim = c(1, 5), at  = c( 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7))
# box(lwd = 2)
# axis(1, at = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7), labels =( expression(~bold("0"), ~bold("100"), ~bold("200"), ~bold("400"), ~bold("800"), ~bold("1,600"), ~bold("3,200"), ~bold("6,400"), ~bold("12,800"), ~bold("25,600"), ~bold("250,000"))), lty = 1, lwd = 2, tck = 0.03, cex.axis = 1.3, las =2)
# box(lwd = 2)  
# axis(2, lty =1, lwd = 2, tck = 0.03, cex.axis = 1.5, font = 2, las =1)
# axis(3, at = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7), tck = 0.03, labels = FALSE, lwd = 2)
# axis(4, lty = 1, lwd = 2, tck = 0.02, labels = FALSE)
# title(ylab = expression(bold ("Location of Max. Probability, " ~  bolditalic((LWR)))), line = 3, cex.lab = 1.5)
# title(xlab=expression(bold ("Landslide Area Bin, "~bolditalic(A[L]) ~(m^2))) , line = 5, cex.lab = 1.5)
# 
# 
# dev.new()
# par(mar= c(8,6,2,2))
# boxplot(locationTaiwanLW1, locationTaiwanLW2, locationTaiwanLW3, locationTaiwanLW4, locationTaiwanLW5, locationTaiwanLW6, locationTaiwanLW7, locationTaiwanLW8, locationTaiwanLW9, locationTaiwanLW10, lty =1, lwd =1.5, range = 0, pch =1, axes = FALSE, xaxs = "i", yaxs = "i")
# #boxplot(locationTaiwanLW3, locationTaiwanLW4, locationTaiwanLW5, locationTaiwanLW6, locationTaiwanLW7, locationTaiwanLW8, locationTaiwanLW9, locationTaiwanLW10, lty =1, lwd =1.5, range = 0, pch =1, axes = FALSE, xlim = c(0, 11), xaxs = "i", yaxs = "i",  at  = c( 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7))
# box(lwd = 2)
# axis(1, at = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7), labels =( expression(~bold("0"), ~bold("100"), ~bold("200"), ~bold("400"), ~bold("800"), ~bold("1,600"), ~bold("3,200"), ~bold("6,400"), ~bold("12,800"), ~bold("25,600"), ~bold("250,000"))), lty = 1, lwd = 2, tck = 0.03, cex.axis = 1.3, las =2)
# box(lwd = 2)  
# axis(2, lty =1, lwd = 2, tck = 0.03, cex.axis = 1.5, font = 2, las =1)
# axis(3, at = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7), tck = 0.03, labels = FALSE, lwd = 2)
# axis(4, lty = 1, lwd = 2, tck = 0.02, labels = FALSE)
# title(ylab = expression(bold ("Location Parameter Value, " ~bolditalic(s))), line = 4, cex.lab = 1.5)
# title(xlab=expression(bold ("Landslide Area Bin, "~bolditalic(A[L]) ~(m^2))) , line = 5, cex.lab = 1.5)
# 
# dev.new()
# par(mar= c(8,6,2,2))
# #boxplot(scaleTaiwanLW1, scaleTaiwanLW2, scaleTaiwanLW3, scaleTaiwanLW4, scaleTaiwanLW5, scaleTaiwanLW6, scaleTaiwanLW7, scaleTaiwanLW8, scaleTaiwanLW9, scaleTaiwanLW10, lty =1, lwd =1.5, range = 0, pch =1, axes = FALSE, xaxs = "i", yaxs = "i", ylim = c(0,40))
# boxplot(scaleTaiwanLW3, scaleTaiwanLW4, scaleTaiwanLW5, scaleTaiwanLW6, scaleTaiwanLW7, scaleTaiwanLW8, scaleTaiwanLW9, scaleTaiwanLW10, lty =1, lwd =1.5, range = 0, pch =1, axes = FALSE, ylim = c(0,50), xlim = c(0, 11), xaxs = "i", yaxs = "i",  at  = c( 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7))
# box(lwd = 2)
# axis(1, at = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7), labels =( expression(~bold("0"), ~bold("100"), ~bold("200"), ~bold("400"), ~bold("800"), ~bold("1,600"), ~bold("3,200"), ~bold("6,400"), ~bold("12,800"), ~bold("25,600"), ~bold("250,000"))), lty = 1, lwd = 2, tck = 0.03, cex.axis = 1.3, las =2, add= TRUE)
# box(lwd = 2)  
# axis(2, lty =1, lwd = 2, tck = 0.03, cex.axis = 1.5, font = 2, las =1)
# axis(3, at = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7), tck = 0.03, labels = FALSE, lwd = 2)
# axis(4, lty = 1, lwd = 2, tck = 0.02, labels = FALSE)
# title(ylab = expression(bold ("Scale Parameter Value, " ~bolditalic(a))), line = 4, cex.lab = 1.5)
# title(xlab=expression(bold ("Landslide Area Bin, "~bolditalic(A[L]) ~(m^2))) , line = 6, cex.lab = 1.5)
# 
# dev.new()
# par(mar= c(8,6,2,2))
# #boxplot(shapeTaiwanLW1, shapeTaiwanLW2, shapeTaiwanLW3, shapeTaiwanLW4, shapeTaiwanLW5, shapeTaiwanLW6, shapeTaiwanLW7, shapeTaiwanLW8, shapeTaiwanLW9, shapeTaiwanLW10, lty =1, lwd =1.5, range = 0, pch =1, axes = FALSE, xaxs = "i", yaxs = "i", ylim = c(0,15))
# boxplot(shapeTaiwanLW3, shapeTaiwanLW4, shapeTaiwanLW5, shapeTaiwanLW6, shapeTaiwanLW7, shapeTaiwanLW8, shapeTaiwanLW9, shapeTaiwanLW10, lty =1, lwd =1.5, range = 0, pch =1, axes = FALSE, ylim = c(0,20), xlim = c(0, 11), xaxs = "i", yaxs = "i",  at  = c( 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7))
# box(lwd = 2)
# axis(1, at = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7), labels =( expression(~bold("0"), ~bold("100"), ~bold("200"), ~bold("400"), ~bold("800"), ~bold("1,600"), ~bold("3,200"), ~bold("6,400"), ~bold("12,800"), ~bold("25,600"), ~bold("250,000"))), lty = 1, lwd = 2, tck = 0.03, cex.axis = 1.3, las =2)
# box(lwd = 2)  
# axis(2, lty =1, lwd = 2, tck = 0.03, cex.axis = 1.5, font = 2, las =1)
# axis(3, at = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7), tck = 0.03, labels = FALSE, lwd = 2)
# axis(4, lty = 1, lwd = 2, tck = 0.02, labels = FALSE)
# title(ylab = expression(bold ("shape Parameter Value, " ~bolditalic(rho))), line = 4, cex.lab = 1.5)
# title(xlab=expression(bold ("Landslide Area Bin, "~bolditalic(A[L]) ~(m^2))) , line = 6, cex.lab = 1.5)
# 
# 
# Taiwan_median_shape = numeric(10)
# Taiwan_median_location = numeric(10)
# Taiwan_median_scale = numeric(10)
# Taiwan_median_rollover = numeric(10)
# for (j in 1:10){
#   bina = paste("shapeTaiwan","LW", j, sep = "")
#   binb = paste("locationTaiwan","LW", j, sep = "")
#   binc = paste("scaleTaiwan","LW", j, sep = "")
#   bind = paste("rolloverTaiwan","LW", j, sep = "")
#   Taiwan_median_shape[j] = median(get(bina))
#   Taiwan_median_location[j] = median(get(binb))
#   Taiwan_median_scale[j] = median(get(binc))
#   Taiwan_median_rollover[j] = median(get(bind))}
# 
# Taiwan_q25_shape = numeric(10)
# Taiwan_q25_location = numeric(10)
# Taiwan_q25_scale = numeric(10)
# Taiwan_q25_rollover = numeric(10)
# for (j in 1:10){
#   bina = paste("shapeTaiwan","LW", j, sep = "")
#   binb = paste("locationTaiwan","LW", j, sep = "")
#   binc = paste("scaleTaiwan","LW", j, sep = "")
#   bind = paste("rolloverTaiwan","LW", j, sep = "")
#   Taiwan_q25_shape[j] = quantile(get(bina), 0.25)
#   Taiwan_q25_location[j] = quantile(get(binb), 0.25)
#   Taiwan_q25_scale[j] = quantile(get(binc), 0.25)
#   Taiwan_q25_rollover[j] = quantile(get(bind), 0.25)}
# 
# Taiwan_q75_shape = numeric(10)
# Taiwan_q75_location = numeric(10)
# Taiwan_q75_scale = numeric(10)
# Taiwan_q75_rollover = numeric(10)
# for (j in 1:10){
#   bina = paste("shapeTaiwan","LW", j, sep = "")
#   binb = paste("locationTaiwan","LW", j, sep = "")
#   binc = paste("scaleTaiwan","LW", j, sep = "")
#   bind = paste("rolloverTaiwan","LW", j, sep = "")
#   Taiwan_q75_shape[j] = quantile(get(bina), 0.75)
#   Taiwan_q75_location[j] = quantile(get(binb), 0.75)
#   Taiwan_q75_scale[j] = quantile(get(binc), 0.75)
#   Taiwan_q75_rollover[j] = quantile(get(bind), 0.75)}
# 
# dev.new()
# par(mar= c(8,6,2,2))
# #boxplot(shapeTaiwanLW1, shapeTaiwanLW2, shapeTaiwanLW3, shapeTaiwanLW4, shapeTaiwanLW5, shapeTaiwanLW6, shapeTaiwanLW7, shapeTaiwanLW8, shapeTaiwanLW9, shapeTaiwanLW10, lty =1, lwd =1.5, range = 0, pch =1, axes = FALSE, xaxs = "i", yaxs = "i", ylim = c(0,15))
# boxplot(shapeTaiwanLW3, shapeTaiwanLW4, shapeTaiwanLW5, shapeTaiwanLW6, shapeTaiwanLW7, shapeTaiwanLW8, shapeTaiwanLW9, shapeTaiwanLW10, lty =1, lwd =1.5, range = 0, pch =1, axes = FALSE, ylim = c(0,20), xlim = c(0, 11), xaxs = "i", yaxs = "i",  at  = c( 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7))
# box(lwd = 2)
# axis(1, at = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7), labels =( expression(~bold("0"), ~bold("100"), ~bold("200"), ~bold("400"), ~bold("800"), ~bold("1,600"), ~bold("3,200"), ~bold("6,400"), ~bold("12,800"), ~bold("25,600"), ~bold("250,000"))), lty = 1, lwd = 2, tck = 0.03, cex.axis = 1.3, las =2)
# box(lwd = 2)  
# axis(2, lty =1, lwd = 2, tck = 0.03, cex.axis = 1.5, font = 2, las =1)
# axis(3, at = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7), tck = 0.03, labels = FALSE, lwd = 2)
# axis(4, lty = 1, lwd = 2, tck = 0.02, labels = FALSE)
# title(ylab = expression(bold ("shape Parameter Value, " ~bolditalic(rho))), line = 4, cex.lab = 1.5)
# title(xlab=expression(bold ("Landslide Area Bin, "~bolditalic(A[L]) ~(m^2))) , line = 6, cex.lab = 1.5)
# 
# 
# errbar(1:10, Northridge_median_shape, Northridge_q25_shape, Northridge_q75_shape, col = "aquamarine4", pch = 0, ylim = c(1,7),xaxt='n', yaxt='n',xlab=NA,ylab=NA,xaxs = "i", yaxs = "i", xlim = c(0, 11), cex = 1.3, lwd = 4 )
# par(new = TRUE)
# errbar(1:10, Guatemala_median_shape, Guatemala_q25_shape, Guatemala_q75_shape, col = "orange", pch = 2,ylim = c(1,7),xaxt='n', yaxt='n',xlab=NA,ylab=NA,xaxs = "i", yaxs = "i", xlim = c(0, 11), cex = 1.3, lwd = 4 )
# par(new = TRUE)
# errbar(1:10, Umbria_median_shape, Umbria_q25_shape, Umbria_q75_shape, col = "DarkOrchid", pch = 4,ylim = c(1,7),xaxt='n', yaxt='n',xlab=NA,ylab=NA,xaxs = "i", yaxs = "i", xlim = c(0, 11), cex = 1.3, lwd = 4 )
# axis(1, at = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7), labels =( expression(~bold("0"), ~bold("100"), ~bold("200"), ~bold("400"), ~bold("800"), ~bold("1,600"), ~bold("3,200"), ~bold("6,400"), ~bold("12,800"), ~bold("25,600"), ~bold("250,000"))), lty = 1, lwd = 2, tck = 0.03, cex.axis = 1.3, las =2)
# box(lwd = 2)  
# axis(2, lty =1, lwd = 2, tck = 0.03, cex.axis = 1.5, font = 2, las =1)
# axis(3, at = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.7), tck = 0.03, labels = FALSE, lwd = 2)
# axis(4, lty = 1, lwd = 2, tck = 0.02, labels = FALSE)
# title(ylab = expression(bold ("shape Parameter Value, " ~bolditalic(rho))), line = 4, cex.lab = 1.5)
# title(xlab=expression(bold ("Landslide Area Bin, "~bolditalic(A[L]) ~(m^2))) , line = 6, cex.lab = 1.5)
# legend("topright", c("Northridge", "Guatemala", "Umbria"), col = c("aquamarine4", "orange", "DarkOrchid"), pch = c(0,2,4),bty = 'o', cex = 1.1, box.lwd = 2, inset = 0.01)