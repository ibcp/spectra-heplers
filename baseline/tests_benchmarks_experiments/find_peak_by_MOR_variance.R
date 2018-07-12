library(ggplot2)
library(hyperSpec)

Rcpp::sourceCpp('baseline/baseline.mor.cpp')
source('baseline/baseline.mor.R')

load("data/testspc.RData")

wl <- as.numeric(colnames(testspc))
test_spc <- testspc[75, wl>499 & wl<3501,drop=F]
wl <- as.numeric(colnames(test_spc))
plotspc(as.hyperSpec(test_spc))

test_spc_smoothed <- prospectr::savitzkyGolay(test_spc, m = 0, p=2, w=15)
plotspc(as.hyperSpec(test_spc_smoothed))
test_spc <- test_spc_smoothed
wl <- as.numeric(colnames(test_spc_smoothed))

# I-MOR BASELINE
imorbl <- baseline_imor(test_spc)

# GET MOR BASELINES FROM 1 to wopt+10
wopt <- mor_getwopt(test_spc)
wopt <- 50

bl <- sapply(1:wopt, function(w) {baseline.mor(test_spc, w)})
bl <- t(bl)

# PLOT ALL MOR BASELINES
plotspc(as.hyperSpec(test_spc))
for (i in 1:wopt) {
  plotspc(as.hyperSpec(bl,wl=wl), col='red',add=T)
}


# PLOT I-MOR, MEAN(MOR), MEDIAN(MOR) BASELINES
mormedian <- apply(bl, 2, median)
plotspc(as.hyperSpec(test_spc))
plotspc(as.hyperSpec(matrix(mormedian,nrow = 1), wl=wl), col = 'red', add=T)
plotspc(as.hyperSpec(matrix(colMeans(bl),nrow = 1), wl=wl), col = 'green', add=T)
plotspc(as.hyperSpec(matrix(imorbl,nrow = 1), wl=wl), col = 'blue', add=T)


# PLOT MEAN(MOR) PM SD
morsd <- apply(bl, 2, sd)
qplot(morsd,geom = 'histogram', bins=51)

mean.pm.sd <- aggregate (as.hyperSpec(bl,wl=wl), factor(rep('all',nrow(bl))), mean_pm_sd)
plot (mean.pm.sd, col = matlab.palette (1), fill = ".aggregate")
plotspc(as.hyperSpec(test_spc), add=T, col=alpha('black',0.4))
points(x=wl[morsd<1], y=test_spc[1,morsd<1],pch=16)

# USE LOW SD REGIONS AS ANCHOR POINTS
plot (mean.pm.sd[,,min~1000], col = matlab.palette (1), fill = ".aggregate")
plotspc(as.hyperSpec(test_spc), add=T, col=alpha('black',0.4))
# points(x=wl[c(18,24)], y=test_spc[1,c(18,24)],pch=16,col='darkgreen')

indx <- which(morsd<10)
x <- wl[indx]
y <- test_spc[1,indx]

diff_left  <- diff(y)[1:(length(indx)-2)]
diff_right <- diff(y)[2:(length(indx)-1)]

indx_ancors <- which(diff_left<0 & diff_right>0)+1 
x_ancors <- x[indx_ancors]
y_ancors <- y[indx_ancors]
points(x=x, y=y, pch=16)
points(x=x_ancors, y=y_ancors, pch=16, col='darkgreen')

blspline_fmm      <- spline(x = x_ancors, y = y_ancors, xout = wl[wl>=min(x_ancors) & wl<=max(x_ancors)], method = "fmm")
blspline_natural  <- spline(x = x_ancors, y = y_ancors, xout = wl[wl>=min(x_ancors) & wl<=max(x_ancors)], method = "natural")

plotspc(as.hyperSpec(test_spc))
lines(x=blspline_fmm$x, y=blspline_fmm$y, pch=16 ,col='red')
lines(x=blspline_natural$x, y=blspline_natural$y, pch=16 ,col='green')

# SMOOTH I-MOR
imor_blspline_fmm      <- spline(x = as.numeric(colnames(test_spc)), y = imorbl, xout = as.numeric(colnames(test_spc)), method = "fmm")
imor_blspline_natural  <- spline(x = as.numeric(colnames(test_spc)), y = imorbl, xout = as.numeric(colnames(test_spc)), method = "natural")

plotspc(as.hyperSpec(test_spc))
lines(x=imor_blspline_fmm$x,     y=imor_blspline_fmm$y,    pch=16 ,col='red')
lines(x=imor_blspline_natural$x, y=imor_blspline_natural$y, pch=16 ,col='green')
imorbl_rollmedian <- zoo::rollmedian(x = imorbl, k=51)
lines(x=as.numeric(names(imorbl_rollmedian)), y=imorbl_rollmedian, pch=15, col='blue')


# DENOISING BY FFT
Y1 <- fft(test_spc)
plot(abs(Y1), type="h", ylim = c(0,1))
for (i in 1:length(Y1)) {
  if (abs(Y1[i]) <1) {
    Y1[i] <- 0
  }
}
plot(Arg(Y1), type="h")


# t <- seq(0,1,1/1000)
# y1 <- sin(2*pi*100*t) + 0.5*cos(2*pi*200*t) + 0.25*sin(2*pi*250*t + pi/3)
# Y1 <- fft(y1)
# plot(abs(Y1), type="h")
# plot(Arg(Y1), type="h")


t <- seq(0, 1, len = 100) # 1 second sample
x <- sin(2*pi*t*2.3) + 0.25*rnorm(length(t)) # 2.3 Hz sinusoid+noise
z <- fftfilt(rep(1, 10)/10, x) # apply 10-point averaging filter
plot(t, x, type = "l")
lines(t, z, col = "red")




bf <- butter(3, 0.1) # 10 Hz low-pass filter
t <- seq(0, 1, len = 100) # 1 second sample
x <- sin(2*pi*t*2.3) + 0.25*rnorm(length(t)) # 2.3 Hz sinusoid+noise
z <- filter(bf, x) # apply filter
plot(t, x, type = "l")
lines(t, z, col = "red")

wl <- as.numeric(colnames(testspc))
test_spc <- testspc[75, wl>499 & wl<3501,drop=F]
wl <- as.numeric(colnames(test_spc))
plotspc(as.hyperSpec(test_spc))
z <- filter(bf, test_spc[1,]) # apply filter
#plot(x = as.numeric(colnames(test_spc)), y=test_spc[1,], type='l')
lines(x = as.numeric(colnames(test_spc)), z, col = "red")