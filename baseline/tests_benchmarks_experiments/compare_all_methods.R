library(ggplot2)
library(hyperSpec)

source('baseline/baseline.Goldindec.R')
source('baseline/baseline.wavelet.R')
source('baseline/baseline.mor.R')

load("data/testspc.RData")

# VISUAL COMPARATION ===============================
s <- as.hyperSpec(testspc[75,,drop=F])
s <- s[,,500~3500]
plotspc(s)


bl.Goldindec <- baseline.Goldindec(s, 4, 0.5)          # Goldinec
bl.wavelet   <- baseline.wavelet(s)                  # Wavelet
bl.mor       <- baseline.mor(s, mor_getwopt(s$spc[1,]))   # simple MOR
bl.imor      <- baseline.imor(s)                         #  I-MOR
bl.polynom4  <- hyperSpec::spc.fit.poly.below(s, poly.order = 4) # 4 order polynom
bl.polynom5  <- hyperSpec::spc.fit.poly.below(s, poly.order = 5) # 5 order polynom

# IMOR+Smoothing
s_smoothed <- prospectr::savitzkyGolay(s$spc, m = 0, p=2, w=9)
plotspc(as.hyperSpec(s_smoothed))
bl.imor_smooth <- s
bl.imor_smooth$spc[,colnames(bl.imor_smooth$spc) %in% colnames(s_smoothed)] <- baseline.imor(s_smoothed)
plotspc(bl.imor_smooth,add=T,col = 'red')


plotdata <- rbind(s, bl.Goldindec, bl.wavelet, bl.mor, bl.imor, bl.imor_smooth, bl.polynom4, bl.polynom5)
plotdata$method = factor( c('Original spectrum', 'Goldindec', 'Wavelet', 'MOR', 'I-MOR', 'I-MOR + Smooth', '4th order polynom', '5th order polynom'), 
                  levels=c('Original spectrum', 'Goldindec', 'Wavelet', 'MOR', 'I-MOR', 'I-MOR + Smooth', '4th order polynom', '5th order polynom'))

n <- nrow(plotdata)
hues = seq(15, 375, length = n)
hcl <- hcl(h = hues, l = 65, c = 100)[1:(n-1)]

qplotspc(plotdata, spc.nmax = 1000, mapping = aes(x = .wavelength, y = spc, colour = method, group = .rownames)) + 
  scale_color_manual(values=c('black',hcl)) +
  scale_alpha_manual(values=c(1,rep(0.75,n-1))) +
  scale_x_continuous(breaks = seq(500,3500,by=200), minor_breaks = seq(600,3400,by=100)) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=16),
    legend.position = c(0.5,0.1),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size=14),
    legend.background = element_blank(),
    legend.key.size = unit(1,"cm")
    ) +
  guides(color = guide_legend(override.aes = list(size = 2))) 
  
  
qplotspc(plotdata[,,1600~2300], spc.nmax = 1000, mapping = aes(x = .wavelength, y = spc, colour = method, group = .rownames)) + 
  scale_color_manual(values=c('black',hcl)) +
  scale_alpha_manual(values=c(1,rep(0.75,n-1))) +
  scale_x_continuous(breaks = seq(1600,2300,by=100), minor_breaks = seq(1650,2300,by=50)) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=16),
    legend.position = c(0.5,0.1),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size=14),
    legend.background = element_blank(),
    legend.key.size = unit(1,"cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2)))

## BENCHMARKS =================================
library(microbenchmark)
# Start from 10 spectra
s <- testspc[1:10, ]
res <- microbenchmark(
  'Goldinec' = baseline.Goldindec(s, NULL, 4, 0.5),
  'Wavelet'  = baseline.wavelet(s),
  'MOR'      = baseline.mor(s, 30),
  'I-MOR'    = baseline.imor(s),
  'Polynom'  = hyperSpec::spc.fit.poly.below(as.hyperSpec(s), poly.order = 4),
  times=20
)

print(res,'s')
# Unit: seconds
# expr         min          lq         mean       median           uq         max      neval cld
# Goldinec 41.61330384 41.70685780 41.949858492 41.783434968 42.234039327 42.72612477    20   c
# Wavelet   5.05897899  5.15703974  5.208787424  5.187343574  5.218990699  5.65962530    20   b 
# MOR       0.00735328  0.00750078  0.008395001  0.007592707  0.009285892  0.01345498    20   a  
# I-MOR     0.08289533  0.08412170  0.086789738  0.086890019  0.088734400  0.09280998    20   a  
# Polynom   0.01710310  0.01904076  0.022106710  0.019728261  0.020335955  0.05906408    20   a  

tmp_res <- res
tmp_res$time <- res$time / 1000000000
y_min <- 0
y_max = 1.05 * max(tmp_res$time)
plt <- ggplot2::ggplot(tmp_res, ggplot2::aes_string(x = "expr", y = "time"))
plt <- plt + ggplot2::coord_cartesian(ylim = c(y_min, y_max))
plt <- plt + ggplot2::geom_violin(width=8)
plt <- plt + ggplot2::scale_x_discrete(name = "")
plt <- plt + ggplot2::scale_y_log10(name = 'Time [seconds]', breaks = c(0,0.01,0.05,0.1,1,5,10,50,100))
plt <- plt + ggplot2::coord_flip()
plt + theme(axis.text = element_text(size=14))
