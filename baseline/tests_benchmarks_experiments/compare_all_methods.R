library(ggplot2)
library(hyperSpec)

source('baseline/Goldinec.R')
source('baseline/baseline.wavelet.R')
# source('baseline/mor.R')
Rcpp::sourceCpp('baseline/mor.cpp')

load("data/testspc.RData")

wl <- as.numeric(colnames(testspc))
test_spc <- testspc[75, wl>499 & wl<3501,drop=F]
wl <- as.numeric(colnames(test_spc))
plotspc(as.hyperSpec(test_spc))

goldindecbl <- goldindec(test_spc, 4, 0.5)                  # Goldinec
waveletbl   <- baseline.wavelet(test_spc)                   # Wavelet
morbl       <- baseline_mor(test_spc, mor_getwopt(test_spc))   # simple MOR
imorbl      <- baseline_imor(test_spc)                         #  I-MOR
polybl4     <- spc.fit.poly.below(as.hyperSpec(test_spc), poly.order = 4) # 4 order polynom
polybl5     <- spc.fit.poly.below(as.hyperSpec(test_spc), poly.order = 5) # 5 order polynom

test_spc_smoothed <- prospectr::savitzkyGolay(test_spc, m = 0, p=2, w=9)
plotspc(as.hyperSpec(test_spc_smoothed))
imorsmoothbl <-test_spc 
imorsmoothbl[colnames(imorsmoothbl) %in% colnames(test_spc_smoothed)] <- baseline_imor(test_spc_smoothed)

plotdata <- as.hyperSpec(
  rbind(
    test_spc,
    goldindecbl,
    t(waveletbl),
    matrix(morbl,nrow=1),
    matrix(imorbl,nrow=1),
    matrix(imorsmoothbl,nrow=1),
    polybl4$spc,
    polybl5$spc
  ),
  wl = wl,
  data = data.frame(method = factor(       c('Original spectrum', 'Goldindec', 'Wavelet', 'MOR', 'I-MOR', 'I-MOR + Smooth', '4th order polynom', '5th order polynom'), 
                                    levels=c('Original spectrum', 'Goldindec', 'Wavelet', 'MOR', 'I-MOR', 'I-MOR + Smooth', '4th order polynom', '5th order polynom')))
)

n <- nrow(plotdata)
hues = seq(15, 375, length = n)
hcl <- hcl(h = hues, l = 65, c = 100)[1:(n-1)]

qplotspc(plotdata, spc.nmax = 1000, mapping = aes(x = .wavelength, y = spc, colour = method, group = .rownames)) + 
  scale_color_manual(values=c('black',hcl)) +
  scale_alpha_manual(values=c(1,rep(0.75,n-1))) +
  scale_x_continuous(breaks = seq(500,3500,by=200), minor_breaks = seq(600,3400,by=100)) +
  theme(
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