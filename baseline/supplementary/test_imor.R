library(hyperSpec)
library(dplyr)
load("data/testspc.RData")

wl <- as.numeric(colnames(testspc))
test_spc <- testspc[sample(1:nrow(testspc),2), wl>499 & wl<3001]

wopt <- mor.getwopt(test_spc$spc)
morbl <- mor.morbl(test_spc$spc,wopt) # simple MOR
imorbl <- mor.imorbl(test_spc$spc)    #  I-MOR
polybl <- spc.fit.poly.below(test_spc, poly.order = 5) # Current default
spc_sm <- prospectr::savitzkyGolay(test_spc$spc, m = 0, p=2, w=min(wopt))
imorbl_sm <- mor.imor(spc_sm)
plotspc(test_spc)
plotspc(as.hyperSpec(morbl,wl=wl(test_spc)), col = 'red', add=T)
plotspc(as.hyperSpec(imorbl,wl=wl(test_spc)), col = 'green', add=T)
plotspc(polybl, col = 'blue', add=T)
#plotspc(as.hyperSpec(spc_sm))
plotspc(as.hyperSpec(imorbl_sm),col='darkgreen', add=T) 



n <- 100
diff <- rep(NA,n)
s <- hspc[[51,,400~max,drop=FALSE]]
prevO <- mor.opening(s,1)
for (w in 2:n) {
  O <- mor.opening(s, w)
  diff[w] <- sum(abs(O-prevO))
  prevO <- O
}
which(diff==0)
plotspc(hspc[51,,400~max])
plotspc(as.hyperSpec(mor.morbl(hspc[[51,,400~max,drop=F]], 27)), add=T, col='red')
plotspc(as.hyperSpec(mor.morbl(hspc[[51,,400~max,drop=F]], 38)), add=T, col='red')
plotspc(as.hyperSpec(mor.morbl(hspc[[51,,400~max,drop=F]], 49)), add=T, col='red')
#plotspc(as.hyperSpec(mor.morbl(hspc[[51,,400~max,drop=F]], 53)), add=T, col='red')
plotspc(as.hyperSpec(mor.morbl(hspc[[51,,400~max,drop=F]], 61)), add=T, col='blue')
plotspc(as.hyperSpec(mor.morbl(hspc[[51,,400~max,drop=F]], 100)), add=T, col='blue')
plotspc(as.hyperSpec(mor.imorbl(hspc[[51,,400~max,drop=F]])), add=T, col='green')
#================ BENCHMARKS: ======================
library(rbenchmark)
test_spc <- testspc[sample(1:nrow(testspc),100),,500~3000]
benchmark(
  mor.getwopt(test_spc),
  mor.morbl(test_spc,50),
  mor.imorbl(test_spc),
  spc.fit.poly.below(as.hyperSpec(test_spc), poly.order = 4),
  replications = 5
)
#test replications elapsed relative user.self sys.self user.child sys.child
#1                    mor.getwopt(hspc$spc)            5 188.286  224.953   188.246    0.056          0         0
#3                     mor.imorbl(hspc$spc)            5 297.798  355.792   297.687    0.048          0         0
#2                  mor.morbl(hspc$spc, 50)            5  21.561   25.760    21.557    0.004          0         0
#4 spc.fit.poly.below(hspc, poly.order = 4)            5   0.837    1.000     0.837    0.000          0         0


#================ BENCHMARKS: apply vs for ======================
p <- ncol(spc)
spc <- matrix(rep(test_spc$spc,20),ncol = p)
w <- 26
y <- function (s) {
  news <- s
  for (j in 1:p) {
    news[j] <- min(s[max(1,j-w,na.rm = T):min(p,j+w,na.rm = T)])
  }
  return(news)
}
z <- function(spc){
  newspc <- spc
  for(i in 1:nrow(spc)){
    newspc[i,] <- y(spc[i,])
  }
}
benchmark(
apply(spc,1,y),
z(spc),
replications = 20
)