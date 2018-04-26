find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

myfindPeaks <- 
  function (x, thresh=0.05, span=0.25, lspan=0.05, noisey=TRUE)
  {
    n <- length(x)
    y <- x
    mu.y.loc <- y
    if(noisey)
    {
      mu.y.loc <- (x[1:(n-2)] + x[2:(n-1)] + x[3:n])/3
      mu.y.loc <- c(mu.y.loc[1], mu.y.loc, mu.y.loc[n-2])
    }
    y.loess <- loess(x~I(1:n), span=span)
    y <- y.loess[[2]]
    sig.y <- var(y.loess$resid, na.rm=TRUE)^0.5
    DX.1 <- sign(diff(mu.y.loc, na.pad = FALSE))
    pks <- which(diff(DX.1, na.pad = FALSE) < 0 & DX.1[-(n-1)] > 0) + 1
    out <- pks
    if(noisey)
    {
      n.w <- floor(lspan*n/2)
      out <- NULL
      for(pk in pks)
      {
        inner <- (pk-n.w):(pk+n.w)
        outer <- c((pk-2*n.w):(pk-n.w),(pk+2*n.w):(pk+n.w))
        mu.y.outer <- mean(y[outer])
        if(!is.na(mu.y.outer)) 
          if (mean(y[inner])-mu.y.outer > thresh*sig.y) out <- c(out, pk)
      }
    }
    out
  }

findpeaks <- function(vec,bw=1,x.coo=c(1:length(vec)))
{
  pos.x.max <- NULL
  pos.y.max <- NULL
  pos.x.min <- NULL
  pos.y.min <- NULL
  for(i in 1:(length(vec)-1)) 	{ 		
    if((i+1+bw)>length(vec)){
      sup.stop <- length(vec)}
    else{
      sup.stop <- i+1+bw
    }
    if((i-bw)<1){inf.stop <- 1}else{inf.stop <- i-bw}
    subset.sup <- vec[(i+1):sup.stop]
    subset.inf <- vec[inf.stop:(i-1)]
    
    is.max   <- sum(subset.inf > vec[i]) == 0
    is.nomin <- sum(subset.sup > vec[i]) == 0
    
    no.max   <- sum(subset.inf > vec[i]) == length(subset.inf)
    no.nomin <- sum(subset.sup > vec[i]) == length(subset.sup)
    
    if(is.max & is.nomin){
      pos.x.max <- c(pos.x.max,x.coo[i])
      pos.y.max <- c(pos.y.max,vec[i])
    }
    if(no.max & no.nomin){
      pos.x.min <- c(pos.x.min,x.coo[i])
      pos.y.min <- c(pos.y.min,vec[i])
    }
  }
  return(list('pos.x.max' = pos.x.max, 'pos.y.max' = pos.y.max, 'pos.x.min' = pos.x.min, 'pos.y.min' = pos.y.min))
}