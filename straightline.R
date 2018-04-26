straightline <- 
  function(spc, wl.left, wl.right) {
    indx.left <- wl2i(s,  wl.left)[1]
    indx.right <- wl2i(s, wl.right)[1]
    x <- c(wl(s)[indx.left], wl(s)[indx.right])
    
    for (i in 1:nrow(s)){
      y <- s$spc[i,c(indx.left, indx.right)]
      k <- (y[2]-y[1])/(x[2]-x[1])
      b <- (x[2]*y[1]-x[1]*y[2])/(x[2]-x[1])
      s$spc[i,indx.left:indx.right] <- wl(s)[indx.left:indx.right]*k + b
    }
    
    return(s)
  }