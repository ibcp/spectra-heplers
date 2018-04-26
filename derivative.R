derivative <- 
  function (spc, order = 1, p = 2, w = 9, normalize = TRUE) {
    spc.der <- savitzkyGolay(spc$spc, m = order, p = p, w = w)
    spc.der <- as.hyperSpec(X = spc.der)
    if (normalize) {
      spc.der <- normalize(spc.der, type = 'vector')
    }
    return(spc.der)
  }