normalize <- function (spc, type = c('scale', 'area', 'vector', 'minmax', 'max', 'peak'), 
            peak.range = wl2i(spc, min~max), 
            scale.options = list(center = TRUE, scale = TRUE)) {
  switch(
    type,
      scale = {
        new.spc <- scale(spc, center = scale.options$center, scale = scale.options$scale)
      },
      area = {
        m <- apply(spc, 1, mean, na.rm = TRUE)
        new.spc <- sweep (spc, 1, m, "/")
      },
      vector = {
        norm.to <- apply(spc, 1, function (x) {sqrt(sum(x[!is.na(x)]*x[!is.na(x)]))})
        new.spc <- sweep (spc, 1, norm.to, "/")
      },
      minmax = {
        m <- apply(spc, 1, min, na.rm = TRUE)
        new.spc <- sweep (spc, 1, m, `-`)
        m <- apply(new.spc, 1, max, na.rm = TRUE)
        new.spc <- sweep (new.spc, 1, m, `/`)
      },
      max = {
        m <- apply(spc, 1, function(x) {max(abs(x),na.rm = TRUE)})
        new.spc <- sweep (spc, 1, m, `/`)
      },
      peak = {
        if (missing(peak.range)) {
          stop('wl.range must be setted when peak normalization is used!')
        }
        norm.to <- apply(spc[,, peak.range, wl.index = TRUE], 1, max, na.rm = TRUE)
        new.spc <- sweep (spc, 1, norm.to, "/")
      },
      {
        stop('Unknown type of normalization!')
      }
  )
  return(new.spc)
}