## BULK LOADING
# Examples: 
#   s <- read.spcfiles(Sys.glob('path_to_dir/*.spc'), func = read.spc)
#   s <- read.spcfiles(Sys.glob('path_to_dir/*.spc'), func = read.jdx)
read.spcfiles <- 
  function (files, 
            func = function(f) { 
                      read.csv(f, header = F, colClasses = 'numeric', row.names = 1) %>% 
                        t %>% as.hyperSpec %>% orderwl }, 
            v = FALSE) {
    for (i in 1:length(files)) {
      filename <-  basename(files[i])
      spc <-func(files[i])
      if (v) {
        if (is (spc, "hyperSpec"))
          cat (filename, ": ", nrow (spc), " spectrum(a), ", nwl (spc), " data pts / spc.\n", sep = "")
        else
          cat (filename, ": list of ", length (spc), " spectra, ",
               paste (range (sapply (spc, nwl)), collapse = " - "),
               " data pts / spc\n", sep = "")
      }
      
      spc$filename <- rep(filename, nrow(spc))
      spc$filename_full <- rep(files[i], nrow(spc))
      
      if (i == 1) {
        spc_arr <-  spc
      }
      else {
        spc_arr <- hyperSpec::collapse(spc_arr, spc, wl.tolerance = 0.1)
      }
    }
    spc_arr <- spc_arr[,c('spc','filename','filename_full'),]
    spc_arr <- orderwl(spc_arr)
    return(spc_arr)
  }