library(hexView)
library(hyperSpec)

read.opus <- 
    function(file.name, sp=NULL, speclib="ICRAF", MID) {
        if(!(speclib=="ICRAF"|speclib=="New")){ stop("'speclib' must be one of the following: 'ICRAF' or 'New'") }
        if(file.exists(file.name)){
            ## Read metainfo
            try( pa <- hexView::readRaw(file.name, offset = 0, nbytes = file.info(file.name)$size, human = "char", size = 1, endian = "little"), silent=TRUE )
            if(!class(.Last.value)[1]=="try-error"){
                pr <- pa$fileRaw
                ## Get source of instrument
                ins <- grepRaw("INS", pr, all=TRUE)
                ins <- readRaw(file.name, offset = ins[length(ins)]+7, nbytes = 3, human = "char", size = 1, endian = "little")
                ins <- blockString(ins)
                ## Get source of infrared to know if NIR or MIR
                src <- grepRaw("SRC", pr, all=TRUE)
                src <- readRaw(file.name, offset = src[length(src)]+7, nbytes = 3, human = "char", size = 1, endian = "little")
                src <- blockString(src)
                instr.range <- tolower(paste(ins, src, sep="-"))
                ## Get Beam Splitter
                bms <- grepRaw("BMS", pr, all=TRUE)
                bms <- readRaw(file.name, offset = bms[length(bms)]+7, nbytes = 4, human = "char", size = 1, endian = "little")
                bms <- blockString(bms)
                ## Wavenumbers for MIR spectra from Tensor are assigned prefix "m", MIR spectra from Alpha prefixed "a"; for NIR MPA "n"
                if(instr.range=="ten-off"){ instr.range="ten-mir"} ## AS: Old ten-mir written as tensor-27
                # if(speclib=="ICRAF"){ 
                #     icraf.htsxt <- c(3578, 7497.964, 599.76)
                #     icraf.alpha <- c(2542, 3998.12872, 399.387991)
                #     icraf.mpa   <- c(2307, 12493.2, 3598.69)
                #     if(instr.range=="ten-mir"){
                #         wb <- rev(seq(icraf.htsxt[3], icraf.htsxt[2], (icraf.htsxt[2]-icraf.htsxt[3])/(icraf.htsxt[1]-1)))
                #     }
                #     if(instr.range=="alp-mir"){
                #         wb <- rev(seq(icraf.alpha[3], icraf.alpha[2], (icraf.alpha[2]-icraf.alpha[3])/(icraf.alpha[1]-1)))
                #     }
                #     if(instr.range=="mpa-nir"){
                #         wb <- rev(seq(icraf.mpa[3], icraf.mpa[2], (icraf.mpa[2]-icraf.mpa[3])/(icraf.mpa[1]-1)))
                #     }
                #     if(bms=="ZnSe"){
                #         pref="a"
                #         wb <- rev(seq(499.8151, 3996.4810, (3996.4810-499.8151)/(1715-1)))
                #     }
                # }
                if(!(instr.range=="ten-mir"|instr.range=="alp-mir"|instr.range=="mpa-nir"|instr.range=="ten-off"|bms=="ZnSe")){ stop("Unknown file format. See '?read.opus' for more info.") }  
                
                ## Get positions where the following parameters are found in the file  
                z    <- grepRaw("ZFF",pr,all=TRUE)[1]+5
                re   <- grepRaw("RES",pr,all=TRUE)[1]+5
                snm  <- grepRaw("SNM",pr,all=TRUE)[1]+7
                dat  <- grepRaw("DAT",pr,all=TRUE)[1]+7
                lwn  <- grepRaw("LWN",pr,all=TRUE)[1]+7
                fx   <- grepRaw("FXV",pr,all=TRUE)[1]+7 #[3]
                lx   <- grepRaw("LXV",pr,all=TRUE)[1]+7 #[3]
                npt0 <- grepRaw("NPT",pr,all=TRUE)[1]+3 #[2]
                npt1 <- grepRaw("NPT",pr,all=TRUE)[1]+7 #[3]
                mxy  <- grepRaw("MXY",pr,all=TRUE)[1]+7 
                mny  <- grepRaw("MNY",pr,all=TRUE)[1]+7  #[3]
                end  <- grepRaw("END",pr,all=TRUE)+11
                tim  <- grepRaw("TIM",pr,all=TRUE)+11
                ## calculate end and start of each block:
                offs <- sapply(5:10, function(x){end[x]})
                byts <- diff(offs)
                ZFF <- readRaw(file.name, offset=z, nbytes=4, human="int", size=2)[[5]][1]
                RES <- readRaw(file.name, offset=re, nbytes=4, human="int", size=2)[[5]][1]
                snm.lab.material <- blockString(readRaw(file.name, offset = snm, nbytes = 22, human = "char", size = 1, endian = "little"))
                if(!nzchar(snm.lab.material)){
                    SSN <- ""
                    Material <- ""
                    warning("Product name not found inside OPUS file...")
                } else {
                    if(!length(grep(snm.lab.material, pattern=";"))==0){
                        snm.lab.material <- as.vector(strsplit(snm.lab.material,";"))[[1]]
                        SSN <- paste0(snm.lab.material[2], snm.lab.material[1])
                        Material <- snm.lab.material[3]
                    } 
                    else {
                        if(!length(grep(snm.lab.material, pattern="_"))==0){ ## AS: ICR_02182
                            SSN <- sub("_", "", snm.lab.material)
                            Material <- ""          
                        }
                        else {
                            if(!length(snm.lab.material)==0){ 
                                SSN <- snm.lab.material
                                Material <- ""          
                            }
                        }
                    }   
                }
                ## Set three SSN first three characters to lower
                SSN <- paste0(tolower(substr(SSN,1,3)), substr(SSN,4,20))
                Scandate <- blockString(readRaw(file.name, offset = dat, nbytes = 10, human = "char", size = 1, endian = "little"))
                Scantime <- blockString(readRaw(file.name, offset = tim-4, nbytes = 8, human = "char", size = 1, endian = "little"))
                Scandate <- paste(Scandate,Scantime)
                LWN <- readRaw(file.name, offset=lwn, nbytes=8, human="real", size=8)[[5]][1]
                
                ## Get number of data points for each spectra data block   		
                NPT0 <- readRaw(file.name, offset=npt0, nbytes=12, human="int", size=4)[[5]][2]
                NPT1 <- readRaw(file.name, offset=npt1, nbytes=4, human="int", size=4)[[5]][1]
                fxv <- readRaw(file.name, offset=fx, nbytes=16, human="real", size=8)[[5]][1] ## fxv:	Frequency of first point
                lxv <- readRaw(file.name, offset=lx, nbytes=16, human="real", size=8)[[5]][1] ## lxv:	Frequency of last point
                Wavenumbers <- rev(seq(lxv, fxv, (fxv-lxv)/(NPT1-1)))
                ## Read all through all the data blocks inside the OPUS file:
                nbytes1 <- NPT0*4 ## initial parameters
                smxa <- c()
                smna <- c()
                nbytes.f <- NPT1*4
                # if(offs[1]<2000){
                #     offs.f <- offs[3]
                # }
                # 
                # if(offs[1]>20000){
                #     offs.f<-offs[2]
                # }
                offs.f <- end[5]
                ## Selected spectra block
                # opus.p <- readRaw(file.name,width=NULL,offset=offs.f-4,nbytes=nbytes.f,human="real",size=4,endian="little")
                opus.p <- readRaw(file.name,width=NULL,offset=offs.f, nbytes=nbytes.f,human="real",size=4,endian="little") 
                spectra <- opus.p[[5]]
                
                ## Make compatible to ICRAF spectra:
                # if(speclib=="ICRAF"){
                #     ## TH: Is spline fitting necessary? AS: Yes, this standardizes all spectral data points to conform to ICRAF spectral library
                #     spectra <- spline(Wavenumbers, spectra, xout=wb, method="natural")$y
                #     Wavenumbers <- wb
                # }
                
                ## Add meta ID
                if(missing(MID)){ MID <- paste(Sys.getenv(c("USERDNSDOMAIN"))[[1]], instr.range, sep="_") } 
                
                ## create data.frames:
                ## Combine the above parameters
                out <- new ("hyperSpec", spc = matrix(spectra, nrow = 1), 
                            wavelength = Wavenumbers, 
                            data = data.frame(  File.Name = file.name,
                                                SAMPLEID  = SSN, 
                                                Material  = Material, 
                                                Zero.Filing = ZFF, 
                                                Resolution  = RES, 
                                                LWN         = LWN, 
                                                DateTime    = as.POSIXct(Scandate, format="%d/%m/%Y %H:%M:%S "), 
                                                MID = MID),
                            label = list(.wavelength = expression("Wavenumbers, cm"^-1), spc = 'Absorabance'))
                return(out)
            }
        } else {
            warning(paste("File",file.name,"does not exist"))
        }
    }