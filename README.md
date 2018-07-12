# spectra-heplers

Set of useful tools for processing vibrational (Raman / SERS / IR) or other similar spectra.

## Baseline correction

### `baseline.mor.R`
A method based on morphological operations (MOR) **and** fully automated method based on iterative morphological operations (I-MOR). I-MOR shows good results for "realworld" spectra. However, the output baseline has stairs-like shape which can harm spectral information. Both MOR and I-MOR code is C++ enhanced, which makes them quite fast. 

**References:**

- Liankui Dai and YunliangChen. *EXPRESS: An automated baseline correction method based on iterative morphological operations.* DOI:10.1177/0003702817752371

### `baseline.Goldinec.R`
Goldindec semiautomated baseline correction method based on regression with cost fucntion. The code of the function is optimized and transformed version of original Matlab code (see refs.). The function has been tested to give the same output as original Matlab code. 

However, on "real world" spectra the method does not perform much better that just polynom. 

**References:**

- Juntao Liu, Jianyang Sun, Xiuzhen Huang, Guojun Li, Binqiang Liu. *Goldindec: A Novel Algorithm for Raman Spectrum Baseline Correction.* DOI:10.1366/14-07798
- https://sourceforge.net/projects/transcriptomeassembly/files/Baseline_Correction/Goldindec.rar/

### `baseline.wavelet.R` 
This is simple wrap around original code (see refs.) to make baseline coorection done by one command. The method shows good results as well. However, there are about 8 parameters must be set. So, the method is far away from automatic baseline correction.

**References:**

- Z.M. Zhang, S. Chen, Y.Z. Liang, et al., An intelligent background-correction algorithm for highly fluorescent samples in Raman spectroscopy. Journal of Raman Spectroscopy 41 (6), 659 (2010).
- Original repo: https://github.com/zmzhang/baselineWavelet

### Examples:

Same works for vectors and hyperSpec objects
```{r}
# load('data/testspc.RData')

# source('baseline/baseline.mor.R')
bl.mor <- baseline.mor(testspc[1:10,], w=20)
bl.imor <- baseline.imor(testspc[1:10,])

# source('baseline/baseline.wavelet.R')
# by default scales=1:70, wavelet='mexh', gapTh=3, skip=2, SNR.Th=1, ridgeLength=5,lambda=1000, differences=1
bl.wavelet <- baseline.wavelet(testspc[1:10,])

# source('baseline/baseline.Goldindec.R')
# by default p (i.e. polynom order) = 4 and peak_ratio = 0.5
bl.Goldindec <- baseline.Goldindec(testspc[1:10,])

# Plot
wl <- as.numeric(colnames(testspc))
plot(x=wl, y=testspc[1,], type='l')
lines(x=wl, y=bl.mor[1,], col='red')
lines(x=wl, y=bl.imor[1,], col='darkred')
lines(x=wl, y=bl.wavelet[1,], col='green')
lines(x=wl, y=bl.Goldindec[1,], col='blue')
```

### Comparation:
| Method   | Speed[^1]         | Number of parameters[^2] | Pros | Cons |
|----------|-------------------|---|-----|------|
| Goldinec | 41.783434968 sec. | 2 | TBD | TBD |
|  Wavelet |  5.187343574 sec. | 8 | TBD | TBD |
|      MOR |  0.007592707 sec. | 1 | TBD | TBD |
|    I-MOR |  0.086890019 sec. | **0** | TBD | TBD |
|  Polynom |  0.019728261 sec. | 1 | TBD | TBD |

[^1]: median time for calculating baseline for 10 spectrum with 1770 wavelength points

[^2]: except tolerance levels

![Compare baseline correction methods][bl_correction_comparation]
![Compare baseline correction methods (zoom)][bl_correction_comparation_zoom]
![Compare baseline correction methods speed][bl_correction_comparation_speed]

## Other spectra preprocessing:

## `normalize.R` 

Normalize spectra using different type of normalizations: area, max, minmax (i.e. 01), vector, and peak.

**Examples:**
```{r}
# library(hyperSpec)
flu.normalized <- normalize(flu, 'area') # Makes area of each spectrum to be equal 1
flu.normalized <- normalize(flu, 'vector') # Makes 2-norm = sqrt(sum(x*x)) of each spectrum to be equal 1
laser.normalized <- normalize(laser, 'peak', peak.range=c(404,405.15)) # Makes intensity of the peak in range [404,405.15] to be equal 1
```
### `straightline.R`

The values in given data are replaced by straight line. 
This is one of ways that cosmic particles cotribution in Raman/SERS spectra usually handled.

**Examples**

Pretend we need to replace range between 900 and 1100 for the first ten spectra.
```{r}
# library(hyperSpec)
res <- straightline(chondro[1:10], from = 990, to = 1102)
plotspc(res)
```
See more description in the file.

[bl_correction_comparation]: https://github.com/rguliev/spectra-heplers/blob/master/baseline/tests_benchmarks_experiments/compare_all_methods.png
[bl_correction_comparation_zoom]: https://github.com/rguliev/spectra-heplers/blob/master/baseline/tests_benchmarks_experiments/compare_all_methods(1600-2300).png	
[bl_correction_comparation_speed]:  https://github.com/rguliev/spectra-heplers/blob/master/baseline/tests_benchmarks_experiments/compare_all_methods_speed.png	
