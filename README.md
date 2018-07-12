# spectra-heplers

Set of useful tools for processing vibrational (Raman / SERS / IR) or other similar spectra.

## Baseline correction

### `baseline.mor.R` and `baseline.imor.R`
A method based on morphological operations (MOR) and fully automated method based on iterative morphological operations (I-MOR). I-MOR shows good results for "realworld" spectra. However, the output baseline has stairs-like shape which can harm spectral information.

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

### `baseline.anchor.R`

# Finished functions:
## `normalize.R` 

Normalize spectra using different type of normalizations: area, max, minmax (i.e. 01), vector, and peak.

**Usage:**
```{r}
# library(hyperSpec)
flu.normalized <- normalize(flu, 'area') # Makes area of each spectrum to be equal 1
flu.normalized <- normalize(flu, 'vector') # Makes 2-norm = sqrt(sum(x*x)) of each spectrum to be equal 1
laser.normalized <- normalize(laser, 'peak', peak.range=c(404,405.15)) # Makes intensity of the peak in range [404,405.15] to be equal 1
```
See more description in the file.


## `straightline.R`

The values in given data are replaced by straight line. 
This is one of ways that cosmic particles cotribution in Raman/SERS spectra usually handled.

**Usage**

Pretend we need to replace range between 900 and 1100 for the first ten spectra.
```{r}
# library(hyperSpec)
res <- straightline(chondro[1:10], from = 990, to = 1102)
plotspc(res)
```
See more description in the file.

# TODO:
- add descriptions for each function
- standardize input format of spectra data
