# spectra-heplers

Set of useful tools for processing Raman / SERS / IR or other similar spectra.

**NOTE: For now it's raw and messy repo. Some functions use hyperSpec as input, some use matrix.**

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
