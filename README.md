# SIDR
A statistical framework designed to assess the reproducibility of Hi-C chromatin interactions identified by peak callers such as **FitHiC2** and **HiC-DC+**. 

# Installation

1. Install the package from Github:
 
```r
devtools::install_github("qunhualilab/SIDR")
```
2. Install the package from source:

```r
install.packages("/PATH/TO/SOURCE/SIDR_1.0.tar.gz",repos = NULL, type = "source")
```

# Input

**SIDR** takes a pair of Hi-C replicate files as input. Each Hi-C file contains **p-values** derived from peak-calling methods, which measure the significance of interactions between paired bins at a fixed resolution.  

- Each bin is represented by its **chromosome** and **start-point coordinate**.  
- Each file should have **five columns**, separated by **space** or **tab**.  

**Column specification:**

| Column | Description |
|--------|-------------|
| 1      | Chromosome name of the first bin |
| 2      | Coordinate of the first bin |
| 3      | Chromosome name of the second bin |
| 4      | Coordinate of the second bin |
| 5      | P-value of interaction between the two bins |

**Example:** (the header should be included)

| chromosomeI | fragmentI | chromosomeJ | fragmentJ | pvalue      |
|-------------|-----------|-------------|-----------|------------|
| chr18       | 50000     | chr18       | 100000    | 1.0000000  |
| chr18       | 50000     | chr18       | 150000    | 0.9999981  |
| ...         | ...       | ...         | ...       | ...        |


# Example Run

This section demonstrates a typical workflow using SIDR with two Hi-C replicate files.

---

## Step 1: Merge Hi-C Replicates

To assess reproducibility, **SIDR** first merges the two Hi-C replicate files.  
Use the `mergeHiC` function, specifying the path to each replicate file and the chromosomes to analyze.

```r
# Merge two Hi-C replicate files for chromosome 18
data <- mergeHiC(
  "/PATH/TO/SOURCE/rep1.txt",
  "/PATH/TO/SOURCE/rep2.txt",
  chrI = "chr18",
  chrJ = "chr18"
)
```
Only interactions present in both replicates for the specified chromosomes are retained.

## Step 2: Stratify Hi-C Data by Genomic Distance

Divide the merged Hi-C data into distance strata. This is required for SIDR analysis.

```r
# Divide Hi-C data into 5 strata by genomic distance
ind <- stratifyData(data, ns = 5)
```

## Step 3: Determine the initial parameters 

Compute initial parameter values for the SIDR model. You can skip this step if you would use different initial values.
```r
# Determine the initial values of the parameters
init <- choosePara(
  data, 
  ns = 5,      # number of strata
  ind = ind,   # strata indices
  prp = 0.6,   # initial proportion of reproducible interactions
  rho = 0.15   # initial correlation
)
```
## Step 4: Fit the SIDR model

Fit the stratified IDR model using the initialized parameters.

```r
# Prepare initial parameter vector
thet.ini <- c(
  logit(init$mixps), 
  init$mus, 
  init$sigma, 
  logit(init$rho), 
  init$omega
)

# Fit the SIDR model
result <- fit_IDR_stratified(
  data, 
  ns = 5, 
  ind = ind, 
  thet.ini = thet.ini
)
```
