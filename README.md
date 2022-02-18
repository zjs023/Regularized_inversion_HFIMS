# Matlab tool for Regularized inversion of aerosol hygroscopic growth factor probability density function (GF-PDF)

This program, published with [Zhang et al. (2022)](https://amt.copernicus.org/preprints/amt-2021-334/), is developed to invert aerosol GF-PDF using regularized inversion methods including Tikhonov, and Twomey's method, from the humidity-controlled fast integrated mobility spectrometer measurements (HFIMS) measurements. For more details, please refer to the paper.

## Table of contents
[Getting started](#getting-started)

First, download the whole program package which includes the MATLAB [code](https://github.com/zjs023/Regularized_inversion_HFIMS/tree/master/m%20files) and related [data](https://github.com/zjs023/Regularized_inversion_HFIMS/tree/master/data/FIMS). If use git, clone the repository by:
```shell
git clone git://github.com/zjs023/Regularized_inversion_HFIMS
```
Change the current folder to `Regularized_inversion_HFIMS/` using the function `cd(...)` and add the MATLAB path in the Command Window using
```Matlab
addpath(genpath('./'));
```
Run the program by calling `GFPDF_inv_cmp` in Command Window, and it will automatically run the inversions of three predefined GF-PDFs and compare the results from different methods.

[Packages](#Packages)

##  **[data/FIMS/](https://github.com/zjs023/Regularized_inversion_HFIMS/tree/master/data/FIMS)
