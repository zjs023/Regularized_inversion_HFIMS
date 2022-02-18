# Matlab tool for Regularized inversion of aerosol hygroscopic growth factor probability density function (GF-PDF)

This program, published with [Zhang et al. (2022)](https://amt.copernicus.org/preprints/amt-2021-334/), is developed to invert aerosol GF-PDF using regularized inversion methods including Tikhonov regularization and Twomey's method, from the humidity-controlled fast integrated mobility spectrometer measurements (HFIMS) measurements. For more details, please refer to the paper.

### _Getting started_

First, download the whole program package which includes the MATLAB [code](https://github.com/zjs023/Regularized_inversion_HFIMS/tree/master/m%20files) and related [data](https://github.com/zjs023/Regularized_inversion_HFIMS/tree/master/data/FIMS). If use git, clone the repository by:
```shell
git clone git://github.com/zjs023/Regularized_inversion_HFIMS
```
Change the current folder to `Regularized_inversion_HFIMS/` using the function `cd(...)` and add the MATLAB path in the Command Window using
```Matlab
addpath(genpath('./'));
```
Run the program by calling `GFPDF_inv_cmp` in Command Window, and it will automatically run the inversions of three predefined GF-PDFs and compare the results from different methods. Pay attention to the submodule switches (i.e., `option.` and `action.`) if adjustments are needed. 

### _Packages_

#### [GFPDF_inversion/](https://github.com/zjs023/Regularized_inversion_HFIMS/tree/master/m%20files/GFPDF_inversion) 
includes the main scripts for conducting the inversions. Six inversion methods were incorporated here, categorized as parametric and nonparametric inversion methods. The parametric methods assume a prior known distribution of the GF-PDF (i.e., multiple lognormal and piecewise linear function), and use least-squares fittings to search for the best solution, as detailed in [Wang et al. (2019)](https://www.tandfonline.com/doi/full/10.1080/02786826.2019.1628917). While nonparametric methods use a inverse mode introduced in [Zhang et al. (2022)](https://amt.copernicus.org/preprints/amt-2021-334/) and do not require a prior knowledge of the functional form of the GF-PDF. Nonparametric methods includes unregularized least-squares and regularized methods, such as Tikhonov regularization and Twomey's method. 

#### [utilities/](https://github.com/zjs023/Regularized_inversion_HFIMS/tree/master/m%20files/utilities)
includes subroutines called by the main scripts.

#### [plot/](https://github.com/zjs023/Regularized_inversion_HFIMS/tree/master/m%20files/plot)
includes the `confplot` function to plot the inverted GF-PDF with standard deviation bounds (code from [Michele Giugliano](https://www.mathworks.com/matlabcentral/fileexchange/2683-confplot)).

#### [regutools/](https://github.com/zjs023/Regularized_inversion_HFIMS/tree/master/m%20files/regutools)
includes the functions (e.g., `l_curve`, `ncsolve`, and ect.) for Tikhonov regularization, from a MATLAB package of Regularization Tools (code from [Per Christian Hansen](https://www.mathworks.com/matlabcentral/fileexchange/52-regtools?s_tid=srchtitle)).

----------------------------------------------------------------------
#### _Contact info_

This program was mainly written by Jiaoshi Zhang ([jiaoshi@wustl.edu](mailto:jiaoshi@wustl.edu)) with contributions from Dr. [@JianWang](https://scholar.google.com/citations?user=0yE2tSMAAAAJ&hl=en) and Dr. [@YangWang](https://scholar.google.com/citations?user=dkU1FrMAAAAJ&hl=en)
while at Washington University in St. Louis. Feel free to contact if you have any questions in running the program or want to contribute.

#### _References_

[Zhang, J., Wang, Y., Spielman, S., Hering, S., and Wang, J.: Regularized inversion of aerosol hygroscopic growth factor probability density function: Application to humidity-controlled fast integrated mobility spectrometer measurements, Atmos. Meas. Tech. Discuss. [preprint], https://doi.org/10.5194/amt-2021-334, in review, 2021.](https://amt.copernicus.org/preprints/amt-2021-334/)

[Wang, Y., Zheng, G., Spielman, S. R., Pinterich, T., Hering, S. V., and Wang, J.: Retrieval of high time resolution growth factor probability density function from a humidity-controlled fast integrated mobility spectrometer, Aerosol Science and Technology, 53, 1092-1106, 2019.](https://www.tandfonline.com/doi/full/10.1080/02786826.2019.1628917)

[Zhang, J., Spielman, S., Wang, Y., Zheng, G., Gong, X., Hering, S., and Wang, J.: Rapid measurement of RH-dependent aerosol hygroscopic growth using a humidity-controlled fast integrated mobility spectrometer (HFIMS), Atmos. Meas. Tech., 14, 5625-5635, 10.5194/amt-14-5625-2021, 2021.](https://amt.copernicus.org/articles/14/5625/2021/)

[Pinterich, T., Spielman, S. R., Hering, S., and Wang, J.: A water-based fast integrated mobility spectrometer (WFIMS) with enhanced dynamic size range, Aerosol Science and Technology, 51, 1212-1222, 10.1080/02786826.2017.1338664, 2017a.](https://www.tandfonline.com/doi/full/10.1080/02786826.2017.1338664)

[Pinterich, T., Spielman, S. R., Wang, Y., Hering, S. V., and Wang, J.: A humidity-controlled fast integrated mobility spectrometer (HFIMS) for rapid measurements of particle hygroscopic growth, Atmos. Meas. Tech., 10, 4915-4925, 10.5194/amt-10-4915-2017, 2017b.](https://amt.copernicus.org/articles/10/4915/2017/)

[Olfert, J. S., Kulkarni, P., and Wang, J.: Measuring aerosol size distributions with the fast integrated mobility spectrometer, Journal of Aerosol Science, 39, 940-956, https://doi.org/10.1016/j.jaerosci.2008.06.005, 2008.](https://www.sciencedirect.com/science/article/abs/pii/S002185020800116X)

[Sipkens, T. A., Olfert, J. S., and Rogak, S. N.: Inversion methods to determine two-dimensional aerosol mass-mobility distributions: A critical comparison of established methods, Journal of Aerosol Science, 140, 105484, https://doi.org/10.1016/j.jaerosci.2019.105484, 2020.](https://www.sciencedirect.com/science/article/abs/pii/S0021850219305889)

[Gysel, M., McFiggans, G. B., and Coe, H.: Inversion of tandem differential mobility analyser (TDMA) measurements, Journal of Aerosol Science, 40, 134-151, https://doi.org/10.1016/j.jaerosci.2008.07.013, 2009.](https://www.sciencedirect.com/science/article/abs/pii/S0021850208001778)

[Stolzenburg, M. R., and McMurry, P. H.: Equations Governing Single and Tandem DMA Configurations and a New Lognormal Approximation to the Transfer Function, Aerosol Science and Technology, 42, 421-432, 10.1080/02786820802157823, 2008.](https://www.tandfonline.com/doi/full/10.1080/02786820802157823)

[Hansen, P. C.: REGULARIZATION TOOLS: A Matlab package for analysis and solution of discrete ill-posed problems, Numerical Algorithms, 6, 1-35, 10.1007/BF02149761, 1994.](https://link.springer.com/article/10.1007/BF02149761)

[Markowski, G. R.: Improving Twomey's Algorithm for Inversion of Aerosol Measurement Data, Aerosol Science and Technology, 7, 127-141, 10.1080/02786828708959153, 1987.](https://www.tandfonline.com/doi/abs/10.1080/02786828708959153)

#### _How to cite_

Please cite [Zhang et al. (2022)](https://amt.copernicus.org/preprints/amt-2021-334/) if only the nonparametric inversion methods are used. One should also cite [Wang et al. (2019)](https://www.tandfonline.com/doi/full/10.1080/02786826.2019.1628917) if the parametric inversion methods are used. 
