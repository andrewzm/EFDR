EFDR
====

Enhanced False Detection Rate (EFDR) is a non-parameteric hypothesting testing procedure used for anomaly detection in noisy spatial signals. The approach, published in Shen et al. (2002), is an extension on standard multiple hypothesis testing approaches (for example those employing the Bonferroni correction, or the standard FDR approach) by reducing the number of tests carried out. A reduced number of tests results in an approach which has more power (a lower Type II error), and hence an increase ability to detect anomalies. EFDR has proven to be vastly superior to the Bonferroni correction and standard FDR: the interested reader is referred to Shen et al. (2002) for a detailed study.

The package `EFDR` contains the required functions to carry out EFDR in a practical setting. It allows for the possibility of a parallel backend (since the computations are relatively intensive), contains basic interpolation methods to grid data which is spatially irregular, and also contains the more standard methods such as detection using the Bonferroni correction, the standard FDR and the Largest Order Statistic method. For a tutorial on how to use this package, please visit http://nbviewer.ipython.org/github/andrewzm/EFDR/blob/master/ipython/EFDR.ipynb

This package is still under development and is not available on CRAN. To install you will need to have `devtools` installed and loaded. Then type in  `install_github('andrewzm/EFDR')`. Please do not hesitate to contact me if you have any queries.

Installation
------------

- To install this package, first load `devtools`, then type `install_github('andrewzm/EFDR', dependencies=T)`. Installation of most dependencies is straightforward. If on Linux, for the package `gstat` you will also need to install `libgeos-dev` which is widely available in linux repos. In Ubuntu this is available using `apt-get install libgeos-dev`.

Open issues
---------------

- Neighbourhood detection is currently sped up using a parallel backend. This, however, will remain slow on a standard desktop. Suggest implementing in C.
- The package `waveslim` produces erroneous results when carrying out a DWT on an image which is not square. As a result, `EFDR` will not work on images which are not square.  Suggest altering the `dwt.2d` interface.


Package details
---------------

Package: EFDR

Type: Package

Title: Carried out Enhanced False Detection Rate in the wavelet domain

Version: 1.0

Date: 2014-10-22

Author: Andrew Zammit Mangion and Hsin-Cheng Huang

Maintainer: Andrew Zammit Mangion <azm [@] uow.edu.au>

Description: EFDR is a tool to detect anomalies in an image. The image is first
    transformed into the wavelet domain, the coefficients at each resolution
    are standardised, and standard statistical tests (in a multiple hypothesis
    testing setting) are carried out to find the anomalies. The power of EFDR
    exceeds that of standard FDR, which would carry out tests on every wavelet
    coefficient: EFDR choose which wavelets to test based on a criteria
    described in Shen, Xiaotong, Hsin-Cheng Huang, and Noel Cressie.
    "Nonparametric hypothesis testing for a spatial signal." Journal of the
    American Statistical Association 97.460 (2002): 1122-1140.

Imports:
    Matrix,
    waveslim,
    foreach,
    doMC,
    parallel,
    gstat,
    tidyr

License: GPL (>= 2)

NeedsCompilation: no
