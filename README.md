EFDR
====

Enhanced False Discovery Rate (EFDR) is a non-parameteric hypothesting testing procedure used for anomaly detection in noisy spatial signals. The approach, published in Shen et al. (2002), is an extension on standard multiple hypothesis testing approaches (for example those employing the Bonferroni correction, or the standard FDR approach) by reducing the number of tests carried out. A reduced number of tests results in an approach which has more power (a lower Type II error), and hence an increase ability to detect anomalies. EFDR has proven to be vastly superior to the Bonferroni correction and standard FDR: the interested reader is referred to Shen et al. (2002) for a detailed study.

The package `EFDR` contains the required functions to carry out EFDR in a practical setting. It allows for the possibility of a parallel backend (since the computations are relatively intensive), contains basic interpolation methods to grid data which is spatially irregular, and also contains the more standard methods such as detection using the Bonferroni correction, the standard FDR and the Largest Order Statistic method. 

This package will soon be available on CRAN but until then, to install you will need to have `devtools` installed and loaded. See the section *Installation* for more details. Please do not hesitate to contact me if you have any queries.

Examples
--------

For a tutorial on how to use this package, please click [here](http://htmlpreview.github.io/?https://github.com/andrewzm/EFDR/blob/master/docs/EFDR.html). For a case study on the detection of the El Niño and El Niña phenomena, please click [here](http://htmlpreview.github.io/?https://github.com/andrewzm/EFDR/blob/master/docs/EFDR_SST.html).

Installation
------------

- To install this package, first load `devtools`, then type `install_github('andrewzm/EFDR', dependencies=T, build_vignettes=F)`. Installation of most dependencies is straightforward. If you are on Linux, for the package `gstat` you will also need to install `libgeos-dev` which is widely available in linux repos. In Ubuntu this is available using `apt-get install libgeos-dev`.


Open issues
---------------

- Neighbourhood detection is currently sped up using a parallel backend. This, however, will remain slow on a standard desktop. Suggest implementing in C.
- Image padding should be integrated, possibly using conditional simulation, see Pavlicová et al. (2008).


References
----------

Pavlicová, Martina, Thomas J. Santner, and Noel Cressie. "Detecting signals in FMRI data using powerful FDR procedures." Statistics and its interface 1 (2008): 23-32.

Shen, Xiaotong, Hsin-Cheng Huang, and Noel Cressie. "Nonparametric hypothesis testing for a spatial signal." Journal of the American Statistical Association 97.460 (2002): 1122-1140.

