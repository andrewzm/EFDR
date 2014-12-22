EFDR
====

Enhanced False Discovery Rate (EFDR) is a non-parameteric hypothesting testing procedure used for anomaly detection in noisy spatial signals. The approach, published in Shen et al. (2002), is an extension on standard multiple hypothesis testing approaches (for example those employing the Bonferroni correction, or the standard FDR approach) by reducing the number of tests carried out. A reduced number of tests results in an approach which has more power (a lower Type II error), and hence an increase ability to detect anomalies. EFDR has proven to be vastly superior to the Bonferroni correction and standard FDR: the interested reader is referred to Shen et al. (2002) for a detailed study.

The package `EFDR` contains the required functions to carry out EFDR in a practical setting. It allows for the possibility of a parallel backend (since the computations are relatively intensive), contains basic interpolation methods to grid data which is spatially irregular, and also contains the more standard methods such as detection using the Bonferroni correction, the standard FDR and the Largest Order Statistic method. 

Examples
--------

For tutorials on how to use this package, please visit the [project development page](http://www.github.com/andrewzm/EFDR) or visit the project index HTML page by typing `help(package = "EFDR")` after installation.

Installation
------------

- Installation of most dependencies is straightforward. If you are on Linux, for the package `gstat` you will also need to install `libgeos-dev` which is widely available in linux repos. In Ubuntu this is available using `apt-get install libgeos-dev`. 
- To install the latest (development) version please install and load `devtools` and in an R console type `install_github("EFDR","andrewzm")`.


Open issues
---------------

- Neighbourhood detection is currently sped up using a parallel backend. This, however, will remain slow on a standard desktop. Suggest implementing in C.
- Image padding should be integrated, possibly using conditional simulation, see Pavlicová et al. (2008).


References
----------

Pavlicová, M., Santner,  T. J., and Cressie, N. "Detecting signals in FMRI data using powerful FDR procedures." Statistics and its interface 1 (2008): 23-32.

Shen, X., Huang, H.-C., and Cressie, N. "Nonparametric hypothesis testing for a spatial signal." Journal of the American Statistical Association 97.460 (2002): 1122-1140.

