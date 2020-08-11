# Overview

This repository contains data and scripts for reproducing the results accompanying the manuscript

### Fitness inference from complex evolutionary histories with genetic linkage  
Muhammad S. Sohail<sup>1,</sup>\*, Raymond H. Y. Louie<sup>1-4,</sup>\*, Matthew R. McKay<sup>1,5,#</sup> and John P. Barton<sup>6,#</sup>

<sup>1</sup> Department of Electronic and Computer Engineering, Hong Kong University of Science and Technology  
<sup>2</sup> Institute for Advanced Study, Hong Kong University of Science and Technology  
<sup>3</sup> Kirby Institute, University of New South Wales  
<sup>4</sup> School of Medical Sciences, University of New South Wales  
<sup>5</sup> Department of Chemical and Biological Engineering, Hong Kong University of Science and Technology  
<sup>6</sup> Department of Physics and Astronomy, University of California, Riverside  
\* Equal contributions  
<sup>#</sup> correspondence to [m.mckay@ust.hk](mailto:m.mckay@ust.hk) and [john.barton@ucr.edu](mailto:john.barton@ucr.edu)  

# Contents

Scripts for generating and analyzing the simulation data can be found in the `simulation-analysis.ipynb` notebook, and references therein. Scripts for processing and analyzing the HIV-1 data are contained in the `HIV-analysis.ipynb` notebook. Finally, scripts for analysis and figures contained in the manuscript are located in the `figures.ipynb` notebook.  

Due to the large size and number of some files generated by simulations and by the interim analysis of HIV-1 data, some data has been stored in a compressed format using Zenodo. To access the full set of data, navigate to the [Zenodo record](https://zenodo.org/record/3979632). Then download and extract the contents of the archives `src-MPL-HIV.tar.gz`, `src-wfsim-data.tar.gz`, and `src-Matlab-data.tar.gz` into the folders `src/MPL/HIV`, `src/wfsim/data`, and `src/Matlab`, respectively.

# MPL

This repository includes code for inferring selection coefficients from temporal genetic data using the marginal path likelihood (MPL) method. Code implementing MPL in C++ is located in the `src/MPL` directory.

### Software dependencies

MPL makes use of C++11 and the [GNU Scientific Library](https://www.gnu.org/software/gsl/). The code has been tested using GSL version 2.6. Compiling MPL should require a few seconds on a typical desktop computer.

### Demonstration

The shell script `src/MPL/run-test.sh` will compile and run MPL on an example set of simulation data (see Figure 1 in the manuscript). Estimated selection coefficients are stored in `src/MPL/test/example_0_T400_ns1000_dt1_MPL.dat`. This file should match the file of the same name in the `src/MPL/out` directory. Running this demonstration should require a few seconds on a typical desktop computer.

# License

This repository is dual licensed as [GPL-3.0](LICENSE-GPL) (source code) and [CC0 1.0](LICENSE-CC0) (figures, documentation, and our presentation of the data).
