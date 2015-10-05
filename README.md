# PiPlusAnalysis

All the steps to make Pi Plus Analysis (Multilpicity Ratio and Braodening) are detailed here.

To compile the library please make sure you have installed make and imake (xutils-dev in Ubuntu)

First it is necessary to set the enviorment variables

setenv CLASTOOL ClasToolDirectory
setenv OS_NAME Linux
setenv CLAS_PACK CLAS_PACKDirectory
setenv CLAS_LIB CLAS_LIBDirectory
setenv ANALYSER ANALYSERDirectory/analysis_lib
setenv MYSQL_INCLUDE $MYSQLINC
setenv MYSQL_LIB $MYSQLIB

To compile HAPRAD_CPP it is necessary to have libMathMore in root, to do this libgsl need to be installed and
in configure set --enable-gsl-shared
(make sure configure output is "Checking whether to build libMathMore ... yes") 
It is also necessary to install cernlib package
