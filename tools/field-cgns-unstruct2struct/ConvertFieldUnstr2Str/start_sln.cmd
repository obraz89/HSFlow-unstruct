
:: copy this file to your case dir and launch !!!

@set HSFlowLib=E:\science\devel\HSFlow-unstruct\lib

@set HSFLOW_CASE_DIR=%CD% 

:: set paths to required libs

@set MPI_DIR=E:\lib\msmpi_10.1.2

@set CGNS_DIR=%HSFlowLib%/cgns-vc14-rlz-win64

@set PETSC_DIR=%HSFlowLib%/petsc_3.8.2-vc14-rlz

@set MKL_DIR=%HSFlowLib%/lib/mkl/lib/intel64

@start E:\science\devel\HSFlow-unstruct\tools\field-cgns-unstruct2struct\ConvertFieldUnstr2Str\ConvertFieldUnstr2Str.sln

===============================
