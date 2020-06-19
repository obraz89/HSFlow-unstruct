
:: copy this file to your case dir and launch !!!

@set HSFlowRoot=E:\science\devel\HSFlow-unstruct

@set HSFLOW_CASE_DIR=%CD% 

:: set paths to required libs

@set MPI_DIR=E:\lib\msmpi_10.1.2

@set CGNS_DIR=%HSFlowRoot%/lib/cgns-vc14-rlz-win64

@set PETSC_DIR=%HSFlowRoot%/lib/petsc_3.8.2-vc14-rlz

@set MKL_DIR=%HSFlowRoot%/lib/mkl/lib/intel64

@start %HSFlowRoot%/HSFlow-unstruct.sln

===============================
