
:: copy this file to your case dir and launch !!!

@set HSFlowRoot=E:\science\devel\HSFlow-unstruct

@if %HSFLOW_CASE_DIR%.==. set HSFLOW_CASE_DIR=%HSFlowRoot%/test_case 

:: set paths to required libs

@set MPI_DIR=%HSFlowRoot%/lib/msmpi_v8

@set CGNS_DIR=%HSFlowRoot%/lib/cgns-vc14-rlz-win64

@set PATH=%PATH%;%INTEL_COMPILERS%/bin/intel64;

@start %HSFlowRoot%/HSFlow-unstruct.sln

===============================
