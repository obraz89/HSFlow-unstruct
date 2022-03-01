@set HSFlowRoot=S:\Science\HSFlow-unstruct\x64\Debug

@set HSFLOW_CASE_DIR=%CD%

@set PATH=%PATH%; 

start mpiexec -n 1 %HSFlowRoot%\HSFlow-unstruct.exe -l log.txt -w 60 %CD%

::mpiexec -n 2 S:\Science\HSFlow-unstruct\x64\Debug\HSFlow-unstruct.exe -l log.txt %CD%
