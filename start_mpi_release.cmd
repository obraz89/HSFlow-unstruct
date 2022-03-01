@set HSFlowRoot=S:\Science\HSFlow-unstruct\x64\Release

@set HSFLOW_CASE_DIR=%CD%

mpiexec -n 1 %HSFlowRoot%\HSFlow-unstruct.exe -l log.txt %CD%
