========================================================================
    CONSOLE APPLICATION : HSFlow-unstruct Project Overview
========================================================================

HSFlow-unstruct is a solver for 3d Euler and Navier-Stokes equations 
on unstructured grids.

It is an educational project to show main features of finite-volume Godunov-type schemes


/////////////////////////////////////////////////////////////////////////////
Installation notes:

1. Clone the repository to your local machine (simple way is to use GitHub Desktop)

2. Copy lib directory inside the root of your repository (provided separately)

From lib/msmpi_v8 run MSMpiSetup.exe and msmpisdk.msi
This will install Microsoft mpi on your machine (binaries and sdk)

3. Create folder for your test case, copy run scripts to your case folder: 
start_sln.cmd, start_mpi_debug.cmd, start_mpi_release.cmd

4. Configure your project (settings dir, grids, bcs etc), 
run using scripts

/////////////////////////////////////////////////////////////////////////////
Other notes:

AppWizard uses "TODO:" comments to indicate parts of the source code you
should add to or customize.

/////////////////////////////////////////////////////////////////////////////
