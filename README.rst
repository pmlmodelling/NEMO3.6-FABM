===================================================================
Fork of NEMO 3.6 stable branch with modifcations for shelf systems.
===================================================================

You are on the shelf-enabled branch `feat/shelf-enabled`.

Compilation using ``makenemo``
==============================

To compile on a typical PML workstation using the ``makenemo`` tool provided with NEMO::

  module load mpi #required on fedora
  export XIOS_HOME #set to the basepath of XIOS-1 (compile before, 2 doesn't work!)
  # typical global configuration:
  ./makenemo -m GCC_PMPC -n AMM7

To compile on *archer* using the intel compiler::

  module unload PrgEnv-cray PrgEnv-gnu
  module load PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  export XIOS_HOME #set to the basepath of XIOS-1 (compile before, 2 doesn't work!)
  #typical global:
  ./makenemo -m XC_ARCHER_INTEL_NOSIGNEDZERO -n ORCA2_LIM_FABM

Shelf-enabled code with FABM coupler is in the ``feat/fabm`` branch.

XIOS-1 library bug
==================

Note that NEMO 3.6 stable is incompatible with XIOS-1 versions more recent than September 2015,
due to what is supposed to be a bug-fix, that is incompatible with NEMO 3.6 (XIOS-1 commit of 1st October 2015).
Use the following repository for a NEMO 3.6 compatible version of XIOS-1:

https://gitlab.ecosystem-modelling.pml.ac.uk/momm/XIOS1/tree/nemo3.6-fix
