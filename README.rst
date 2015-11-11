===================================================================
Fork of NEMO 3.6 stable branch with modifcations for shelf systems.
===================================================================

Cloned with:

::

   git svn clone --username momme http://forge.ipsl.jussieu.fr/nemo/svn/branches/2015/nemo_v3_6_STABLE .

Shelf-enabled code is in the ``feat/shelf-enabled`` branch, currently based on svn revisioni ``r5569``.

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
  module load cray-netcdf-hdfparallel
  export XIOS_HOME #set to the basepath of XIOS-1 (compile before, 2 doesn't work!)
  #typical global:
  ./makenemo -m XC_ARCHER_INTEL_NOSIGNEDZERO -n ORCA2_LIM_FABM

