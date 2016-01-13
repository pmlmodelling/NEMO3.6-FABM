============================================================================
Fork of NEMO 3.6 stable branch with modifcations for FABM and shelf systems.
============================================================================

You are on the shelf-enabled FABM branch `feat/fabm`.

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

Specific notes on the **FABM** implementation
=============================================

FABM needs to be compiled separately before the compilation of NEMO.
Usually, the following suffices for this::

   mkdir -p ~/build/nemo && cd ~/build/nemo
   cmake <FABMDIR>/src/ -DFABM_HOST=nemo -DFABM_ERSEM_BASE=<ERSEMDIR>
   make install

In the above, replace `<FABMDIR>` with the directory with the FABM source code, e.g., `~/fabm-git` and `<ERSEMMDIR>` with the directory with the ERSEM source code, e.g., `~/ersem-git`.

This will create the library in the standard folder `~/local/fabm/nemo/lib` where NEMO-FABM will look for linking to NEMO.

The FABM_ coupler for **NEMO** is added in a sub-folder ``FABM`` in ``NEMOGCM/NEMO/TOP_SRC``.
Changes to existing code in order to accomadate FABM_ within **NEMO** are restricted to the ``NEMOGCM/NEMO/TOP_SRC`` and ``NEMOGCM/TOOLS/COMPILE`` folder and shall be marked in the code in the following way:

Additions are encapsulated using the tags::

   ! +++>>> FABM
   ... Some new code...
   ! FABM <<<+++

Removed sections are encapsulated as::

   ! --->>> FABM
   ! ...Some commented code...
   ! FABM <<<---

(In the FCM scripts the ``!`` is replaced by ``#``.)

The initial FABM_ implementation in this branch is carried over from the NEMO-FABM_ repository developed by M. Butenschön and J. Bruggeman, as a single patch commit.

.. _FABM: http://fabm.net
.. _NEMO-FABM: https://gitlab.ecosystem-modelling.pml.ac.uk/momm/NEMO-FABM

XIOS-1 library bug
==================

Note that NEMO 3.6 stable is incompatible with XIOS-1 versions more recent than September 2015,
due to what is supposed to be a bug-fix, that is incompatible with NEMO 3.6 (XIOS-1 commit of 1st October 2015).
Use the following repository for a NEMO 3.6 compatible version of XIOS-1:

https://gitlab.ecosystem-modelling.pml.ac.uk/momm/XIOS1/tree/nemo3.6-fix
