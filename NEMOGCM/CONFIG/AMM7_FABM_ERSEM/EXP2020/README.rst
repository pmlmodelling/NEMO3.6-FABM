==========================
AMM7_FABM_ERSEM experiment
==========================

This folder contains files for running a hindcast of the AMM7 NEMO-FABM-ERSEM configuration on the UK HPC Archer facility using FABM1 enabled ERSEM edge repository.
This is an update of EXP05, to account for new code (mostly FABM1), new inputs and changes in configurations.

FABM-NEMO should be built with the following cmake configuration::

   cmake ~/git/FABM/src -DFABM_HOST=nemo -DFABM_ERSEM_BASE=~/git/edge -DFABM_EMBED_VERSION=ON -DCMAKE_Fortran_FLAGS:STRING=-O3 -fp-model source -traceback

