==========================
AMM7_FABM_ERSEM experiment
==========================

This folder contains files for a test run of the AMM7 NEMO-FABM-ERSEM configuration on the UK HPC Archer facility.

FABM-NEMO should be built with the following cmake configuration::

   cmake ~/git/FABM/src -DFABM_HOST=nemo -DFABM_ERSEM_BASE=~/git/edge -DFABM_EMBED_VERSION=ON -DCMAKE_Fortran_FLAGS:STRING=-O3 -fp-model source -traceback

