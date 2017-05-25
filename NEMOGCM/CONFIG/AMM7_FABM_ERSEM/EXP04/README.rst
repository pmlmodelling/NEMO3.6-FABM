==========================
AMM7_FABM_ERSEM experiment
==========================

This folder contains files for running a hindcast of the AMM7 NEMO-FABM-ERSEM configuration on the UK HPC Archer facility using v1 of FABM enabled ERSEM from the Shelf Seas Biogeochemistry Program with variable atmospheric pCO2.

FABM-NEMO should be built with the following cmake configuration::

   cmake ~/git/FABM/src -DFABM_HOST=nemo -DFABM_ERSEM_BASE=~/git/edge -DFABM_EMBED_VERSION=ON -DCMAKE_Fortran_FLAGS:STRING=-O3 -fp-model source -traceback

