===================================================================
Fork of NEMO 3.6 stable branch with modifcations for shelf systems.
===================================================================

Cloned with:

::

   git svn clone --username momme http://forge.ipsl.jussieu.fr/nemo/svn/branches/2015/nemo_v3_6_STABLE .

The `master` branch is for maintenance with the Paris trunk only.
Shelf-enabled code is in the ``feat/shelf-enabled`` branch, currently based on svn revision ``r5569``.
Shelf-enabled code with FABM coupler is in the ``feat/fabm`` branch.

XIOS-1 library bug
==================

Note that NEMO 3.6 stable is incompatible with XIOS-1 versions more recent than September 2015,
due to what is supposed to be a bug-fix, that is incompatible with NEMO 3.6 (XIOS-1 commit of 1st October 2015).
Use the following repository for a NEMO 3.6 compatible version of XIOS-1:

https://gitlab.ecosystem-modelling.pml.ac.uk/momm/XIOS1/tree/nemo3.6-fix
