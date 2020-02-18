MODULE trcnam_fabm
   !!======================================================================
   !!                      ***  MODULE trcnam_fabm  ***
   !! TOP :   initialisation of some run parameters for FABM bio-model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_fabm'   :                                       FABM model
   !!----------------------------------------------------------------------
   !! trc_nam_fabm      : FABM initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables

   USE par_fabm
   USE trcsms_fabm


   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_fabm   ! called by trcnam.F90 module
   PUBLIC   trc_nam_fabm_override ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_fabm
   END SUBROUTINE trc_nam_fabm

   SUBROUTINE trc_nam_fabm_override(sn_tracer)
      TYPE(PTRACER), DIMENSION(jpmaxtrc), INTENT(INOUT) :: sn_tracer

      INTEGER :: jn
      CHARACTER(LEN=3) :: index

      DO jn=1,jp_fabm
         IF (sn_tracer(jn)%clsname /= 'NONAME' .AND. sn_tracer(jn)%clsname /= model%state_variables(jn)%name) THEN
            WRITE (index,'(i0)') jn
            CALL ctl_stop('Tracer name mismatch in namtrc: '//TRIM(sn_tracer(jn)%clsname)//' found at sn_tracer('//TRIM(index)//') where '//TRIM(model%state_variables(jn)%name)//' was expected.')
         END IF
         sn_tracer(jn)%clsname = TRIM(model%state_variables(jn)%name)
         sn_tracer(jn)%cllname = TRIM(model%state_variables(jn)%long_name)
         sn_tracer(jn)%clunit = TRIM(model%state_variables(jn)%units)
         sn_tracer(jn)%llinit = .FALSE.
      END DO
   END SUBROUTINE trc_nam_fabm_override

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                             No FABM
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_fabm                      ! Empty routine
   END  SUBROUTINE  trc_nam_fabm

   SUBROUTINE trc_nam_fabm_override
   END SUBROUTINE trc_nam_fabm_override
#endif  

   !!======================================================================
END MODULE trcnam_fabm
