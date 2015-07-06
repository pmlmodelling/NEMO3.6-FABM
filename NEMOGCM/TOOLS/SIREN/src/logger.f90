!----------------------------------------------------------------------
! NEMO system team, System and Interface for oceanic RElocable Nesting
!----------------------------------------------------------------------
!
! MODULE: logger
!
! DESCRIPTION:
!> @brief This module create logger file and allow to fill it depending of verbosity.
!> @details
!> verbosity could be choosen between :
!>    - trace : Most detailed information.
!>    - debug : Detailed information on the flow through the system.
!>    - info  : Interesting runtime events (startup/shutdown).
!>    - warning: Use of deprecated APIs, poor use of API, 'almost' errors, 
!> other runtime situations that are undesirable or unexpected, 
!> but not necessarily "wrong". 
!>    - error : Other runtime errors or unexpected conditions.
!>    - fatal : Severe errors that cause premature termination.<br />
!>  default verbosity is warning
!
!> If total number of error exceeded maximum number 
!> authorized, program stop.
!>
!> to open/create logger file:<br/>
!> @code
!>    CALL logger_open(cd_file, [cd_verbosity,] [id_loggerid,] [id_maxerror])
!> @endcode
!> - cd_file is logger file name
!> - cd_verbosity is verbosity to be used [optional, default 'warning']
!> - id_loggerid is file id [optional, use only to flush]
!> - id_maxerror is the maximum number of error authorized before program stop [optional, default 5]
!>
!> to close logger file:<br/>
!> @code
!> CALL logger_close()
!> @endcode
!>
!> to write header in logger file:<br/>
!> @code
!> CALL logger_header()
!> @endcode
!>
!> to write footer in logger file:<br/>
!> @code
!> CALL logger_footer()
!> @endcode
!>
!> to flushing output:<br/>
!> @code
!> CALL logger_flush()
!> @endcode
!>
!> to write TRACE message in logger file:<br/>
!> @code
!> CALL logger_trace(cd_msg [,ld_flush])
!> @endcode
!>    - cd_msg is TRACE message
!>    - ld_flush to flush output [optional]
!>
!> to write DEBUG message in logger file:<br/>
!> @code
!> CALL logger_debug(cd_msg [,ld_flush])
!> @endcode
!>    - cd_msg is DEBUG message
!>    - ld_flush to flush output [optional]
!>
!> to write INFO message in logger file:<br/>
!> @code
!> CALL logger_info(cd_msg [,ld_flush])
!> @endcode
!>    - cd_msg is INFO message
!>    - ld_flush to flush output [optional]
!>
!> to write WARNING message in logger file:<br/>
!> @code
!> CALL logger_warn(cd_msg [,ld_flush])
!> @endcode
!>    - cd_msg is WARNING message
!>    - ld_flush to flush output [optional]
!>
!> to write ERROR message in logger file:<br/>
!> @code
!> CALL logger_error(cd_msg [,ld_flush])
!> @endcode
!>    - cd_msg is ERROR message
!>    - ld_flush to flush output [optional]
!>
!> to write FATAL message in logger file:<br/>
!> @code
!> CALL logger_fatal(cd_msg)
!> @endcode
!>    - cd_msg is FATAL message
!>
!> Examples :<br />
!> @code
!>   CALL logger_open('loggerfile.txt','info')
!>
!>   CALL logger_header()
!>   CALL logger_debug('une info de debug')
!>   CALL logger_info('une info')
!>   CALL logger_warn('un warning')
!>   CALL logger_error('une erreur')
!>   CALL logger_footer()
!>   CALL logger_close()
!> @endcode
!>
!> @code
!>   CALL logger_open('loggerfile.txt')
!>
!>   CALL logger_header()
!>   CALL logger_debug('une info de debug')
!>   CALL logger_info('une info')
!>   CALL logger_warn('un warning')
!>   CALL logger_error('une erreur')
!>   CALL logger_footer()
!>   CALL logger_close()
!> @endcode
!
!> @author
!> J.Paul
! REVISION HISTORY:
!> @date November, 2013- Initial Version
!>
!> @note Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
!----------------------------------------------------------------------
MODULE logger
   USE kind                            ! F90 kind parameter
   USE fct                             ! basic useful function
   IMPLICIT NONE
   ! NOTE_avoid_public_variables_if_possible

   ! type and variable
   PRIVATE :: TLOGGER            !< logger structure

   PRIVATE :: tm_logger          !< logger structure
   PRIVATE :: im_nverbosity      !< number of log level
   PRIVATE :: cm_verbosity       !< verbosity array

   ! function and subroutine
   PUBLIC :: logger_open        !< create a log file with given verbosity
   PUBLIC :: logger_close       !< close log file
   PUBLIC :: logger_header      !< write header on log file
   PUBLIC :: logger_footer      !< write footer on log file
   PUBLIC :: logger_flush       !< flushing output
   PUBLIC :: logger_trace       !< write trace    message in log file
   PUBLIC :: logger_debug       !< write debug    message in log file 
   PUBLIC :: logger_info        !< write info     message in log file
   PUBLIC :: logger_warn        !< write warning  message in log file
   PUBLIC :: logger_error       !< write error    message in log file
   PUBLIC :: logger_fatal       !< write fatal    message in log file, and stop

   PRIVATE :: logger__write     ! cut message to get maximum of 80 character by line in log file

   TYPE TLOGGER   !< logger structure
      INTEGER(i4)       :: i_id = 0                 !< log file id
      CHARACTER(LEN=lc) :: c_name                   !< log file name
      CHARACTER(LEN=lc) :: c_verbosity = "warning"  !< verbosity choose
      CHARACTER(LEN=lc) :: c_verb = ""              !< array of "verbosities" to used 
      INTEGER(i4)       :: i_nerror   = 0           !< number of error
      INTEGER(i4)       :: i_nfatal   = 0           !< number of fatal error
      INTEGER(i4)       :: i_maxerror = 5           !< maximum number of error before stoping program
   END TYPE TLOGGER   

   !  module variable
   INTEGER(i4), PARAMETER :: im_nverbosity=6     !< number of log level
   CHARACTER(len=*), DIMENSION(im_nverbosity), PARAMETER :: cm_verbosity= & !< verbosity array 
   &               (/ 'trace   ',&
   &                  'debug   ',&
   &                  'info    ',& 
   &                  'warning ',&
   &                  'error   ',&
   &                  'fatal   '/)

   TYPE(TLOGGER), SAVE :: tm_logger      !< logger structure
                                                 
CONTAINS
   !-------------------------------------------------------------------
   !> @brief This subroutine create a log file with default verbosity
   !> ('warning').
   !> @details
   !> Optionally verbosity could be change to 
   !> ('trace','debug','info',warning','error','fatal').<br/>
   !> Optionally maximum number of error allowed could be change.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] cd_file      log file name
   !> @param[in] cd_verbosity log file verbosity
   !> @param[in] id_logid     log file id (use to flush)
   !> @param[in] id_maxerror  maximum number of error
   !-------------------------------------------------------------------
   SUBROUTINE logger_open(cd_file, cd_verbosity, id_logid, id_maxerror)
      IMPLICIT NONE
      ! Argument
      CHARACTER(len=*), INTENT(IN) :: cd_file                ! log file name
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: cd_verbosity ! log file verbosity
      INTEGER(i4),      INTENT(IN), OPTIONAL :: id_logid     ! log file id
      INTEGER(i4),      INTENT(IN), OPTIONAL :: id_maxerror  ! log max error

      ! local variable
      INTEGER(i4) :: il_status

      ! loop
      INTEGER(i4) :: ji
      !----------------------------------------------------------------
      ! get id if not already define
      IF( PRESENT(id_logid) )THEN
         tm_logger%i_id=id_logid
      ELSE
         tm_logger%i_id=fct_getunit()
      ENDIF

      ! open log file
      OPEN( tm_logger%i_id, &
      &     STATUS="unknown",    &
      &     FILE=TRIM(ADJUSTL(cd_file)),  &
      &     ACTION="write",      &
      &     POSITION="append",   &
      &     IOSTAT=il_status)
      CALL fct_err(il_status)

      ! keep filename
      tm_logger%c_name=TRIM(ADJUSTL(cd_file))

      ! if present, change verbosity value
      IF( PRESENT(cd_verbosity) )THEN
         tm_logger%c_verbosity=TRIM(ADJUSTL(cd_verbosity))
      ENDIF

      ! compute "tab" of verbosity to be used
      IF( TRIM(ADJUSTL(tm_logger%c_verb)) == "" )THEN
         DO ji=im_nverbosity,1,-1
            tm_logger%c_verb = &
            &  TRIM(tm_logger%c_verb)//" "//TRIM(ADJUSTL(cm_verbosity(ji)))
            IF( TRIM(tm_logger%c_verbosity) == TRIM(cm_verbosity(ji)) )THEN
               EXIT
            ENDIF
         ENDDO
      ENDIF

      IF( PRESENT(id_maxerror) )THEN
         tm_logger%i_maxerror=id_maxerror
      ENDIF

   END SUBROUTINE logger_open
   !-------------------------------------------------------------------
   !> @brief This subroutine close a log file.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !-------------------------------------------------------------------
   SUBROUTINE logger_close()
      IMPLICIT NONE
      ! local variable
      INTEGER(i4) :: il_status
      !----------------------------------------------------------------
      IF( tm_logger%i_id /= 0 )THEN
         tm_logger%i_id = 0
         CLOSE( tm_logger%i_id, &
         &      IOSTAT=il_status)      
         CALL fct_err(il_status)
      ELSE
          CALL logger_open('logger.log')
          CALL logger_header()
          CALL logger_fatal('you must have create logger to use logger_close')
      ENDIF

   END SUBROUTINE logger_close
   !-------------------------------------------------------------------
   !> @brief This subroutine flushing output into log file.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !-------------------------------------------------------------------
   SUBROUTINE logger_flush()
      IMPLICIT NONE
      !----------------------------------------------------------------
      IF( tm_logger%i_id /= 0 )THEN
         CALL logger_close()
         CALL logger_open( tm_logger%c_name, tm_logger%c_verbosity, tm_logger%i_id, &
         &              tm_logger%i_maxerror )     
      ELSE
          CALL logger_open('logger.log')
          CALL logger_header()
          CALL logger_fatal('you must have create logger to use logger_flush')
      ENDIF

   END SUBROUTINE logger_flush
   !-------------------------------------------------------------------
   !> @brief This subroutine write header on log file.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !-------------------------------------------------------------------
   RECURSIVE SUBROUTINE logger_header()
      IMPLICIT NONE
      ! local variable
      INTEGER(i4)       :: il_status
      !----------------------------------------------------------------
      IF( tm_logger%i_id /= 0 )THEN
         WRITE( tm_logger%i_id,    &
            &   FMT='(4(a/))',     &
            &   IOSTAT=il_status ) &
            &   "--------------------------------------------------",&
            &   "INIT     : verbosity "//TRIM(tm_logger%c_verbosity),&
            &   "INIT     : max error "//TRIM(fct_str(tm_logger%i_maxerror)), &
            &   "--------------------------------------------------"
         CALL fct_err(il_status)
      ELSE
          CALL logger_open('logger.log')
          CALL logger_header()
          CALL logger_fatal('you must have create logger to use logger_header')
      ENDIF

   END SUBROUTINE logger_header
   !-------------------------------------------------------------------
   !> @brief This subroutine write footer on log file.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !-------------------------------------------------------------------
   SUBROUTINE logger_footer()
      IMPLICIT NONE
      ! local variable
      INTEGER(i4)       :: il_status
      !----------------------------------------------------------------
      IF( tm_logger%i_id /= 0 )THEN
         WRITE( tm_logger%i_id,    &
            &   FMT='(4(/a))',     &
            &   IOSTAT=il_status ) &
            &   "--------------------------------------------------",&
            &   "END      : log ended ",              &
            &   "END      : "//TRIM(fct_str(tm_logger%i_nerror))//   &
            &   " ERROR detected ",                                  &
            &   "END      : "//TRIM(fct_str(tm_logger%i_nfatal))//   &
            &   " FATAL detected ",                                  &
            &   "--------------------------------------------------"
         CALL fct_err(il_status)
      ELSE
          CALL logger_open('logger.log')
          CALL logger_header()
          CALL logger_fatal('you must have create logger to use logger_footer')
      ENDIF
   END SUBROUTINE logger_footer
   !-------------------------------------------------------------------
   !> @brief This subroutine write trace message on log file.
   !> @details
   !> Optionally you could flush output.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] cd_msg    message to write
   !> @param[in] ld_flush  flushing ouput
   !-------------------------------------------------------------------
   SUBROUTINE logger_trace(cd_msg, ld_flush)
      IMPLICIT NONE
      ! Argument
      CHARACTER(LEN=*), INTENT(IN)  :: cd_msg
      LOGICAL,          INTENT(IN), OPTIONAL :: ld_flush
      !----------------------------------------------------------------
      IF( tm_logger%i_id /= 0 )THEN
         IF( INDEX(TRIM(tm_logger%c_verb),'trace')/=0 )THEN

            CALL logger__write("TRACE   :",cd_msg)

            IF( PRESENT(ld_flush) )THEN
               IF( ld_flush )THEN
                  CALL logger_flush()
               ENDIF
            ENDIF      
         ENDIF
      ELSE
          CALL logger_open('logger.log')
          CALL logger_header()
          CALL logger_fatal('you must have create logger to use logger_trace')
      ENDIF
   END SUBROUTINE logger_trace
   !-------------------------------------------------------------------
   !> @brief This subroutine write debug message on log file.
   !> @details
   !> Optionally you could flush output.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] cd_msg    message to write
   !> @param[in] ld_flush  flushing ouput
   !-------------------------------------------------------------------
   SUBROUTINE logger_debug(cd_msg, ld_flush)
      IMPLICIT NONE
      ! Argument
      CHARACTER(LEN=*), INTENT(IN)  :: cd_msg
      LOGICAL,          INTENT(IN), OPTIONAL :: ld_flush
      !----------------------------------------------------------------
      IF( tm_logger%i_id /= 0 )THEN
         IF( INDEX(TRIM(tm_logger%c_verb),'debug')/=0 )THEN

            CALL logger__write("DEBUG   :",cd_msg)

            IF( PRESENT(ld_flush) )THEN
               IF( ld_flush )THEN
                  CALL logger_flush()
               ENDIF
            ENDIF      
         ENDIF
      ELSE
          CALL logger_open('logger.log')
          CALL logger_header()
          CALL logger_fatal('you must have create logger to use logger_debug')
      ENDIF
   END SUBROUTINE logger_debug
   !-------------------------------------------------------------------
   !> @brief This subroutine write info message on log file.
   !> @details
   !> Optionally you could flush output.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] cd_msg    message to write
   !> @param[in] ld_flush  flushing ouput
   !-------------------------------------------------------------------
   SUBROUTINE logger_info(cd_msg, ld_flush)
      IMPLICIT NONE
      ! Argument
      CHARACTER(LEN=*), INTENT(IN)  :: cd_msg
      LOGICAL,          INTENT(IN), OPTIONAL :: ld_flush
      !----------------------------------------------------------------
      IF( tm_logger%i_id /= 0 )THEN
         IF( INDEX(TRIM(tm_logger%c_verb),'info')/=0 )THEN

            CALL logger__write("INFO    :",cd_msg)

            IF( PRESENT(ld_flush) )THEN
               IF( ld_flush )THEN
                  CALL logger_flush()
               ENDIF
            ENDIF      
         ENDIF
      ELSE
          CALL logger_open('logger.log')
          CALL logger_header()
          CALL logger_fatal('you must have create logger to use logger_info')
      ENDIF
   END SUBROUTINE logger_info
   !-------------------------------------------------------------------
   !> @brief This subroutine write warning message on log file.
   !> @details
   !> Optionally you could flush output.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] cd_msg    message to write
   !> @param[in] ld_flush  flushing ouput
   !-------------------------------------------------------------------
   SUBROUTINE logger_warn(cd_msg, ld_flush)
      IMPLICIT NONE
      ! Argument
      CHARACTER(LEN=*), INTENT(IN)  :: cd_msg
      LOGICAL,          INTENT(IN), OPTIONAL :: ld_flush
      !----------------------------------------------------------------
      IF( tm_logger%i_id /= 0 )THEN
         IF( INDEX(TRIM(tm_logger%c_verb),'warn')/=0 )THEN

            CALL logger__write("WARNING :",cd_msg)

            IF( PRESENT(ld_flush) )THEN
               IF( ld_flush )THEN
                  CALL logger_flush()
               ENDIF
            ENDIF      
         ENDIF
      ELSE
          CALL logger_open('logger.log')
          CALL logger_header()
          CALL logger_fatal('you must have create logger to use logger_warn')
      ENDIF
   END SUBROUTINE logger_warn
   !-------------------------------------------------------------------
   !> @brief This subroutine write error message on log file.
   !> @details
   !> Optionally you could flush output.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] cd_msg    message to write
   !> @param[in] ld_flush  flushing ouput
   !-------------------------------------------------------------------
   SUBROUTINE logger_error(cd_msg, ld_flush)
      IMPLICIT NONE
      ! Argument
      CHARACTER(LEN=*), INTENT(IN)  :: cd_msg
      LOGICAL,          INTENT(IN), OPTIONAL :: ld_flush

      ! local variable
      CHARACTER(LEN=lc) :: cl_nerror
      !----------------------------------------------------------------
      IF( tm_logger%i_id /= 0 )THEN
         ! increment the error number
         tm_logger%i_nerror=tm_logger%i_nerror+1

         IF( INDEX(TRIM(tm_logger%c_verb),'error')/=0 )THEN

            CALL logger__write("ERROR   :",cd_msg)

            IF( PRESENT(ld_flush) )THEN
               IF( ld_flush )THEN
                  CALL logger_flush()
               ENDIF
            ENDIF      
         ENDIF

         IF( tm_logger%i_nerror >= tm_logger%i_maxerror )THEN
            WRITE(cl_nerror,*) tm_logger%i_maxerror
            CALL logger_fatal(&
            &  'Error count reached limit of '//TRIM(ADJUSTL(cl_nerror)) )
         ENDIF
      ELSE
          CALL logger_open('logger.log')
          CALL logger_header()
          CALL logger_fatal('you must have create logger to use logger_error')
      ENDIF

   END SUBROUTINE logger_error
   !-------------------------------------------------------------------
   !> @brief This subroutine write fatal error message on log file, 
   !> close log file and stop process.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] cd_msg message to write
   !-------------------------------------------------------------------
   RECURSIVE SUBROUTINE logger_fatal(cd_msg)
      IMPLICIT NONE
      ! Argument
      CHARACTER(LEN=*),           INTENT(IN) :: cd_msg
      !----------------------------------------------------------------
      IF( tm_logger%i_id /= 0 )THEN
         IF( INDEX(TRIM(tm_logger%c_verb),'fatal')/=0 )THEN
            ! increment the error number
            tm_logger%i_nfatal=tm_logger%i_nfatal+1

            CALL logger__write("FATAL   :",cd_msg)

            CALL logger_footer()
            CALL logger_close()

            WRITE(*,*) 'FATAL ERROR'
            STOP
         ENDIF
      ELSE
          CALL logger_open('logger.log')
          CALL logger_header()
          CALL logger_fatal('you must have create logger to use logger_fatal')
      ENDIF
   END SUBROUTINE logger_fatal
   !-------------------------------------------------------------------
   !> @brief This subroutine cut message to get maximum of 80 character 
   !> by line in log file.
   !>
   !> @author J.Paul
   !> - November, 2013- Initial Version
   !
   !> @param[in] cd_verb   verbosity of the message to write
   !> @param[in] cd_msg    message to write
   !-------------------------------------------------------------------
   SUBROUTINE logger__write(cd_verb, cd_msg)
      IMPLICIT NONE
      ! Argument
      CHARACTER(LEN=*),           INTENT(IN) :: cd_verb
      CHARACTER(LEN=*),           INTENT(IN) :: cd_msg

      ! local variable
      INTEGER(i4)       :: il_status
      INTEGER(i4)       :: il_verb
      INTEGER(i4)       :: il_msg
      CHARACTER(LEN=lc) :: cl_verb
      CHARACTER(LEN=lc) :: cl_msg
      CHARACTER(LEN=lc) :: cl_tmp

      !----------------------------------------------------------------
      cl_verb=TRIM(ADJUSTL(cd_verb))
      cl_msg=TRIM(ADJUSTL(cd_msg))

      il_verb=LEN_TRIM(cl_verb)
      il_msg=LEN_TRIM(cl_msg)
      DO WHILE( il_verb + il_msg > 78 )
         cl_tmp=TRIM(cl_verb)//' '//TRIM(cl_msg(1:78-il_verb))

         WRITE( tm_logger%i_id,  &
         &      FMT=*,           &
         &      IOSTAT=il_status &
         &      ) TRIM(cl_tmp)
         CALL fct_err(il_status)


         cl_msg=cl_msg(78-il_verb+1:il_msg)
         cl_verb="        :"

         il_msg=LEN_TRIM(cl_msg)

      ENDDO

      cl_tmp=TRIM(cl_verb)//' '//TRIM(cl_msg)
      WRITE( tm_logger%i_id,  &
      &      FMT=*,           &
      &      IOSTAT=il_status &
      &      ) TRIM(cl_tmp)
      CALL fct_err(il_status)

   END SUBROUTINE logger__write
END MODULE logger

