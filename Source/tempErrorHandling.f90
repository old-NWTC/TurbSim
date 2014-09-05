module ts_errors
   use NWTC_Library
contains

!=======================================================================
SUBROUTINE TS_Abort ( Message )
   ! This routine outputs fatal warning messages and ends the program.

      ! Argument declarations.

   CHARACTER(*), INTENT(IN)      :: Message                                      ! Warning message.


      ! Write the message to the summary file

   !WRITE (p%US, "(/'ERROR:  ', A / )") Message
   !WRITE (p%US, "('ABORTING PROGRAM.')" )

      ! Write the message to the screen
   CALL ProgAbort ( Message, .FALSE., 5.0_ReKi )

RETURN

END SUBROUTINE TS_Abort
!=======================================================================

end module ts_errors