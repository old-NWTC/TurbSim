module ts_errors
   use NWTC_Library
contains

SUBROUTINE TS_Abort ( Message )
   ! This routine outputs fatal warning messages and ends the program.

   USE                              TSMods, ONLY:   US                           ! unit number for summary file

      ! Argument declarations.

   CHARACTER(*), INTENT(IN)      :: Message                                      ! Warning message.


      ! Write the message to the summary file

   WRITE (US, "(/'ERROR:  ', A / )") Message
   WRITE (US, "('ABORTING PROGRAM.')" )

      ! Write the message to the screen
   CALL ProgAbort ( Message )

RETURN

END SUBROUTINE TS_Abort
!=======================================================================
SUBROUTINE TS_Warn ( Message, UnSum )

   ! This routine outputs non-fatal warning messages and returns to the calling routine.

   USE                              TSMods, ONLY:   US                           ! unit number for summary file

      ! Argument declarations.

   CHARACTER(*),   INTENT(IN)      :: Message                                      ! Warning message.
   INTEGER(IntKi), INTENT(IN)      :: UnSum                                        ! Summary unit number: if > 0, will print a message in the sum file.

      ! Write the message to the screen

   CALL WrScr( '' )
   CALL ProgWarn ( Message )
   CALL WrScr( '' )


      ! Write the message to the summary file if requested

   IF ( UnSum > 0 ) THEN
      WRITE (UnSum, "(/'WARNING:  ', A / )") Message
   ENDIF


RETURN
END SUBROUTINE TS_Warn ! ( Message )
!=======================================================================

end module ts_errors