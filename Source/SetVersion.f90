!=======================================================================
SUBROUTINE SetVersion


!  This routine sets the version number.  By doing it this way instead
!  of the old way of setting it in a DATA statement in Modules.f90, we
!  will no longer have to recompile everything every time we change
!  versions.

USE               NWTC_Library, ONLY : ProgName
USE               NWTC_Library, ONLY : ProgVer

IMPLICIT          NONE

ProgName = 'TurbSim'                     ! The name of this program.
ProgVer  = ' (v1.05.01c, 24-Feb-2012)'
!BJJ: make sure that DEBUG_OUT and PSD_OUT and COH_OUT are .FALSE.

RETURN
END SUBROUTINE SetVersion
