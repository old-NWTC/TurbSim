!=======================================================================
MODULE TurbSim_Types

use NWTC_Library

   INTEGER(IntKi), PARAMETER :: MaxMsgLen = 1024 ! Maximum length of error messages

   INTEGER(IntKi), PARAMETER :: SpecModel_NONE   =  0  ! No turbulence
   INTEGER(IntKi), PARAMETER :: SpecModel_IECKAI =  1  ! IEC Kaimal
   INTEGER(IntKi), PARAMETER :: SpecModel_IECVKM =  2  ! IEC von Karman 
   INTEGER(IntKi), PARAMETER :: SpecModel_GP_LLJ =  3  ! Great Plains Low-Level Jet
   INTEGER(IntKi), PARAMETER :: SpecModel_NWTCUP =  4  ! NWTC (upwind)
   INTEGER(IntKi), PARAMETER :: SpecModel_SMOOTH =  5  ! Risoe Smooth-Terrain   
   INTEGER(IntKi), PARAMETER :: SpecModel_WF_UPW =  6  ! Wind Farm Upwind
   INTEGER(IntKi), PARAMETER :: SpecModel_WF_07D =  7  ! Wind Farm  7 rotor diameters downwind
   INTEGER(IntKi), PARAMETER :: SpecModel_WF_14D =  8  ! Wind Farm 14 rotor diameters downwind
   INTEGER(IntKi), PARAMETER :: SpecModel_TIDAL  =  9  ! Tidal (Hydro)
   INTEGER(IntKi), PARAMETER :: SpecModel_RIVER  = 10  ! River (Hydro)
   INTEGER(IntKi), PARAMETER :: SpecModel_API    = 11  ! API
   INTEGER(IntKi), PARAMETER :: SpecModel_MODVKM = 12  ! user-specified scaling in von Karman model
   INTEGER(IntKi), PARAMETER :: SpecModel_USRVKM = 13  ! user-specified scaling in von Karman model
   INTEGER(IntKi), PARAMETER :: SpecModel_USER   = 14  ! User-defined spectra from file
   
   
   type :: RandNum_ParameterType
   
      integer(IntKi)                  :: pRNG
      INTEGER(IntKi)                  :: RandSeed   (3)                           ! The array that holds the initial random seeds for the 3 components.
      INTEGER(IntKi),    ALLOCATABLE  :: RandSeedAry(:)                           ! The array that holds the random seeds.
            
   end type RandNum_ParameterType

   type :: RandNum_OtherStateType
      INTEGER(IntKi), ALLOCATABLE     :: NextSeed   (:)                           ! The array that holds the next random seed for the 3 components.            
   end type RandNum_OtherStateType

   
   type :: TurbSim_GridParameterType
      
      REAL(ReKi)                   :: GridHeight                               ! Grid height
      REAL(ReKi)                   :: GridRes_Y                                ! Distance between two consecutive horizontal points on the grid (Horizontal resolution)
      REAL(ReKi)                   :: GridRes_Z                                ! Distance between two consecutive vertical points on the grid (Vertical resolution)
      REAL(ReKi)                   :: GridWidth                                ! Grid width.
      REAL(ReKi)                   :: Zbottom                                  ! The height of the lowest point on the grid (before tower points are added), equal to Z(1)
      INTEGER(IntKi)               :: HubIndx                                  ! Index that tells where the hub point is in the V matrix
      
      
   end type TurbSim_GridParameterType
   
   
END MODULE TurbSim_Types
   
MODULE TSMods

USE                             NWTC_Library

use TurbSim_Types

IMPLICIT                        NONE
SAVE
type(RandNum_ParameterType)      :: p_RandNum                   ! parameters for random numbers
type(TurbSim_GridParameterType)  :: p_grid                      ! parameters for TurbSim
type(RandNum_OtherStateType)     :: OtherSt_RandNum             ! other states for random numbers (next seed, etc)


REAL(ReKi), PARAMETER        :: profileZmax = 140.                       ! Upper height limit for extrapolating GP_LLJ profiles of ustar and zl
REAL(ReKi), PARAMETER        :: profileZmin =  50.                       ! Lower height limit for extrapolating GP_LLJ profiles of ustar and zl
REAL(ReKi), PARAMETER        :: Omega     = 7.292116E-05                 ! Angular speed of rotation of the earth (rad/s)
REAL(ReKi), PARAMETER        :: Tolerance = 0.0001                       ! The largest difference between two numbers that are assumed to be equal


INTEGER,    PARAMETER        :: IEC_ETM        = 1                       ! Number to indicate the IEC Normal Turbulence Model
INTEGER,    PARAMETER        :: IEC_EWM1       = 2                       ! Number to indicate the IEC Extreme Wind speed Model (50-year)
INTEGER,    PARAMETER        :: IEC_EWM50      = 3                       ! Number to indicate the IEC Extreme Wind speed Model ( 1-year)
INTEGER,    PARAMETER        :: IEC_EWM100     = 5                       ! Number to indicate the IEC Extreme Wind speed Model ( 100-year)
INTEGER,    PARAMETER        :: IEC_NTM        = 4                       ! Number to indicate the IEC Extreme Turbulence Model


INTEGER,    PARAMETER        :: UACT     = 14                            ! I/O unit for AeroDyn coherent turbulence
INTEGER,    PARAMETER        :: UACTTS   = 15                            ! I/O unit for coherent turbulence time step history file
INTEGER,    PARAMETER        :: UAFFW    = 9                             ! I/O unit for AeroDyn FF data (*.bts file).
INTEGER,    PARAMETER        :: UAHH     = 10                            ! I/O unit for AeroDyn HH data (*.hh  file).
INTEGER,    PARAMETER        :: UATWR    = 13                            ! I/O unit for AeroDyn tower data (*.twr file).
INTEGER,    PARAMETER        :: UBFFW    = 16                            ! I/O unit for BLADED FF data (*.wnd file).
INTEGER,    PARAMETER        :: UFFF     = 4                             ! I/O unit for formatted FF data.
INTEGER,    PARAMETER        :: UFTP     = 12                            ! I/O unit for formatted HH turbulence properties.
INTEGER,    PARAMETER        :: UGTP     = 11                            ! I/O unit for GenPro HH turbulence properties.
INTEGER,    PARAMETER        :: UI       = 1                             ! I/O unit for input file.
INTEGER,    PARAMETER        :: US       = 3                             ! I/O unit for summary file.
INTEGER,    PARAMETER        :: USpec    = 17                            ! I/O unit for user-defined spectra

INTEGER,    PARAMETER        :: UC       = 22                            ! I/O unit for Coherence debugging file.
INTEGER,    PARAMETER        :: UD       = 20                            ! I/O unit for debugging data.
INTEGER,    PARAMETER        :: UP       = 21                            ! I/O unit for PSD debugging file.

LOGICAL,    PARAMETER        :: COH_OUT  = .FALSE.                       ! This parameter has been added to replace the NON-STANDARD compiler directive previously used
LOGICAL,    PARAMETER        :: DEBUG_OUT= .FALSE.                       ! This parameter has been added to replace the NON-STANDARD compiler directive previously used
LOGICAL,    PARAMETER        :: PSD_OUT  = .FALSE. !                     ! This parameter has been added to replace the NON-STANDARD compiler directive previously used
LOGICAL,    PARAMETER        :: MVK      = .FALSE.                       ! This parameter has been added to replace the NON-STANDARD compiler directive previously used

LOGICAL,    PARAMETER        :: PeriodicY = .FALSE. !.TRUE.
CHARACTER(1), parameter      ::  Comp (3) = (/ 'u', 'v', 'w' /)  ! The names of the wind components







CHARACTER( 23)               :: IECeditionStr (3) = &   ! BJJ not a parameter because may be using -2 or -3 standards
                                (/'IEC 61400-1 Ed. 1: 1993', &
                                  'IEC 61400-1 Ed. 2: 1999', &
                                  'IEC 61400-1 Ed. 3: 2005'/)            ! The string description of the IEC 61400-1 standard being used



REAL(ReKi)                   :: AnalysisTime                             ! Analysis Time. (amount of time for analysis, allows user to perform analysis using one time length, but output UsableTime
REAL(ReKi)                   :: ChebyCoef_WS(11)                         ! The Chebyshev coefficients for wind speed
REAL(ReKi)                   :: ChebyCoef_WD(11)                         ! The Chebyshev coefficients for wind direction
REAL(ReKi)                   :: COHEXP                                   ! Coherence exponent
REAL(ReKi)                   :: CTLy                                     ! Fractional location of tower centerline from right (looking downwind) to left side of the dataset.
REAL(ReKi)                   :: CTLz                                     ! Fractional location of hub height from the bottom of the dataset.
REAL(ReKi)                   :: CTStartTime                              ! Minimum time to add coherent structures
REAL(ReKi)                   :: DistScl                                  ! Disturbance scale for AeroDyn coherent turbulence events
REAL(ReKi)                   :: ETMc                                     ! The c parameter in IEC ETM, 61400-1, Ed 3. Section 6.3.2.3, Eq. 19.  Variable per last sentence in section 7.4.1
REAL(ReKi), ALLOCATABLE      :: EventLen   (:)                           ! The length of each event stored in EventStart() (non-dimensional time)
REAL(ReKi)                   :: EventTimeStep                            ! The average length of timesteps in output events
REAL(ReKi)                   :: Fc                                       ! Coriolis parameter in units (1/sec)
REAL(ReKi), ALLOCATABLE      :: Freq       (:)                           ! The array of frequencies (NumFreq).
REAL(ReKi), ALLOCATABLE      :: Freq_USR(:)                              ! frequencies for the user-defined spectra
REAL(ReKi)                   :: h                                        ! Boundary layer depth
REAL(ReKi)                   :: HFlowAng                                 ! Horizontal flow angle.
REAL(ReKi)                   :: HH_HFlowAng                              ! Horizontal flow angle at the hub (may be different than HFlowAng if using direction profile).
REAL(ReKi)                   :: HubHt                                    ! Hub height.
REAL(ReKi)                   :: InCDec     (3)                           ! Contains the coherence decrements
REAL(ReKi)                   :: InCohB     (3)                           ! Contains the coherence b/L (offset) parameters
REAL(ReKi)                   :: L                                        ! M-O length
REAL(ReKi), ALLOCATABLE      :: L_USR      (:)                           ! User-specified von Karman length scale, varying with height
REAL(ReKi)                   :: Latitude                                 ! The site latitude in radians
REAL(ReKi)                   :: PerTurbInt                               ! Percent Turbulence Intensity
REAL(ReKi)                   :: PC_UW                                    ! u'w' cross-correlation coefficient
REAL(ReKi)                   :: PC_UV                                    ! u'v' cross-correlation coefficient
REAL(ReKi)                   :: PC_VW                                    ! v'w' cross-correlation coefficient
REAL(ReKi), ALLOCATABLE      :: pkCTKE     (:)                           ! Array containing the peak CTKE of each coherent event
REAL(ReKi)                   :: PLExp                                    ! Rotor disk power law exponent
REAL(ReKi), ALLOCATABLE      :: PhaseAngles (:,:,:)                           ! The array that holds the random phases [number of points, number of frequencies, number of wind components=3].
REAL(ReKi)                   :: RICH_NO                                  ! Gradient Richardson number
REAL(ReKi)                   :: RotorDiameter                            ! The assumed diameter of the rotor
REAL(ReKi), ALLOCATABLE      :: S          (:,:,:)                       ! The turbulence PSD array (NumFreq,NTot,3).
REAL(ReKi), ALLOCATABLE      :: SDary      (:)                           ! The array of standard deviations (NumGrid_Z,NumGrid_Y).
REAL(ReKi)                   :: SigmaIEC                                 ! IEC standard deviation.
REAL(ReKi), ALLOCATABLE      :: Sigma_USR  (:)                           ! User-specified standard deviation of the wind speed components (isotropic), varying with height
REAL(ReKi)                   :: StdScale   (3)                           ! Scaling for the user-specified standard deviation
REAL(ReKi)                   :: TimeStep                                 ! Time step.
REAL(ReKi)                   :: Sigma_U2                                 ! Standard Deviation of U velocity, squared.
REAL(ReKi)                   :: Sigma_V2                                 ! Standard Deviation of V velocity, squared.
REAL(ReKi)                   :: Sigma_W2                                 ! Standard Deviation of W velocity, squared.
REAL(ReKi)                   :: TurbIntH20                               ! Turbulence intensity used for HYDRO module.
REAL(ReKi), ALLOCATABLE      :: TRH        (:)                           ! The transfer function  matrix (NumSteps).
REAL(ReKi)                   :: TsclFact                                 ! Scale factor for time (h/U0) in coherent turbulence events
REAL(ReKi), ALLOCATABLE      :: U          (:)                           ! The steady u-component wind speeds for the grid (ZLim).
REAL(ReKi)                   :: H_ref                                    ! Height for reference wind speed.
REAL(ReKi), ALLOCATABLE      :: DUDZ       (:)                           ! The steady u-component wind shear for the grid (ZLim).
REAL(ReKi), ALLOCATABLE      :: U_USR      (:)                           ! User-specified total wind speed, varying with height
REAL(ReKi)                   :: UHub                                     ! Hub-height (total) wind speed (m/s)
REAL(ReKi)                   :: UJetMax                                  ! The (horizontal) wind speed at the height of the jet maximum (m/s)
REAL(ReKi)                   :: UsableTime                               ! Usable time.  Program adds GridWidth/MeanHHWS.
REAL(ReKi), ALLOCATABLE      :: Uspec_USR(:)                             ! user-defined u-component spectrum
REAL(ReKi)                   :: Ustar                                    ! Shear or friction velocity (m/s) -- rotor-disk average
REAL(ReKi), ALLOCATABLE      :: Ustar_profile(:)                         ! A profile of ustar (measure of friction velocity with height)
REAL(ReKi)                   :: UstarDiab                                ! The diabatic ustar value
REAL(ReKi)                   :: UstarOffset                              ! A scaling/offset value used with the Ustar_profile to ensure that the mean hub u'w' and ustar inputs agree with the profile values
REAL(ReKi)                   :: UstarSlope                               ! A scaling/slope value used with the Ustar_profile to ensure that the mean hub u'w' and ustar inputs agree with the profile values
REAL(ReKi)                   :: U_Ref                                    ! The input wind speed at the reference height.  (Added by M. Buhl for API profiles)
REAL(ReKi), ALLOCATABLE      :: V          (:,:,:)                       ! An array containing the summations of the rows of H (NumSteps,NTot,3).
REAL(ReKi)                   :: VFlowAng                                 ! Vertical flow angle.
REAL(ReKi)                   :: Vave                                     ! The IEC Vave for ETM
REAL(ReKi)                   :: Vref                                     ! The IEC Vref for ETM
REAL(ReKi), ALLOCATABLE      :: Vspec_USR(:)                             ! user-defined v-component spectrum

REAL(ReKi), ALLOCATABLE      :: WindDir_profile(:)                       ! A profile of horizontal wind angle (measure of wind direction with height)
REAL(ReKi), ALLOCATABLE      :: WindDir_USR    (:)                       ! User-specified wind direction profile, varying with height
REAL(ReKi), ALLOCATABLE      :: Work       (:,:)                         ! A temporary work array (NumSteps+2,3).
REAL(ReKi), ALLOCATABLE      :: Wspec_USR(:)                             ! user-defined w-component spectrum
REAL(ReKi), ALLOCATABLE      :: Y          (:)                           ! The lateral locations of the points (YLim).
REAL(ReKi)                   :: Ym_max                                   ! The nondimensional lateral width of the coherent turbulence dataset
REAL(ReKi), ALLOCATABLE      :: Z          (:)                           ! The vertical locations of the points (ZLim).
REAL(ReKi), ALLOCATABLE      :: Z_USR      (:)                           ! Heights of user-specified variables
REAL(ReKi)                   :: Z0                                       ! Surface roughness length, meters
REAL(ReKi)                   :: ZI                                       ! Mixing layer depth
REAL(ReKi)                   :: ZJetMax                                  ! The height of the jet maximum (m)
REAL(ReKi)                   :: ZL                                       ! A measure of stability
REAL(ReKi), ALLOCATABLE      :: ZL_profile(:)                            ! A profile of z/l (measure of stability with height)
REAL(ReKi)                   :: ZLoffset                                 ! An offset to align the zl profile with the mean zl input parameter
REAL(ReKi)                   :: Zm_max                                   ! The nondimensional vertical height of the coherent turbulence dataset

!REAL(ReKi)                   :: RefHt                                    ! Reference height. ADDED BY Y.G.
!REAL(ReKi)                   :: URef                                     ! Wind Speed at Reference Height. ADDED BY Y.G.
 REAL(ReKi)                   :: U0_1HR

INTEGER,    ALLOCATABLE      :: EventName  (:)                           ! The timestep where the event starts, which determines the name of the event file
INTEGER,    ALLOCATABLE      :: EventTS    (:)                           ! The length of each event stored in EventStart() (number of timesteps)
INTEGER                      :: IECedition                               ! The edition number of the IEC 61400-1 standard that is being used (determines the scaling)
INTEGER                      :: IECstandard                              ! The standard number (x) of the IEC 61400-x that is being used
INTEGER                      :: IEC_WindType                             ! Number to indicate the IEC wind type
INTEGER,    ALLOCATABLE      :: IYmax      (:)                           ! A temporary variable holding the maximum number of horizontal positions at each z
INTEGER                      :: MaxDims                                  ! Maximum number of time steps plus 2.
INTEGER                      :: NTot                                     ! Number of points in grid, plus the hub center.
INTEGER                      :: NumCTt                                   ! Number of data points in the output coherent event timestep file
INTEGER                      :: NumEvents                                ! Number of events in the event data file
INTEGER                      :: NumFreq                                  ! Number of frequencies (=NumSteps/2).
INTEGER                      :: NumGrid_Y                                ! Grid dimension. (in horizontal direction)
INTEGER                      :: NumGrid_Z                                ! Grid dimension. (in vertical direction)
INTEGER                      :: NumOutSteps                              ! Number of output time steps.
INTEGER                      :: NumSteps                                 ! Number of time steps for the FFT.
INTEGER                      :: NumUSRf                                  ! Number of frequencies in the user-defined spectra
INTEGER                      :: NumUSRz                                  ! Number of heights defined in the user-defined profiles.
INTEGER                      :: ScaleIEC                                 ! Flag to indicate if turbulence should be scaled to target value; 0 = NO scaling; 1 = scale based on hub; 2 = scale each point individually
INTEGER                      :: SpecModel                                ! Integer value of spectral model (see SpecModel enum)
INTEGER                      :: YLim                                     ! Number of horizontal positions in the grid
INTEGER                      :: ZLim                                     ! Number of vertical positions in the grid, plus extra hub point (if necessary), plus tower points


LOGICAL                      :: Clockwise                                ! Flag to indicate clockwise rotation when looking downwind.
LOGICAL                      :: ExtraHubPT                               ! Flag to indicate if the hub is on the regular grid or if an extra point must be added
LOGICAL                      :: ExtraTwrPT                               ! Flag to indicate if the tower is on the regular grid or if an extra point must be added

LOGICAL                      :: KHtest                                   ! Flag to indicate that turbulence should be extreme, to demonstrate effect of KH billows
LOGICAL                      :: NumTurbInp                               ! Flag to indicate if turbulence is user-specified (as opposed to IEC standard A, B, or C)
LOGICAL                      :: Periodic                                 ! Flag to indicate that output files must contain exactly one full (time) period
LOGICAL                      :: UVskip                                   ! Flag to determine if UV cross-feed term should be skipped or used
LOGICAL                      :: UWskip                                   ! Flag to determine if UW cross-feed term should be skipped or used
LOGICAL                      :: VWskip                                   ! Flag to determine if VW cross-feed term should be skipped or used
LOGICAL                      :: WrACT                                    ! Flag to output AeroDyn coherent turbulence
LOGICAL                      :: WrADFF                                   ! Flag to output AeroDyn FF data (binary).
LOGICAL                      :: WrADHH                                   ! Flag to output AeroDyn HH data (formatted).
LOGICAL                      :: WrADTWR                                  ! Flag to output AeroDyn tower data (binary).
LOGICAL                      :: WrBHHTP                                  ! Flag to output binary HH turbulence parameters.
LOGICAL                      :: WrBLFF                                   ! Flag to output BLADED FF data (binary)
LOGICAL                      :: WrFHHTP                                  ! Flag to output formatted HH turbulence parameters.
LOGICAL                      :: WrFmtFF                                  ! Flag to output formatted FF data (Traditional SNLWIND-3D format).

CHARACTER(200)               :: DescStr                                  ! String used to describe the run (and the first line of the summary file)
CHARACTER(200)               :: FormStr                                  ! String used to store format specifiers.
CHARACTER(200)               :: FormStr1                                 ! String used to store format specifiers.
CHARACTER(200)               :: FormStr2                                 ! String used to store format specifiers.
CHARACTER(  1)               :: IECTurbC                                 ! IEC turbulence characteristic.
CHARACTER(  1)               :: IECTurbE                                 ! IEC Extreme turbulence class.
CHARACTER( 35)               :: IEC_WindDesc                             ! The description of the IEC wind type
CHARACTER(200)               :: InFile = 'TurbSim.inp'                   ! Root name of the I/O files.
CHARACTER(  6)               :: RNG_type                                 ! Type of Random Number Generator to use
CHARACTER(197)               :: RootName                                 ! Root name of the I/O files.
CHARACTER( 50)               :: TMName                                   ! Turbulence model name.
CHARACTER(  6)               :: TurbModel                                ! Turbulence model.
CHARACTER(  3)               :: WindProfileType                          ! The wind profile type


TYPE                         :: Event                                    ! Coherent turbulent event to add to the background wind
   INTEGER                   :: EventNum                                 ! The event number (index into EventName() array)
   REAL(ReKi)                :: TStart                                   ! The time at which to add this event
   REAL(ReKi)                :: delt                                     ! The delta time before the event begins (for interpolation in AeroDyn)
   LOGICAL(1)                :: Connect2Prev = .FALSE.                   ! Whether this event is connected to the next, otherwise there is space between them
   TYPE(Event), POINTER      :: Next         => NULL()                   ! The next event to add
END TYPE

TYPE (Event), POINTER        :: PtrHead      => NULL()                   ! Pointer to the first event
TYPE (Event), POINTER        :: PtrTail      => NULL()                   ! Pointer to the last event

TYPE     :: CohStr_ParameterType   
   REAL(ReKi)              ::  ScaleWid                        ! Scaling width for LE coherent turbulence (RotDiam in AeroDyn FD_Wind)
   REAL(ReKi)              ::  ScaleVel                        ! Scaling velocity for LE coherent turbulence, U0.  2*U0 is the difference in wind speed between the top and bottom of the wave.   
   REAL(ReKi)              ::  Uwave                           ! Wind speed at center of the k-h wave 
   REAL(ReKi)              ::  Wsig                            ! Standard deviation of the w-component wind speed
   
   CHARACTER(200)               :: CTEventPath                              ! String used to store the name of the coherent event definition file
   CHARACTER(200)               :: CTEventFile                              ! String used to store the name of the coherent event definition file
   CHARACTER(  3)               :: CTExt                                    ! String used to determine the type of coherent structures ("dns" or "les")
   
END TYPE CohStr_ParameterType

TYPE :: CohStr_OutputType
   REAL(ReKi)                   :: lambda              ! The expected value of interarrival times for the Poisson process
   INTEGER                      :: NumCTEvents                              ! Number of events to be inserted into the .cts file
   REAL(ReKi)                   :: ExpectedTime        ! Amount of time the coherent structures should take
   REAL(ReKi)                   :: EventTimeSum = 0.0  ! Amount of time the coherent structure takes
END TYPE CohStr_OutputType

TYPE(CohStr_ParameterType)   :: p_CohStr
TYPE(CohStr_OutputType)      :: y_CohStr

END MODULE TSMods
!=======================================================================
