#================================================================================#
# This makefile created by B. Jonkman on 6-Apr-2015,                             #
# adapted from FAST v8.10.00.                                                    #
# (c) 2015 National Renewable Energy Laboratory                                  #
#                                                                                #
# This makefile has been tested on Windows 7 with gfortran.                      #
# This makefile works with mingw32-make.exe.                                     #
#                                                                                #
# It was designed to be used with:                                               #
#     TurbSim (v2.00.01)                                                         #
#                                                                                #
# Older versions of the source code may not work with this makefile.             #
#================================================================================#

   # 32-bit or 64-bit?
BITS = 32
#BITS = 64


   # Location of source files for TurbSim and its dependencies.

TurbSim_DIR  = ../Source

NWTC_Lib_DIR = $(TurbSim_DIR)/../NWTC_Library/source
NETLIB_DIR   = $(NWTC_Lib_DIR)/NetLib/fftpack \
               $(NWTC_Lib_DIR)/NetLib/lapack
pRNG_DIR     = $(TurbSim_DIR)/dependencies/RNG

BIN_DIR      = ../bin

   # Instructions for linking with LAPACK.
   # You will probably need to change these for your system.

ifeq ($(OS),Windows_NT)
   LAPACK_LINK  = -llapack -lblas -LC:/LAPACK/win32
else
   LAPACK_LINK  = -llapack -lblas
endif


   # Name of compiler to use and flags to use.
FC      = gfortran

FFLAGS  = -O2 -m$(BITS) -fbacktrace -ffree-line-length-none -x f95-cpp-input
LDFLAGS = -O2 -m$(BITS) -fbacktrace -Wl,--large-address-aware

#FFLAGS  = -O0 -m$(BITS) -fbacktrace -ffree-line-length-none -x f95-cpp-input -g -pg
#LDFLAGS = -O0 -m$(BITS) -fbacktrace -Wl,--stack=999999999,--large-address-aware -g -pg

#-DDOUBLE_PRECISION

# -mthreads
# some useful gfortran options:
#  -DFPE_TRAP_ENABLED
#  -Wl,--large-address-aware   # not necessary when $(BITS) is 64 (64-bit target)
#  -Wl,--stack=999999999       # not necessary when $(BITS) is 64
#  -Wconversion-extra -Wconversion
#  -fdefault-real-8 -fcheck=bounds,do,mem,pointer -std=f2003 -O0 -v -Wall
#  -pg                          # generate debugging info for debugger
# http://gcc.gnu.org/onlinedocs/gfortran/Option-Index.html#Option-Index
#
# makefile:
# --warn-undefined-variables



   # Destination and RootName for executable

OUTPUT_NAME = TurbSim
DEST_DIR    = $(BIN_DIR)

   #==========================================================#
   # You should not need to change anything beyond this point #
   #==========================================================#

   # System-specific settings.

ifeq ($(OS),Windows_NT)
      # Windows
   DEL_CMD   = del
   EXE_EXT   = _gwin$(BITS).exe
   INTER_DIR = Obj_win$(BITS)
   MD_CMD    = @mkdir
   OBJ_EXT   = .obj
   PATH_SEP  = \\
   SYS_FILE  = SysGnuWin
else
      # Linux
   DEL_CMD   = rm -f
   EXE_EXT   = _glin$(BITS)
   INTER_DIR = Obj_lin$(BITS)
   MD_CMD    = @mkdir -p
   OBJ_EXT   = .o
   PATH_SEP  = /
   SYS_FILE  = SysGnuLinux
endif


   # Source files (by module)

LIB_SOURCES =             \
	SingPrec.f90           \
	NWTC_Base.f90          \
	$(SYS_FILE).f90        \
	NWTC_Library_Types.f90 \
	NWTC_IO.f90            \
	NWTC_Num.f90           \
	ModMesh_Types.f90      \
	ModMesh.f90            \
	ModMesh_Mapping.f90    \
	NWTC_Library.f90

NETLIB_SOURCES=          \
	NWTC_FFTPACK.f90      \
	NWTC_LAPACK.f90       \
	fftpack4.1.f

pRNG_SOURCES=            \
	RANLUX.f90

TurbSim_SOURCES =        \
	TurbSim_Types.f90     \
	Root_Searching.f90    \
	RandNum.f90           \
	Profiles.f90          \
	BlankModVKM.f90       \
	CohStructures.f90     \
	TS_FileIO.f90         \
	TSsubs.f90            \
	TurbSim.f90           \
	VelocitySpectra.f90

vpath %.f90 $(NWTC_Lib_DIR) $(NETLIB_DIR) $(pRNG_DIR) $(TurbSim_DIR)
vpath %.f   $(NETLIB_DIR)
vpath %.mod $(INTER_DIR)
vpath %.obj $(INTER_DIR)


ALL_SOURCES = $(TurbSim_SOURCES) $(LIB_SOURCES) $(NETLIB_SOURCES) $(pRNG_SOURCES)
tmp_objs1   = $(ALL_SOURCES:.f90=.obj)
tmp_objs2   = $(tmp_objs1:.F90=.obj)       #note the upper case here (from IceFloe)
ALL_OBJS    = $(tmp_objs2:.f=.obj)


   # Rule to do everything.

all:     default
# use this for "all" for debugging: print variables:
#all:     ; $(info $$IceF_SOURCES is [${IceF_SOURCES}]) echo debugging
#
default:  $(INTER_DIR) $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT)


   # General rules for compliling the files.

%.obj: %.f90
	$(FC) -I $(INTER_DIR) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)

%.obj: %.F90
	$(FC) -I $(INTER_DIR) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)

%.obj: %.f
	$(F77) -I $(INTER_DIR) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)




   #  Dependency rules.
#NWTC Library dependency rules:
NWTC_Base.obj:              SingPrec.obj
$(SYS_FILE).obj:            NWTC_Base.obj
NWTC_Library_Types.obj:     $(SYS_FILE).obj
NWTC_IO.obj:                NWTC_Library_Types.obj
NWTC_Num.obj:               NWTC_IO.obj
ModMesh_Types.obj:          NWTC_Num.obj
ModMesh.obj:                ModMesh_Types.obj
ModMesh_Mapping.obj:        ModMesh.obj NWTC_LAPACK.obj
NWTC_Library.obj:           ModMesh.obj  ModMesh_Mapping.obj

NWTC_LAPACK.obj:            NWTC_Base.obj
NWTC_FFTPACK.obj:           NWTC_Library.obj fftpack4.1.obj

RANLUX.obj:                 NWTC_Library.obj

# TurbSim dependency rules:
TurbSim_Types.obj:          NWTC_Library.obj
BlankModVKM.obj:            NWTC_Library.obj
Root_Searching.obj:         NWTC_Library.obj
RandNum.obj:                TurbSim_Types.obj  RANLUX.obj
Profiles.obj:               TurbSim_Types.obj
VelocitySpectra.obj:        TurbSim_Types.obj
CohStructures.obj:          TurbSim_Types.obj Profiles.obj RandNum.obj
TSsubs.obj:                 TurbSim_Types.obj Profiles.obj RandNum.obj VelocitySpectra.obj BlankModVKM.obj NWTC_FFTPACK.obj NWTC_LAPACK.obj
TS_FileIO.obj:              TurbSim_Types.obj Profiles.obj RandNum.obj TSsubs.obj
TurbSim.obj:                TurbSim_Types.obj Profiles.obj TS_FileIO.obj CohStructures.obj


   # Make sure the destination directory for the intermediate files exist.

$(INTER_DIR):
	$(MD_CMD) $(INTER_DIR)


   # For linking main program.

$(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT): $(ALL_OBJS) | $(INTER_DIR)
	$(FC) $(LDFLAGS) -I $(INTER_DIR) -o $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT) \
	$(foreach src, $(ALL_OBJS), $(addprefix $(INTER_DIR)/,$(src))) $(LAPACK_LINK)

   # Cleanup afterwards.

clean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*.obj

