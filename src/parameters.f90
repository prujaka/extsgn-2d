module parameters
  implicit none
  integer, parameter       :: NX        = 2000
  integer, parameter       :: NY        = 2000
  integer, parameter       :: NEQS      = 5 ! Must be superior or equal to 5
  integer, parameter       :: DP        = kind(1.0d0)
  integer, parameter       :: ITFINAL   = 100000000
  real(kind=dp), parameter :: TIMEFINAL = 40.0d0
  real(kind=dp), parameter :: CFL       = 0.5d0
  real(kind=dp), parameter :: LAMBDA    = 50.0d0
  real(kind=DP), parameter :: DELTA     = 0.2928932188134524d0
  real(kind=DP), parameter :: XLEFT     = -300.0d0
  real(kind=DP), parameter :: XRIGHT    = 300.0d0
  real(kind=DP), parameter :: XMID      = 0.0d0
  real(kind=DP), parameter :: YLEFT     = -300.0d0
  real(kind=DP), parameter :: YRIGHT    = 300.0d0
  real(kind=DP), parameter :: YMID      = 0.0d0
  real(kind=dp), parameter :: GG        = 9.810d0
  real(kind=dp), parameter :: RADIUS    = 20.0d0
  real(kind=DP), parameter :: DX        = (XRIGHT-XLEFT) / DFLOAT(NX)
  real(kind=DP), parameter :: DY        = (YRIGHT-YLEFT) / DFLOAT(NY)
  real(kind=DP), parameter :: DV        = DX*DY
  real(kind=DP), parameter :: DL        = dmin1(DX,DY)
  real(kind=DP)            :: dt        = 1.0d-8
  real(kind=DP), parameter :: PERC_FREQ = 1.0d0
  character(LEN=13)        :: DATA_FORMAT = '(9(e11.4,1x))'

  real(kind=DP), parameter :: HL_INIT   = 1.8d0
  real(kind=DP), parameter :: UL_INIT   = 0.0d0
  real(kind=DP), parameter :: VL_INIT   = 0.0d0
  real(kind=DP), parameter :: ETAL_INIT = 1.8d0
  real(kind=DP), parameter :: WL_INIT   = 0.0d0

  real(kind=DP), parameter :: HR_INIT   = 1.0d0
  real(kind=DP), parameter :: UR_INIT   = 0.0d0
  real(kind=DP), parameter :: VR_INIT   = 0.0d0
  real(kind=DP), parameter :: ETAR_INIT = 1.0d0
  real(kind=DP), parameter :: WR_INIT   = 0.0d0

  real(kind=DP), parameter :: OMEGA     = 0.0d0

  integer, parameter     :: IC_RP_X     = 0
  integer, parameter     :: IC_RP_Y     = 1
  integer, parameter     :: IC_RP_CYL   = 2
  integer, parameter     :: METHOD_GODUNOV = 0
  integer, parameter     :: METHOD_IMEX = 1

  integer, parameter     :: SELECTOR_IC = IC_RP_CYL
  integer, parameter     :: SELECTOR_METHOD = METHOD_IMEX

  real(kind=DP), parameter :: BC_U_LEFT = -1.0d0
  real(kind=DP), parameter :: BC_U_RIGHT= -1.0d0
  real(kind=DP), parameter :: BC_V_LEFT = -1.0d0
  real(kind=DP), parameter :: BC_V_RIGHT= -1.0d0

  character(LEN=7), parameter :: OUTPUT_FILE = 'res.out'
  integer, parameter          :: OUTPUT_FILENAME_LENGTH = 7

end module parameters
