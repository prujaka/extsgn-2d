module parameters
  implicit none
  integer, parameter       :: DP        = kind(1.0d0)
  integer, parameter       :: NEQS      = 5 ! Must be superior or equal to 5

  integer, parameter       :: NX        = 200
  integer, parameter       :: NY        = 200
  real(kind=dp), parameter :: LAMBDA    = 10.0d0
  real(kind=dp), parameter :: TFIN      = 40.0d0
  real(kind=dp), parameter :: CFL       = 0.5d0
  integer, parameter       :: ITFIN     = 100000000

  real(kind=DP), parameter :: XL        = -300.0d0
  real(kind=DP), parameter :: XR        = 300.0d0
  real(kind=DP), parameter :: XMID      = 0.0d0
  real(kind=DP), parameter :: YL        = -300.0d0
  real(kind=DP), parameter :: YR        = 300.0d0
  real(kind=DP), parameter :: YMID      = 0.0d0
  real(kind=dp), parameter :: RADIUS    = 40.0d0
  real(kind=dp), parameter :: GG        = 9.810d0

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

  real(kind=DP), parameter :: DX        = (XR-XL) / DFLOAT(NX)
  real(kind=DP), parameter :: DY        = (YR-YL) / DFLOAT(NY)
  real(kind=DP), parameter :: DV        = DX*DY
  real(kind=DP), parameter :: DL        = dmin1(DX,DY)
  real(kind=DP)            :: dt        = 1.0d-8
  character(LEN=13)        :: DATA_FORMAT = '(9(e11.4,1x))'
  real(kind=DP), parameter :: PERC_FREQ = 1.0d0


  real(kind=DP), parameter :: DELTA     = 0.2928932188134524d0
  real(kind=DP), parameter :: OMEGA     = 0.0d0

  real(kind=DP), parameter :: BC_U_LEFT = -1.0d0
  real(kind=DP), parameter :: BC_U_RIGHT= -1.0d0
  real(kind=DP), parameter :: BC_V_LEFT = -1.0d0
  real(kind=DP), parameter :: BC_V_RIGHT= -1.0d0

  integer, parameter       :: NAMELEN = 7
  integer, parameter       :: VTK_STEP = 1
  character(LEN=NAMELEN), parameter :: FILE_DAT = 'res.dat'
  character(LEN=NAMELEN), parameter :: FILE_VTK = 'res.vtk'

  integer, parameter       :: METHOD_GODUNOV = 0
  integer, parameter       :: METHOD_IMEX = 1
  integer, parameter       :: IC_RP_X     = 0
  integer, parameter       :: IC_RP_Y     = 1
  integer, parameter       :: IC_RP_CYL   = 2
  integer, parameter       :: IC_RP_SQR   = 3

  integer, parameter       :: SELECTOR_IC = IC_RP_SQR
  integer, parameter       :: SELECTOR_METHOD = METHOD_IMEX

end module parameters
