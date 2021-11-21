module parameters
  implicit none
  integer, parameter  :: dp        = kind(1.0d0)
  integer, parameter  :: NEQS      = 5 ! Must be superior or equal to 5

  integer, parameter  :: NX        = 800
  integer, parameter  :: NY        = 800
  real(dp), parameter :: LAMBDA    = 75.0d0
  real(dp), parameter :: TFIN      = 40.0d0
  real(dp), parameter :: CFL       = 0.45d0
  real(dp), parameter :: GG        = 9.810d0
  integer, parameter  :: ITFIN     = 100000000

  real(dp), parameter :: XL        = -300.0d0
  real(dp), parameter :: XR        = 300.0d0
  real(dp), parameter :: XMID      = 0.0d0
  real(dp), parameter :: YL        = -300.0d0
  real(dp), parameter :: YR        = 300.0d0
  real(dp), parameter :: YMID      = 0.0d0
  real(dp), parameter :: RADIUS    = 40.0d0

  real(dp), parameter :: HL_INIT   = 1.8d0
  real(dp), parameter :: UL_INIT   = 0.0d0
  real(dp), parameter :: VL_INIT   = 0.0d0
  real(dp), parameter :: ETAL_INIT = 1.8d0
  real(dp), parameter :: WL_INIT   = 0.0d0

  real(dp), parameter :: HR_INIT   = 1.0d0
  real(dp), parameter :: UR_INIT   = 0.0d0
  real(dp), parameter :: VR_INIT   = 0.0d0
  real(dp), parameter :: ETAR_INIT = 1.0d0
  real(dp), parameter :: WR_INIT   = 0.0d0

  real(dp), parameter :: DX        = (XR-XL) / DFLOAT(NX)
  real(dp), parameter :: DY        = (YR-YL) / DFLOAT(NY)
  real(dp), parameter :: DV        = DX*DY
  real(dp), parameter :: DL        = dmin1(DX,DY)
  real(dp)            :: dt        = 1.0d-8
  character(LEN=13)        :: DATA_FORMAT = '(9(e11.4,1x))'
  real(dp), parameter :: PERC_FREQ = 1.0d0


  real(dp), parameter :: DELTA     = 0.2928932188134524d0
  real(dp), parameter :: OMEGA     = 0.0d0

  real(dp), parameter :: BC_U_LEFT = -1.0d0
  real(dp), parameter :: BC_U_RIGHT= -1.0d0
  real(dp), parameter :: BC_V_LEFT = -1.0d0
  real(dp), parameter :: BC_V_RIGHT= -1.0d0

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
  integer, parameter       :: IC_GENRP_SINX = 4
  integer, parameter       :: IC_GENRP_SINX_SINY = 5

  integer, parameter       :: SELECTOR_IC = IC_RP_SQR
  integer, parameter       :: SELECTOR_METHOD = METHOD_IMEX

end module parameters
