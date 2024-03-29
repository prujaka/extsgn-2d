module m_parameters
  implicit none
  integer, parameter  :: dp        = kind(1.0d0)
  integer, parameter  :: NEQS      = 5 ! Number of equations

  ! Physical m_parameters of the model
  real(dp), parameter :: LAMBDA    = 75.0d0 ! Penalization parameter value
  real(dp), parameter :: GG        = 9.81d0 ! Gravity acceleration

  ! Time-related m_parameters
  real(dp), parameter :: TFIN      = 15.0d0 ! Final time
  real(dp), parameter :: CFL       = 0.5d0 ! CFL value
  integer, parameter  :: ITFIN     = 100000000 ! Max number of iterations
  real(dp)            :: dt        = 1.0d-8 ! initial dt

  ! Computational domain m_parameters: number of cells, and domain boudaries
  integer, parameter  :: NX        = 500
  integer, parameter  :: NY        = 500
  real(dp), parameter :: XL        = -300.0d0
  real(dp), parameter :: XR        = 300.0d0
  real(dp), parameter :: YL        = -300.0d0
  real(dp), parameter :: YR        = 300.0d0

  real(dp), parameter :: DX        = (XR-XL) / DFLOAT(NX)
  real(dp), parameter :: DY        = (YR-YL) / DFLOAT(NY)
  real(dp), parameter :: DV        = DX*DY
  real(dp), parameter :: DL        = dmin1(DX,DY)

  ! Cylinder radius for the 2D cylindric Riemann problem (RP)
  real(dp), parameter :: XMID      = 0.0d0
  real(dp), parameter :: YMID      = 0.0d0
  real(dp), parameter :: RADIUS    = 40.0d0

  ! Riemann problem initial values of the unknowns (left or inner)
  real(dp), parameter :: HL_INIT   = 1.8d0
  real(dp), parameter :: UL_INIT   = 0.0d0
  real(dp), parameter :: VL_INIT   = 0.0d0
  real(dp), parameter :: ETAL_INIT = 1.8d0
  real(dp), parameter :: WL_INIT   = 0.0d0

  ! Riemann problem initial values of the unknowns (right or outer)
  real(dp), parameter :: HR_INIT   = 1.0d0
  real(dp), parameter :: UR_INIT   = 0.0d0
  real(dp), parameter :: VR_INIT   = 0.0d0
  real(dp), parameter :: ETAR_INIT = 1.0d0
  real(dp), parameter :: WR_INIT   = 0.0d0


  ! IMEX Delta parameter (see eq. (25) in the paper)
  real(dp), parameter :: DELTA     = 0.2928932188134524d0
  real(dp), parameter :: OMEGA     = 0.0d0

  ! Boundary conditions velocity coefficients
  real(dp), parameter :: BC_U_LEFT = -1.0d0
  real(dp), parameter :: BC_U_RIGHT= -1.0d0
  real(dp), parameter :: BC_V_LEFT = -1.0d0
  real(dp), parameter :: BC_V_RIGHT= -1.0d0

  ! Output files specifications
  integer, parameter                :: NAMELEN     = 22
  character(LEN=NAMELEN), parameter :: FILE_DAT    = 'out/res.dat'
  character(LEN=NAMELEN), parameter :: FILE_VTK    = 'out/res.vtk'
  character(LEN=13)                 :: DATA_FORMAT = '(9(e11.4,1x))'
  integer, parameter                :: VTK_STEP    = 1
  integer, parameter                :: N_FILES     = 100
  integer, parameter                :: GENERATE_VTK  = 0

  ! Numerical method specifier constants
  integer, parameter :: METHOD_FIRST_ORDER_SPLITTING = 0
  integer, parameter :: METHOD_IMEX_ARS_222          = 1

  ! Initial condition specifier constants
  integer, parameter :: IC_RP_X            = 0
  integer, parameter :: IC_RP_Y            = 1
  integer, parameter :: IC_RP_CYL          = 2
  integer, parameter :: IC_RP_SQR          = 3
  integer, parameter :: IC_GENRP_SINX      = 4
  integer, parameter :: IC_GENRP_SINX_SINY = 5
  integer, parameter :: IC_MATRIX          = 6

  ! Percentage output enumerator
  enum, bind(c)
    enumerator :: PERC_OUTPUT_OFF, PERC_OUTPUT_ON
  end enum

  ! Initial condition and Numerical method specifier flags
  integer, parameter :: SELECTOR_IC = IC_MATRIX
  integer, parameter :: SELECTOR_METHOD = METHOD_IMEX_ARS_222
  integer, parameter :: SELECTOR_PERC_OUTPUT = PERC_OUTPUT_ON

end module m_parameters
