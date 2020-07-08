module parameters
  implicit none
  integer, parameter       :: NX        = 20
  integer, parameter       :: NY        = 20
  integer, parameter       :: NEQS      = 5 ! Must be superior or equal to 3
  integer, parameter       :: DP        = kind(1.0d0)
  real(kind=DP), parameter :: XLEFT     = -300.0d0
  real(kind=DP), parameter :: XRIGHT    = 300.0d0
  real(kind=DP), parameter :: XMID      = 0.0d0
  real(kind=DP), parameter :: YLEFT     = -300.0d0
  real(kind=DP), parameter :: YRIGHT    = 300.0d0
  real(kind=DP), parameter :: YMID      = 0.0d0
  real(kind=DP), parameter :: DX        = (XRIGHT-XLEFT) / DFLOAT(NX)
  real(kind=DP), parameter :: DY        = (YRIGHT-YLEFT) / DFLOAT(NY)
  real(kind=dp), parameter :: GG        = 9.81d0
  real(kind=dp), parameter :: LAMBDA    = 1.0d0
  real(kind=DP)            :: dt        = 1.0d-8

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

  character(LEN=7), parameter :: OUTPUT_FILE = 'res.out'
  integer, parameter          :: OUTPUT_FILENAME_LENGTH = 7

end module parameters
