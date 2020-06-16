module methods
  use parameters
  implicit none
  contains

  subroutine set_mesh(x,y)
    use parameters
    implicit none
    real(kind=DP), intent(out) :: x(0:NX+1)
    real(kind=DP), intent(out) :: y(0:NY+1)
    integer :: i

    do i=0,NX+1
      x(i) = XLEFT + 0.5d0*DX + DFLOAT(i-1)*DX
    enddo
    do i=0,NY+1
      y(i) = YLEFT + 0.5d0*DY + DFLOAT(i-1)*DY
    enddo

    return
  end

  subroutine set_ic_rpx(x,prim)
    use parameters
    implicit none
    real(kind=DP), intent(in)  :: x(0:NX+1)
    real(kind=DP), intent(out) :: prim(NEQS, 0:NX+1, 0:NY+1)
    integer :: i, j

    do j=1, NY
      do i=1, NX
        if (x(i).le.XMID) then
          prim(1,i,j) = HL_INIT
          prim(2,i,j) = UL_INIT
          prim(3,i,j) = VL_INIT
          prim(4,i,j) = ETAL_INIT
          prim(5,i,j) = WL_INIT
        else
          prim(1,i,j) = HR_INIT
          prim(2,i,j) = UR_INIT
          prim(3,i,j) = VR_INIT
          prim(4,i,j) = ETAR_INIT
          prim(5,i,j) = WR_INIT
        endif
      enddo
    enddo

    return
  end

end module methods
