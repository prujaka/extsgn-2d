module model
  use parameters
  implicit none
  contains

  real(8) function pressure(h, eta)
    implicit none
    real(dp), intent(in) :: h, eta

    pressure = 0.5d0*gg*h*h - LAMBDA*eta*(eta/h-1.0d0)/3.0d0
  end

  real(8) function sound_speed(h, eta)
    implicit none
    real(dp), intent(in) :: h, eta

    sound_speed = DSQRT(gg*h + LAMBDA*eta*eta/(h*h)/3.0d0)
  end

  real(8) function compute_dt(h,u,v,eta)
    implicit none
    real(dp), dimension(0:NX+1, 0:NY+1), intent(in) :: h, u, v, eta
    real(dp) :: cmax, c
    integer :: i,j

    cmax = 0.0d0
    do i=1,NX
      do j=1,NY
      c = sound_speed(h(i,j),eta(i,j))
      cmax = max( cmax, (dabs(u(i,j)) + c)/dx + (dabs(v(i,j)) + c)/dy )
      enddo
    enddo
    compute_dt = CFL/cmax

  end

end module model
