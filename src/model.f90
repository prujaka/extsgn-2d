module model
  use parameters
  implicit none
  contains

  real(8) function get_p(h, eta)
    use parameters
    implicit none
    real(kind=DP), intent(in) :: h, eta

    get_p = 0.5d0*gg*h*h - LAMBDA*eta*(eta/h-1.0d0)/3.0d0
    return
  end

  real(8) function get_a(h, eta)
    use parameters
    implicit none
    real(kind=DP), intent(in) :: h, eta

    get_a = DSQRT(gg*h + LAMBDA*eta*eta/(h*h)/3.0d0)
    return
  end

end module model
