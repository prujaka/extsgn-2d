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

  subroutine get_prims_gn(prim,h,u,v,eta,w)
    use parameters
    implicit none
    real(kind=DP), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), dimension(0:NX+1,0:NY+1), intent(out) :: h,u,v,eta,w

    h(:,:) = prim(1,:,:)
    u(:,:) = prim(2,:,:)
    v(:,:) = prim(3,:,:)
    eta(:,:) = prim(4,:,:)
    w(:,:) = prim(5,:,:)

  end subroutine get_prims_gn

  subroutine ode_exact_solution(h,u,v,eta,w)
		use parameters
		implicit none
		real(kind=DP), dimension(0:NX+1,0:NY+1), intent(inout) :: h,u,v,eta,w
		real(kind=DP), dimension(0:NX+1,0:NY+1):: etanew

    h(:,:) = h(:,:)
    u(:,:) = u(:,:)
    etanew(:,:) = (eta(:,:)-h(:,:))*dcos( dsqrt(LAMBDA)*dt/h(:,:) )&
               + w(:,:)*h(:,:)*dsin( dsqrt(LAMBDA)*dt/h(:,:) )/dsqrt(LAMBDA)&
               + h(:,:)
    w(:,:) = -dsqrt(LAMBDA)*(etanew(:,:)-h(:,:))&
                *dsin( dsqrt(LAMBDA)*dt/h(:,:) )&
                /h(:,:) + w(:,:)*dcos( dsqrt(LAMBDA)*dt/h(:,:) )
    eta(:,:) = etanew(:,:)
		return
	end subroutine ode_exact_solution

  subroutine make_sources(h,eta,w,S)
			use parameters
			implicit none
			real(kind=DP), dimension(0:NX+1,0:NY+1), intent(in) :: h,eta,w
			real(kind=DP), intent(out) :: S(NEQS,0:NX+1,0:NY+1)
			integer :: i

			S(1,:,:) = 0.0d0
			S(2,:,:) = 0.0d0
			S(3,:,:) = 0.0d0
			S(4,:,:) = h(:,:)*w(:,:)
			S(5,:,:) = LAMBDA*(1.0d0 - eta(:,:)/h(:,:))

			return
		end subroutine make_sources

end module model
