module methods
  use parameters
  use model
  use aux
  implicit none
  contains

  subroutine initialize_problem(x,y,prim,cons,it,time)
    implicit none
    real(dp), intent(out) :: x(0:NX+1), y(0:NY+1), time
    real(dp), dimension(NEQS, 0:NX+1, 0:NY+1),intent(out) :: prim, cons
    integer, intent(out) :: it
    character(len=20) :: filename = 'img/file.txt'

    call set_mesh(x,y)
    ! call set_ic(x,y,prim)
    call input_matrix_flat(prim, filename)
    call prim_to_cons(prim,cons)
    it = 0
    time = 0.0d0

  end subroutine initialize_problem

  subroutine set_mesh(x,y)
    implicit none
    real(dp), intent(out) :: x(0:NX+1)
    real(dp), intent(out) :: y(0:NY+1)
    integer :: i

    do i=0,NX+1
      x(i) = XL + 0.5d0*DX + DFLOAT(i-1)*DX
    enddo
    do i=0,NY+1
      y(i) = YL + 0.5d0*DY + DFLOAT(i-1)*DY
    enddo

  end subroutine set_mesh

  subroutine set_ic(x,y,prim)
    implicit none
    real(dp), intent(in)  :: x(0:NX+1), y(0:NY+1)
    real(dp), intent(out) :: prim(NEQS, 0:NX+1, 0:NY+1)

    select case(SELECTOR_IC)
      case(IC_RP_X)
        call set_ic_rp_x(x,prim)
      case(IC_RP_Y)
        call set_ic_rp_y(y,prim)
      case(IC_RP_CYL)
        call set_ic_rp_cyl(x,y,prim)
      case(IC_RP_SQR)
        call set_ic_rp_sqr(x,y,prim)
    end select
  end subroutine set_ic

  subroutine set_ic_rp_x(x,prim)
    implicit none
    real(dp), intent(in)  :: x(0:NX+1)
    real(dp), intent(out) :: prim(NEQS, 0:NX+1, 0:NY+1)
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

  end subroutine set_ic_rp_x
  subroutine set_ic_rp_y(y,prim)
    implicit none
    real(dp), intent(in)  :: y(0:NY+1)
    real(dp), intent(out) :: prim(NEQS, 0:NX+1, 0:NY+1)
    integer :: i, j

    do j=1, NY
      do i=1, NX
        if (y(j).le.YMID) then
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

  end subroutine set_ic_rp_y
  subroutine set_ic_rp_cyl(x,y,prim)
    implicit none
    real(dp), intent(in)  :: x(0:NX+1), y(0:NY+1)
    real(dp), intent(out) :: prim(NEQS, 0:NX+1, 0:NY+1)
    integer :: i, j

    do j=1, NY
      do i=1, NX
        if (x(i)*x(i)+y(j)*y(j).le.RADIUS*RADIUS) then
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

  end subroutine set_ic_rp_cyl
  subroutine set_ic_rp_sqr(x,y,prim)
    implicit none
    real(dp), intent(in)  :: x(0:NX+1), y(0:NY+1)
    real(dp), intent(out) :: prim(NEQS, 0:NX+1, 0:NY+1)
    integer :: i, j

    do j=1, NY
      do i=1, NX
        if ((dabs(x(i)-XMID).le.RADIUS) .and. (dabs(y(j)-YMID).le.RADIUS)) then
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

  end subroutine set_ic_rp_sqr

  subroutine prim_to_cons(prim,cons)
    implicit none
    real(dp), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(dp), intent(out) :: cons(NEQS,0:NX+1,0:NY+1)
    integer :: k

    cons(1,:,:) = prim(1,:,:)
    do k=2, NEQS
      cons(k,:,:) = cons(1,:,:)*prim(k,:,:)
    enddo

  end subroutine prim_to_cons

  subroutine cons_to_prim(cons,prim)
    implicit none
    real(dp), intent(in) :: cons(NEQS,0:NX+1,0:NY+1)
    real(dp), intent(out) :: prim(NEQS,0:NX+1,0:NY+1)
    integer :: k

    prim(1,:,:) = cons(1,:,:)
    do k=2, NEQS
      prim(k,:,:) = cons(k,:,:)/prim(1,:,:)
    enddo

  end subroutine cons_to_prim

  subroutine get_solution(prim,cons,it,time)
    implicit none
    real(dp), intent(inout) :: cons(NEQS,0:NX+1,0:NY+1)
    real(dp), intent(inout) :: prim(NEQS,0:NX+1,0:NY+1)
    integer, intent(inout) :: it
    real(dp), intent(inout) :: time
    real(dp) :: milestone=0.0d0, t1

    print*, ''
    print*, 'Calculation started.'

    do while(time<TFIN)
      if (it.ge.ITFIN) exit
      dt = compute_dt(prim(1,:,:),prim(2,:,:),prim(3,:,:),prim(4,:,:))
      if (dt>TFIN-time) dt = TFIN-time

      select case(SELECTOR_METHOD)
        case(METHOD_GODUNOV)
          call timestep_godunov(prim,cons)
        case(METHOD_IMEX)
          call timestep_imex(prim)
      end select

      call print_percentage(time,t1,milestone)

      it=it+1
      time=time+dt
    enddo

  end subroutine get_solution

  subroutine timestep_godunov(prim,cons)
    implicit none
    real(dp), intent(inout) :: cons(NEQS,0:NX+1,0:NY+1)
    real(dp), intent(inout) :: prim(NEQS,0:NX+1,0:NY+1)

    real(dp), allocatable :: F(:,:,:), G(:,:,:), S(:,:,:)
    real(dp), allocatable :: h(:,:), u(:,:), v(:,:), eta(:,:), w(:,:)

    allocate(F(NEQS,0:NX+1,0:NY+1),h(0:NX+1,0:NY+1))
    allocate(G,S,mold=F)
    allocate(u,v,eta,w,mold=h)

    call set_bc(prim)
    call riemann_fluxes_x(prim,F)
    call riemann_fluxes_y(prim,G)
    call godunov(cons,F,G)
    call split_prims_gn(prim,h,u,v,eta,w)
    call ode_exact_solution(h,u,v,eta,w)
    call make_sources(h,eta,w,S)
    call ode_euler_step(S,cons)
    call cons_to_prim(cons,prim)

    deallocate(h, u, v, eta, w, F, G, S)

  end subroutine timestep_godunov

  subroutine timestep_imex(prim)
    implicit none
    real(dp), intent(inout) :: prim(NEQS,0:NX+1,0:NY+1)
    real(dp), allocatable :: F(:,:,:),G(:,:,:)
    real(dp), allocatable :: Fstar(:,:,:),Gstar(:,:,:),S(:,:,:)
    real(dp), allocatable :: h(:,:),u(:,:),v(:,:),eta(:,:),w(:,:)
    real(dp), allocatable :: hstar(:,:),ustar(:,:),vstar(:,:)
    real(dp), allocatable :: etastar(:,:),wstar(:,:)
    real(dp), allocatable :: primstar(:,:,:)

    allocate(F(NEQS,0:NX+1,0:NY+1),h(0:NX+1,0:NY+1))
    allocate(G,S,Fstar,Gstar,primstar,mold=F)
    allocate(u,v,eta,w,hstar,ustar,vstar,etastar,wstar,mold=h)

    call set_bc(prim)
    call muscl_step_x(prim,F)
    call muscl_step_y(prim,G)
    call split_prims_gn(prim,h,u,v,eta,w)
    call imex_step_1(h,u,v,eta,w,hstar,ustar,vstar,etastar,wstar,F,G)
    call make_sources(hstar,etastar,wstar,S)
    call merge_prims_gn(hstar,ustar,vstar,etastar,wstar,primstar)
    call set_bc(primstar)
    call muscl_step_x(primstar,Fstar)
    call muscl_step_y(primstar,Gstar)
    call imex_step_2(h,u,v,eta,w,F,G,Fstar,Gstar,S)
    call merge_prims_gn(h,u,v,eta,w,prim)

    deallocate(F,G,Fstar)
    deallocate(Gstar,S)
    deallocate(h,u,v,eta)
    deallocate(w)

    deallocate(hstar,ustar,vstar)
    deallocate(etastar,wstar)
    deallocate(primstar)

  end subroutine timestep_imex

  subroutine set_bc(prim)
    implicit none
    real(dp), intent(inout) :: prim(NEQS,0:NX+1,0:NY+1)
    integer :: i,j

    do j=1,NY
      prim(:,0,j) = prim(:,1,j)
      prim(2,0,j) = BC_U_LEFT*prim(2,1,j)

      prim(:,NX+1,j) = prim(:,NX,j)
      prim(2,NX+1,j) = BC_U_RIGHT*prim(2,NX,j)
    enddo
    do i=1,NX
      prim(:,i,0) = prim(:,i,1)
      prim(3,i,0) = BC_V_LEFT*prim(3,i,1)

      prim(:,i,NY+1) = prim(:,i,NY)
      prim(3,i,NY+1) = BC_V_RIGHT*prim(3,i,NY)
    enddo


  end subroutine set_bc

  subroutine riemann_fluxes_x(prim,flux)
    implicit none
    real(dp), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(dp), intent(out) :: flux(NEQS,0:NX+1,0:NY+1)
    real(dp), dimension(NEQS) :: priml, primr, F
    integer :: i,j

    do j=1,NY
      do i=0,NX
        priml = prim(:,i,j)
        primr = prim(:,i+1,j)
        call hllc(priml,primr,F)
        flux(:,i,j) = F
      enddo
    enddo

  end subroutine riemann_fluxes_x
  subroutine riemann_fluxes_y(prim,flux)
    implicit none
    real(dp), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(dp), intent(out) :: flux(NEQS,0:NX+1,0:NY+1)
    real(dp), dimension(NEQS) :: priml, primr, F
    integer :: i,j

    do i=1,NX
      do j=0,NY
        priml = prim(:,i,j)
        priml(2) = prim(3,i,j)
        priml(3) = prim(2,i,j)

        primr = prim(:,i,j+1)
        primr(2) = prim(3,i,j+1)
        primr(3) = prim(2,i,j+1)

        call hllc(priml,primr,F)

        flux(:,i,j) = F
        flux(2,i,j) = F(3)
        flux(3,i,j) = F(2)
      enddo
    enddo

  end subroutine riemann_fluxes_y

  subroutine godunov(cons,F,G)
    implicit none

    real(dp), dimension(NEQS,0:NX+1,0:NY+1), intent(in) :: F,G
    real(dp), dimension(NEQS,0:NX+1,0:NY+1), intent(inout) :: cons
    integer :: i,j

    do j=1,NY
      do i=1,NX
        cons(:,i,j)=cons(:,i,j)-dt/dV*( dy*(F(:,i,j)-F(:,i-1,j))&
                                       +DX*(G(:,i,j)-G(:,i,j-1)) )
      enddo
    enddo
  end subroutine godunov

  subroutine split_prims_gn(prim,h,u,v,eta,w)
    implicit none
    real(dp), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(dp), dimension(0:NX+1,0:NY+1), intent(out) :: h,u,v,eta,w

    h(:,:) = prim(1,:,:)
    u(:,:) = prim(2,:,:)
    v(:,:) = prim(3,:,:)
    eta(:,:) = prim(4,:,:)
    w(:,:) = prim(5,:,:)

  end subroutine split_prims_gn

  subroutine merge_prims_gn(h,u,v,eta,w,prim)
    implicit none
    real(dp), dimension(0:NX+1,0:NY+1), intent(in) :: h,u,v,eta,w
    real(dp), intent(out) :: prim(NEQS,0:NX+1,0:NY+1)

    prim(1,:,:) = h(:,:)
    prim(2,:,:) = u(:,:)
    prim(3,:,:) = v(:,:)
    prim(4,:,:) = eta(:,:)
    prim(5,:,:) = w(:,:)

  end subroutine merge_prims_gn

  subroutine ode_exact_solution(h,u,v,eta,w)
    implicit none
    real(dp), dimension(0:NX+1,0:NY+1), intent(inout) :: h,u,v,eta,w
    real(dp) :: etanew(0:NX+1,0:NY+1)

    h(:,:) = h(:,:)
    u(:,:) = u(:,:)
    v(:,:) = v(:,:)
    etanew(:,:) = (eta(:,:)-h(:,:))*dcos( dsqrt(LAMBDA)*dt/h(:,:) )&
               + w(:,:)*h(:,:)*dsin( dsqrt(LAMBDA)*dt/h(:,:) )/dsqrt(LAMBDA)&
               + h(:,:)
    w(:,:) = -dsqrt(LAMBDA)*(etanew(:,:)-h(:,:))&
                *dsin( dsqrt(LAMBDA)*dt/h(:,:) )&
                /h(:,:) + w(:,:)*dcos( dsqrt(LAMBDA)*dt/h(:,:) )
    eta(:,:) = etanew(:,:)
  end subroutine ode_exact_solution

  subroutine make_sources(h,eta,w,S)
    implicit none
    real(dp), dimension(0:NX+1,0:NY+1), intent(in) :: h,eta,w
    real(dp), intent(out) :: S(NEQS,0:NX+1,0:NY+1)

    S(1,:,:) = 0.0d0
    S(2,:,:) = 0.0d0
    S(3,:,:) = 0.0d0
    S(4,:,:) = h(:,:)*w(:,:)
    S(5,:,:) = LAMBDA*(1.0d0 - eta(:,:)/h(:,:))

  end subroutine make_sources

  subroutine ode_euler_step(S,cons)
    implicit none
    real(dp), intent(in) :: S(NEQS,0:NX+1,0:NY+1)
    real(dp), intent(inout) :: cons(NEQS,0:NX+1,0:NY+1)

    cons(:,:,:) = cons(:,:,:) + dt * S(:,:,:)
  end subroutine ode_euler_step

  subroutine muscl_step_x(prim,F)
    implicit none
    real(dp), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(dp), intent(out) :: F(NEQS,0:NX+1,0:NY+1)
    real(dp), dimension(NEQS,0:NX+1,0:NY+1) :: primr,priml

    call data_reconstruction_x(prim,primr,priml)
    call set_bc_muscl_x(primr,priml)
    call riemann_fluxes_muscl_x(primr,priml,F)

  end subroutine muscl_step_x

  subroutine data_reconstruction_x(prim,primr,priml)
    implicit none
    real(dp), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(dp), dimension(NEQS,0:NX+1,0:NY+1), intent(out) :: primr,priml
    real(dp) :: slope(NEQS,0:NX+1,0:NY+1)

    call get_slopes_x(prim,slope)
    primr(:,:,:) = prim(:,:,:) - 0.5d0*slope(:,:,:)
    priml(:,:,:) = prim(:,:,:) + 0.5d0*slope(:,:,:)
  end subroutine data_reconstruction_x

  subroutine get_slopes_x(prim,slope)
    implicit none
    real(dp), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(dp), intent(out) :: slope(NEQS,0:NX+1,0:NY+1)
    integer :: i

    do i=1,NX
      slope(:,i,:) = 0.5d0*(1.0d0+OMEGA)*(prim(:,i,:) - prim(:,i-1,:))&
                   + 0.5d0*(1.0d0-OMEGA)*(prim(:,i+1,:) - prim(:,i,:))
    enddo
  end subroutine get_slopes_x

  subroutine set_bc_muscl_x(primr,priml)
    implicit none
    real(dp), dimension(NEQS,0:NX+1,0:NY+1), intent(inout) :: primr,priml
    integer :: j

    do j=1,NY
      priml(:,0,j) = primr(:,1,j)
      priml(2,0,j) = BC_U_LEFT*primr(2,1,j)

      primr(:,NX+1,j) = priml(:,NX,j)
      primr(2,NX+1,j) = BC_U_RIGHT*priml(2,NX,j)
    enddo

  end subroutine set_bc_muscl_x

  subroutine riemann_fluxes_muscl_x(primr,priml,flux)
    implicit none
    real(dp), dimension(NEQS,0:NX+1,0:NY+1), intent(in) :: primr,priml
    real(dp), intent(out) :: flux(NEQS,0:NX+1,0:NY+1)
    real(dp), dimension(NEQS) :: statel, stater, F
    integer :: i,j

    do j=1,NY
      do i=0,NX
        statel = priml(:,i,j)
        stater = primr(:,i+1,j)
        call hllc(statel,stater,F)
        flux(:,i,j) = F
      enddo
    enddo
  end subroutine riemann_fluxes_muscl_x

  subroutine muscl_step_y(prim,F)
    implicit none
    real(dp), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(dp), intent(out) :: F(NEQS,0:NX+1,0:NY+1)
    real(dp), dimension(NEQS,0:NX+1,0:NY+1) :: primr,priml

    call data_reconstruction_y(prim,primr,priml)
    call set_bc_muscl_y(primr,priml)
    call riemann_fluxes_muscl_y(primr,priml,F)

  end subroutine muscl_step_y

  subroutine data_reconstruction_y(prim,primr,priml)
    implicit none
    real(dp), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(dp), dimension(NEQS,0:NX+1,0:NY+1), intent(out) :: primr,priml
    real(dp) :: slope(NEQS,0:NX+1,0:NY+1)

    call get_slopes_y(prim,slope)
    primr(:,:,:) = prim(:,:,:) - 0.5d0*slope(:,:,:)
    priml(:,:,:) = prim(:,:,:) + 0.5d0*slope(:,:,:)
  end subroutine data_reconstruction_y

  subroutine get_slopes_y(prim,slope)
    implicit none
    real(dp), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(dp), intent(out) :: slope(NEQS,0:NX+1,0:NY+1)
    integer :: j

    do j=1,NY
      slope(:,:,j) = 0.5d0*(1.0d0+OMEGA)*(prim(:,:,j) - prim(:,:,j-1))&
                   + 0.5d0*(1.0d0-OMEGA)*(prim(:,:,j+1) - prim(:,:,j))
    enddo
  end subroutine get_slopes_y

  subroutine set_bc_muscl_y(primr,priml)
    implicit none
    real(dp), dimension(NEQS,0:NX+1,0:NY+1), intent(inout) :: primr,priml
    integer :: i

    do i=1,NX
      priml(:,i,0) = primr(:,i,1)
      priml(3,i,0) = BC_V_LEFT*primr(3,i,1)

      primr(:,i,NY+1) = priml(:,i,NY)
      primr(3,i,NY+1) = BC_V_RIGHT*priml(3,i,NY)
    enddo

  end subroutine set_bc_muscl_y

  subroutine riemann_fluxes_muscl_y(primr,priml,flux)
    implicit none
    real(dp), dimension(NEQS,0:NX+1,0:NY+1), intent(in) :: primr,priml
    real(dp), intent(out) :: flux(NEQS,0:NX+1,0:NY+1)
    real(dp), dimension(NEQS) :: statel, stater, F
    integer :: i,j

    do i=1,NX
      do j=0,NY
        statel = priml(:,i,j)
        statel(2) = priml(3,i,j)
        statel(3) = priml(2,i,j)

        stater = primr(:,i,j+1)
        stater(2) = primr(3,i,j+1)
        stater(3) = primr(2,i,j+1)

        call hllc(statel,stater,F)

        flux(:,i,j) = F
        flux(2,i,j) = F(3)
        flux(3,i,j) = F(2)
      enddo
    enddo

  end subroutine riemann_fluxes_muscl_y

  subroutine imex_step_1(h,u,v,eta,w,hstar,ustar,vstar,etastar,wstar,F,G)
    implicit none
    real(dp), dimension(0:NX+1,0:NY+1), intent(in) :: h,u,v,eta,w
    real(dp), dimension(NEQS,0:NX+1,0:NY+1), intent(in) :: F, G
    real(dp), dimension(0:NX+1,0:NY+1), intent(out) :: hstar,ustar,vstar,&
                                                            etastar,wstar
    real(dp) :: alpha, C4, C5
    integer :: i,j

    do i=1,NX
      do j=1,NY
        hstar(i,j) = h(i,j) - DELTA*(dy*(F(1,i,j)-F(1,i-1,j))&
                                     +DX*(G(1,i,j)-G(1,i,j-1)))*dt/dV
        ustar(i,j) = ( h(i,j)*u(i,j) - DELTA*(dy*(F(2,i,j)-F(2,i-1,j))&
                         +DX*(G(2,i,j)-G(2,i,j-1)))*dt/dV )/hstar(i,j)
        vstar(i,j) = ( h(i,j)*v(i,j) - DELTA*(dy*(F(3,i,j)-F(3,i-1,j))&
                         +DX*(G(3,i,j)-G(3,i,j-1)))*dt/dV )/hstar(i,j)
        C4 = ( h(i,j)*eta(i,j) - DELTA*(dy*(F(4,i,j)-F(4,i-1,j))&
                         +DX*(G(4,i,j)-G(4,i,j-1)))*dt/dV )/hstar(i,j)
        C5 = ( h(i,j)*w(i,j) - DELTA*(dy*(F(5,i,j)-F(5,i-1,j))&
                         +DX*(G(5,i,j)-G(5,i,j-1)))*dt/dV )/hstar(i,j)
        alpha = 1.0d0/(1.0d0 +LAMBDA*DELTA*DELTA*dt*dt/(hstar(i,j)*hstar(i,j)))

        etastar(i,j) = alpha*( C4 + DELTA*dt*C5&
                      + DELTA*DELTA*dt*dt*LAMBDA/hstar(i,j) )

        wstar(i,j) = alpha*( C5 - C4*DELTA*dt*LAMBDA/(hstar(i,j)*hstar(i,j))&
                           + DELTA*dt*LAMBDA/hstar(i,j) )
      enddo
    enddo
  end subroutine imex_step_1

  subroutine imex_step_2(h,u,v,eta,w,F,G,Fstar,Gstar,S)
    implicit none
    real(dp), dimension(0:NX+1,0:NY+1), intent(inout) :: h,u,v,eta,w
    real(dp), dimension(NEQS,0:NX+1,0:NY+1),intent(in):: F,G,Fstar,Gstar,S
    real(dp) :: alpha, C4, C5, hnew
    integer :: i,j

    do i=1,NX
      do j=1,NY
        hnew = h(i,j) - (DELTA-1.0d0)*(dy*(F(1,i,j)-F(1,i-1,j))&
                        + DX*(G(1,i,j)-G(1,i,j-1)))*dt/dV &
                      - (2.0d0 - DELTA)*(dy*(Fstar(1,i,j)-Fstar(1,i-1,j))&
                        + DX*(Gstar(1,i,j)-Gstar(1,i,j-1)))*dt/dV

        u(i,j) = ( h(i,j)*u(i,j)&
                  - (DELTA-1.0d0)*(dy*(F(2,i,j)-F(2,i-1,j))&
                    + DX*(G(2,i,j)-G(2,i,j-1)))*dt/dV &
                  - (2.0d0 - DELTA)*(dy*(Fstar(2,i,j)-Fstar(2,i-1,j))&
                    + DX*(Gstar(2,i,j)-Gstar(2,i,j-1)))*dt/dV )/hnew

        v(i,j) = ( h(i,j)*v(i,j)&
                  - (DELTA-1.0d0)*(dy*(F(3,i,j)-F(3,i-1,j))&
                    + DX*(G(3,i,j)-G(3,i,j-1)))*dt/dV &
                  - (2.0d0 - DELTA)*(dy*(Fstar(3,i,j)-Fstar(3,i-1,j))&
                    + DX*(Gstar(3,i,j)-Gstar(3,i,j-1)))*dt/dV )/hnew

        C4 = ( h(i,j)*eta(i,j)&
                - (DELTA-1.0d0)*(dy*(F(4,i,j)-F(4,i-1,j))&
                  + DX*(G(4,i,j)-G(4,i,j-1)))*dt/dV &
                - (2.0d0 - DELTA)*(dy*(Fstar(4,i,j)-Fstar(4,i-1,j))&
                  + DX*(Gstar(4,i,j)-Gstar(4,i,j-1)))*dt/dV&
                            + (1.0d0 - DELTA)*dt*S(4,i,j) )/hnew

        C5 = ( h(i,j)*w(i,j)&
                - (DELTA-1.0d0)*(dy*(F(5,i,j)-F(5,i-1,j))&
                  + DX*(G(5,i,j)-G(5,i,j-1)))*dt/dV &
                - (2.0d0 - DELTA)*(dy*(Fstar(5,i,j)-Fstar(5,i-1,j))&
                  + DX*(Gstar(5,i,j)-Gstar(5,i,j-1)))*dt/dV&
                            + (1.0d0 - DELTA)*dt*S(5,i,j) )/hnew

        alpha = 1.0d0/( 1.0d0 + DELTA*DELTA*dt*dt*LAMBDA/(hnew*hnew) )

        eta(i,j) = alpha*( C4 + DELTA*dt*C5 + DELTA*DELTA*dt*dt*LAMBDA/hnew )

        w(i,j) = alpha*( C5 - C4*DELTA*dt*LAMBDA/(hnew*hnew)&
                           + DELTA*dt*LAMBDA/hnew )
        h(i,j) = hnew
      enddo
    enddo
  end subroutine imex_step_2

  subroutine hllc(priml,primr,F)
    implicit none
    real(dp), intent(in) :: priml(NEQS), primr(NEQS)
    real(dp), intent(out) :: F(NEQS)
    real(dp) :: rhol, rhor, ul, ur, etal, etar, pl, pr, al, ar
    real(dp) :: sl, sr, sl1, sl2, sr1, sr2, smid
    real(dp) :: rhostarl, rhostarr
    real(dp), dimension(NEQS) :: consl, consr, Fl, Fr
    real(dp), dimension(NEQS) :: consstarl, consstarr, Fstarl, Fstarr
    integer :: k

    rhol = priml(1); rhor = primr(1)
    ul   = priml(2); ur   = primr(2)
    etal = priml(4); etar = primr(4)

    pl = pressure(rhol, etal); pr = pressure(rhor, etar)
    al = sound_speed(rhol, etal); ar = sound_speed(rhor, etar)

    ! Wave speed estimation
    ! Davis S.F.(1988), 'Simplified second order Godunov type methods',
    ! SIAM J. Sci. and Stat. Comput., N3, 445-473.
    sl1=ul-al; sl2=ur-ar
    sr1=ul+al; sr2=ur+ar
    sl=DMIN1(sl1,sl2); sr=DMAX1(sr1,sr2)

    consl(1) = rhol; consr(1) = rhor
    do k=2, NEQS
      consl(k) = rhol*priml(k)
      consr(k) = rhor*primr(k)
    enddo

    Fl(:) = consl(:)*ul
    Fr(:) = consr(:)*ur
    Fl(2) = Fl(2) + pl
    Fr(2) = Fr(2) + pr

    ! Smid
    ! Massoni, J. (1999). Un Modèle Micromécanique pour l'Initiation par
    ! Choc et la Transition vers la Détonation dans les Matériaux Solides
    ! Hautement Energétiques (Thèse).
    smid = (sr*Fr(1)-sl*Fl(1)+Fl(2)-Fr(2))/(sr*rhor-sl*rhol+Fl(1)-Fr(1))

    rhostarl = rhol*(sl-ul)/(sl-smid)
    rhostarr = rhor*(sr-ur)/(sr-smid)
    consstarl(1) = rhostarl
    consstarr(1) = rhostarr
    consstarl(2) = rhostarl*smid
    consstarr(2) = rhostarr*smid
    do k=3, NEQS
      consstarl(k) = rhostarl*priml(k)
      consstarr(k) = rhostarr*primr(k)
    enddo

    Fstarl(:)=Fl(:)+sl*(consstarl(:)-consl(:))
    Fstarr(:)=Fr(:)+sr*(consstarr(:)-consr(:))

    if (0.0d0<sl) then
      F(:) = Fl(:)
    elseif (sr<0.0d0) then
      F(:) = Fr(:)
    elseif(0.0d0<=smid) then
      F(:)=Fstarl(:)
    else
      F(:)=Fstarr(:)
    endif

  end subroutine hllc

end module methods
