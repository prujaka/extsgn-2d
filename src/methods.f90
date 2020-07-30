module methods
  use parameters
  use model
  use aux
  implicit none
  contains

  subroutine initialize_problem(x,y,prim,cons,it,time)
    use parameters
    implicit none
    real(kind=DP), intent(out) :: x(0:NX+1), y(0:NY+1), time
    real(kind=DP), dimension(NEQS, 0:NX+1, 0:NY+1),intent(out) :: prim, cons
    integer, intent(out) :: it

    call set_mesh(x,y)
    call set_ic(x,y,prim)
    call prim_to_cons(prim,cons)
    it = 0
    time = 0.0d0

    return
  end subroutine initialize_problem

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
  end subroutine set_mesh

  subroutine set_ic(x,y,prim)
    use parameters
    implicit none
    real(kind=DP), intent(in)  :: x(0:NX+1), y(0:NY+1)
    real(kind=DP), intent(out) :: prim(NEQS, 0:NX+1, 0:NY+1)

    select case(SELECTOR_IC)
      case(IC_RP_X)
        call set_ic_rp_x(x,prim)
  		case(IC_RP_Y)
  	    call set_ic_rp_y(y,prim)
		end select
    return
  end subroutine set_ic

  subroutine set_ic_rp_x(x,prim)
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
  end subroutine set_ic_rp_x
  subroutine set_ic_rp_y(y,prim)
    use parameters
    implicit none
    real(kind=DP), intent(in)  :: y(0:NY+1)
    real(kind=DP), intent(out) :: prim(NEQS, 0:NX+1, 0:NY+1)
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

    return
  end subroutine set_ic_rp_y

  subroutine prim_to_cons(prim,cons)
    use parameters
    implicit none
    real(kind=DP), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), intent(out) :: cons(NEQS,0:NX+1,0:NY+1)
    integer :: k

    cons(1,:,:) = prim(1,:,:)
    do k=2, NEQS
      cons(k,:,:) = cons(1,:,:)*prim(k,:,:)
    enddo

    return
  end subroutine prim_to_cons

  subroutine cons_to_prim(cons,prim)
    use parameters
    implicit none
    real(kind=DP), intent(in) :: cons(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), intent(out) :: prim(NEQS,0:NX+1,0:NY+1)
    integer :: k

    prim(1,:,:) = cons(1,:,:)
    do k=2, NEQS
      prim(k,:,:) = cons(k,:,:)/prim(1,:,:)
    enddo

    return
  end subroutine cons_to_prim

  subroutine get_solution_godunov(prim,cons,it,time,cmax)
    use parameters
    use model
    implicit none
    real(kind=DP), intent(inout) :: cons(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), intent(inout) :: prim(NEQS,0:NX+1,0:NY+1)
    integer, intent(inout) :: it
    real(kind=DP), intent(inout) :: time
    real(kind=DP), intent(out) :: cmax

    do while(time<TIMEFINAL)
  		if (it.ge.ITFINAL) exit
      call timestep_firstorder(prim,cons,cmax)
  		it=it+1
  		time=time+dt
  	enddo

    return
  end subroutine get_solution_godunov

  subroutine timestep_firstorder(prim,cons,cmax)
    use parameters
    use model
    implicit none
    real(kind=DP), intent(inout) :: cons(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), intent(inout) :: prim(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), intent(out) :: cmax
    real(kind=DP), dimension(NEQS,0:NX+1,0:NY+1) :: Fflux,Gflux,S
    real(kind=DP), dimension(0:NX+1,0:NY+1) :: h,u,v,eta,w
    real(kind=DP) :: cmax1, cmax2

    call set_bc(prim)
    ! call riemann_fluxes_x(prim,Fflux,cmax1)
    ! call riemann_fluxes_y(prim,Gflux,cmax2)
    call muscl_step_x(prim,Fflux,cmax1)
    Gflux = 0.0d0
    cmax = cmax1
    ! cmax = dmax1(cmax1,cmax2)
    dt=CFL*DL/cmax
    call godunov(cons,Fflux,Gflux)
    call get_prims_gn(prim,h,u,v,eta,w)
    call ode_exact_solution(h,u,v,eta,w)
    call make_sources(h,eta,w,S)
    call ode_euler_step(S,cons)
    call cons_to_prim(cons,prim)

    return
  end subroutine timestep_firstorder

  subroutine set_bc(prim)
    use parameters
    implicit none
    real(kind=DP), intent(inout) :: prim(NEQS,0:NX+1,0:NY+1)
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


    return
  end subroutine set_bc

  subroutine riemann_fluxes_x(prim,flux,cmax)
    use parameters
    implicit none
    real(kind=DP), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), intent(out) :: flux(NEQS,0:NX+1,0:NY+1), cmax
    real(kind=DP), dimension(NEQS) :: priml, primr, F
    integer :: i,j

    cmax = 0.0d0
    do j=1,NY
      do i=0,Nx
        priml = prim(:,i,j)
        primr = prim(:,i+1,j)
        call hllc(priml,primr,F,cmax)
        flux(:,i,j) = F
      enddo
    enddo

    return
  end subroutine riemann_fluxes_x
  subroutine riemann_fluxes_y(prim,flux,cmax)
    use parameters
    implicit none
    real(kind=DP), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), intent(out) :: flux(NEQS,0:NX+1,0:NY+1), cmax
    real(kind=DP), dimension(NEQS) :: priml, primr, F
    integer :: i,j

    cmax = 0.0d0
    do i=1,NX
      do j=0,NY
        priml = prim(:,i,j)
        priml(2) = prim(3,i,j)
        priml(3) = prim(2,i,j)

        primr = prim(:,i,j+1)
        primr(2) = prim(3,i,j+1)
        primr(3) = prim(2,i,j+1)

        call hllc(priml,primr,F,cmax)

        flux(:,i,j) = F
        flux(2,i,j) = F(3)
        flux(3,i,j) = F(2)
      enddo
    enddo

    return
  end subroutine riemann_fluxes_y

  subroutine godunov(cons,Fflux,Gflux)
    use parameters
    implicit none

    real(kind=DP), dimension(NEQS,0:NX+1,0:NY+1), intent(in) :: Fflux,Gflux
    real(kind=DP), dimension(NEQS,0:NX+1,0:NY+1), intent(inout) :: cons
    integer :: i,j

    do j=1,Ny
      do i=1,Nx
        cons(:,i,j)=cons(:,i,j)-dt/dV*( dy*(Fflux(:,i,j)-Fflux(:,i-1,j))&
                                       +dx*(Gflux(:,i,j)-Gflux(:,i,j-1)) )
      enddo
    enddo
    return
  end subroutine godunov

	subroutine ode_euler_step(S,cons)
		use parameters
		implicit none
		real(kind=DP), intent(in) :: S(NEQS,0:NX+1,0:NY+1)
		real(kind=DP), intent(inout) :: cons(NEQS,0:NX+1,0:NY+1)

		cons(:,:,:) = cons(:,:,:) + dt * S(:,:,:)
		return
	end subroutine ode_euler_step

  subroutine muscl_step_x(prim,F,cmax)
		use parameters
    use aux
    real(kind=DP), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), intent(out) :: F(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), intent(out) :: cmax
    real(kind=DP), dimension(NEQS,0:NX+1,0:NY+1) :: primr,priml

    call data_reconstruction_x(prim,primr,priml)
    call set_bc_muscl_x(primr,priml)
    call riemann_fluxes_muscl_x(primr,priml,F,cmax)

		return
	end subroutine muscl_step_x

	subroutine data_reconstruction_x(prim,primr,priml)
		use parameters
		implicit none
		real(kind=DP), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
		real(kind=DP), dimension(NEQS,0:NX+1,0:NY+1), intent(out) :: primr,priml
    real(kind=DP) :: slope(NEQS,0:NX+1,0:NY+1)

    call get_slopes_x(prim,slope)
    primr(:,:,:) = prim(:,:,:) - 0.5d0*slope(:,:,:)
    priml(:,:,:) = prim(:,:,:) + 0.5d0*slope(:,:,:)
    return
	end subroutine data_reconstruction_x

	subroutine get_slopes_x(prim,slope)
		use parameters
		implicit none
    real(kind=DP), intent(in) :: prim(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), intent(out) :: slope(NEQS,0:NX+1,0:NY+1)
    integer :: i

    do i=1,NX
      slope(:,i,:) = 0.5d0*(1.0d0+OMEGA)*(prim(:,i,:) - prim(:,i-1,:))&
                   + 0.5d0*(1.0d0-OMEGA)*(prim(:,i+1,:) - prim(:,i,:))
    enddo
		return
	end subroutine get_slopes_x

	subroutine set_bc_muscl_x(primr,priml)
		use parameters
		implicit none
    real(kind=DP), dimension(NEQS,0:NX+1,0:NY+1), intent(inout) :: primr,priml
    integer :: j

    do j=1,NY
      priml(:,0,j) = primr(:,1,j)
      priml(2,0,j) = BC_U_LEFT*primr(2,1,j)

      primr(:,NX+1,j) = priml(:,NX,j)
      primr(2,NX+1,j) = BC_U_RIGHT*priml(2,NX,j)
    enddo

		return
	end subroutine set_bc_muscl_x

	subroutine riemann_fluxes_muscl_x(primr,priml,flux,cmax)
		use parameters
		implicit none
		real(kind=DP), dimension(NEQS,0:NX+1,0:NY+1), intent(in) :: primr,priml
		real(kind=DP), intent(out) :: flux(NEQS,0:NX+1,0:NY+1), cmax
    real(kind=DP), dimension(NEQS) :: statel, stater, F
    integer :: i,j

    cmax = 0.0d0
    do j=1,NY
      do i=0,Nx
        statel = priml(:,i,j)
        stater = primr(:,i+1,j)
        call hllc(statel,stater,F,cmax)
        flux(:,i,j) = F
      enddo
    enddo
		return
	end subroutine riemann_fluxes_muscl_x

  subroutine hllc(priml,primr,F,cmax)
    use parameters
    use model
    implicit none
    real(kind=DP), intent(in) :: priml(NEQS), primr(NEQS)
    real(kind=DP), intent(out) :: F(NEQS)
    real(kind=DP), intent(inout) :: cmax
    real(kind=DP) :: rhol, rhor, ul, ur, etal, etar, pl, pr, al, ar
    real(kind=DP) :: sl, sr, sl1, sl2, sr1, sr2, smid
    real(kind=DP) :: rhostarl, rhostarr
    real(kind=DP), dimension(NEQS) :: consl, consr, Fl, Fr
    real(kind=DP), dimension(NEQS) :: consstarl, consstarr, Fstarl, Fstarr
    integer :: k

    rhol = priml(1); rhor = primr(1)
    ul   = priml(2); ur   = primr(2)
    etal = priml(4); etar = primr(4)

    pl = get_p(rhol, etal); pr = get_p(rhor, etar)
    al = get_a(rhol, etal); ar = get_a(rhor, etar)

    ! Wave speed estimation
  	! Davis S.F.(1988), 'Simplified second order Godunov type methods',
  	! SIAM J. Sci. and Stat. Comput., N3, 445-473.
		sl1=ul-al; sl2=ur-ar
		sr1=ul+al; sr2=ur+ar
		sl=DMIN1(sl1,sl2); sr=DMAX1(sr1,sr2)
		cmax = DMAX1(cmax,sl,sr)

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

    return
  end subroutine hllc

end module methods
