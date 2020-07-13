program extsgn_imex2d
  use parameters
  use aux
  use methods
  use model
  implicit none
  integer :: i,j,k
  real :: cpu_time_start,cpu_time_finish
  real(kind=DP) :: time = 0.0d0, cmax = 0.0d0
  real(kind=DP), allocatable :: x(:), y(:)
  real(kind=DP), allocatable :: prim(:,:,:), cons(:,:,:), sources(:,:,:)
  real(kind=DP), allocatable :: Fflux(:,:,:), Gflux(:,:,:)
  real(kind=DP) :: priml(NEQS), primr(NEQS)
  character(len=20) :: fmt

  allocate(x(0:NX+1), y(0:NY+1))
  allocate(prim(NEQS, 0:NX+1, 0:NY+1), cons(NEQS, 0:NX+1, 0:NY+1))
  allocate(Fflux(NEQS, 0:NX+1, 0:NY+1), Gflux(NEQS, 0:NX+1, 0:NY+1))
  allocate(sources(NEQS, 0:NX+1, 0:NY+1))

  call set_mesh(x,y)
  call set_ic_rpx(x,prim)
  call prim_to_cons(prim,cons)
  call output_solution(OUTPUT_FILE,x,y,prim,time)

  call riemann_fluxes_x(prim,Fflux,cmax)

  data priml /1,2,3,4,5/
  write(fmt,"(A1,I0,A7)") '(',NEQS,'(F7.3))'
  write(*,fmt) (Fflux(k,3,1), k=1,NEQS)

  deallocate(x,y,prim)

end program extsgn_imex2d
