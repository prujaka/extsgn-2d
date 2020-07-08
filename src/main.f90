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

  allocate(x(0:NX+1), y(0:NY+1))
  allocate(prim(NEQS, 0:NX+1, 0:NY+1), cons(NEQS, 0:NX+1, 0:NY+1))
  allocate(Fflux(NEQS, 0:NX+1, 0:NY+1), Gflux(NEQS, 0:NX+1, 0:NY+1))
  allocate(sources(NEQS, 0:NX+1, 0:NY+1))

  call set_mesh(x,y)
  call set_ic_rpx(x,prim)
  call prim_to_cons(prim,cons)
  call output_solution(OUTPUT_FILE,x,y,prim,time)

  call set_bc(prim)
  call output_single_prim_matrix('boundar',prim(1,:,:))

  call hllc(prim(:,1,1),prim(:,2,1),Fflux(:,1,1),cmax)
  write(*,*) (Fflux(k,1,1), k=1, NEQS)

  deallocate(x,y,prim)

end program extsgn_imex2d
