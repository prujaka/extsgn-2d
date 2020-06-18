program extsgn_imex2d
  use parameters
  use aux
  use methods
  implicit none
  integer :: i,j,k
  real :: cpu_time_start,cpu_time_finish
  real(kind=DP) :: time = 0.0d0
  real(kind=DP), allocatable :: x(:), y(:)
  real(kind=DP), allocatable :: prim(:,:,:), cons(:,:,:)

  allocate(x(0:NX+1), y(0:NY+1))
  allocate(prim(NEQS, 0:NX+1, 0:NY+1), cons(NEQS, 0:NX+1, 0:NY+1))

  call set_mesh(x,y)
  call set_ic_rpx(x,prim)
  call prim_to_cons(prim,cons)
  call output_solution(OUTPUT_FILE,x,y,prim,time)

  call set_bc(prim)
  call output_single_prim_matrix('boundar',prim(1,:,:))

  deallocate(x,y,prim)

end program extsgn_imex2d
