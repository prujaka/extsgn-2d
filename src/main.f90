program extsgn_imex2d
  use parameters
  use aux
  use methods
  use model
  implicit none
  integer :: it
  real :: t1,t2
  real(kind=DP) :: time
  real(kind=DP), allocatable :: x(:),y(:)
  real(kind=DP), allocatable :: prim(:,:,:), cons(:,:,:)

  allocate(x(0:NX+1), y(0:NY+1))
  allocate(prim(NEQS, 0:NX+1, 0:NY+1), cons(NEQS, 0:NX+1, 0:NY+1))

  call cpu_time(t1)

  call initialize_problem(x,y,prim,cons,it,time)
  call get_solution(prim,cons,it,time)
  call output_solution(OUTPUT_FILE,x,y,prim,time)

  call cpu_time(t2)

  call print_output_message(it,time,t1,t2)

  deallocate(x,y,prim,cons)

end program extsgn_imex2d
