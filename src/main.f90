program extsgn_imex2d
  use parameters
  use methods
  implicit none
  integer :: i,j,k
  real :: cpu_time_start,cpu_time_finish
  real(kind=DP) :: time = 0.0d0
  real(kind=DP), allocatable :: x(:), y(:)
  real(kind=DP), allocatable :: prim(:,:,:)

  allocate(x(0:NX+1), y(0:NY+1))
  allocate(prim(NEQS, 0:NX+1, 0:NY+1))

  call set_mesh(x,y)
  call set_ic_rpx(x,prim)

  open(unit=10,file=OUTPUT_FILE)
  do i=1, NY
    do j=1, NX
      write(10,*) x(i), y(j), (prim(k,i,j), k=1,NEQS), time
    enddo
  enddo
  close(10)

  deallocate(x,y,prim)

end program extsgn_imex2d
