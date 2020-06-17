module aux
use parameters
implicit none
contains

  subroutine output_solution(filename,x,y,prim,time)
    use parameters
    implicit none
    character(len=7), intent(in) :: filename
    real(kind=DP), intent(in) :: x(0:NX+1), y(0:NY+1), prim(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), intent(in) :: time
    integer :: i,k,j

    open(unit=10,file=filename)
    do i=1, NY
      do j=1, NX
        write(10,*) x(i), y(j), (prim(k,i,j), k=1,NEQS), time
      enddo
    enddo
    close(10)
    return
  end subroutine output_solution

end module aux
