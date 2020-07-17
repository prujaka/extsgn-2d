module aux
use parameters
implicit none
contains

  subroutine output_solution(filename,x,y,prim,time)
    use parameters
    implicit none
    character(len=OUTPUT_FILENAME_LENGTH), intent(in) :: filename
    real(kind=DP), intent(in) :: x(0:NX+1), y(0:NY+1), prim(NEQS,0:NX+1,0:NY+1)
    real(kind=DP), intent(in) :: time
    integer :: i,k,j

    open(unit=10,file=filename)
    do j=1, NY
      do i=1, NX
        write(10,*) x(i), y(j), (prim(k,i,j), k=1,NEQS), time
      enddo
    enddo
    close(10)
    return
  end subroutine output_solution

  subroutine output_single_prim_matrix(filename,array)
    use parameters
    implicit none
    character(len=OUTPUT_FILENAME_LENGTH), intent(in) :: filename
    real(kind=DP), intent(in) :: array(0:NX+1,0:NY+1)
    integer :: i, j

    open(unit=10,file=filename)
    do j=0, NY+1
      write(10,*) (array(i,j), i=0,NX+1)
    enddo
    close(10)
    return
  end subroutine output_single_prim_matrix

  subroutine print_output_message(it,time,t1,t2,cmax)
    use parameters
    implicit none
    integer, intent(in) :: it
    real, intent(in) :: t1,t2
    real(kind=DP), intent(in) :: time,cmax

    print*, ''
    write(*,*) 'NX', NX
    write(*,*) 'Ny', Ny
    write(*,*) 'time', time
    write(*,*) 'it', it
    write(*,*) 'dt', dt
    write(*,*) 'cmax', cmax
    write(*,*) 'cpu_time', t2-t1
    return
  end subroutine print_output_message

end module aux
