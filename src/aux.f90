module aux
use parameters
implicit none
contains

  subroutine print_percentage(time,t1,milestone)
    implicit none
    real(dp), intent(in) :: time
    real(dp), intent(inout) :: t1
    real(dp) :: milestone
    real(dp) :: t2, percentage

    percentage = time*100.0d0/TFIN
    if (percentage > milestone) then
      call cpu_time(t2)
      write(*,'(I4,A,(F8.1),A)') NINT(percentage), ' % completed ', t2-t1, &
        ' sec'
      milestone = milestone + PERC_FREQ
      call cpu_time(t1)
    endif

    return
  end

  subroutine output_dat(x,y,prim,time)
    implicit none
    real(dp), intent(in) :: x(0:NX+1), y(0:NY+1), prim(NEQS,0:NX+1,0:NY+1)
    real(dp), intent(in) :: time
    integer :: i,k,j

    open(unit=10,file=FILE_DAT)
    do i=1, NX
      do j=1, NY
        write(10,FMT=DATA_FORMAT) x(i), y(j), (prim(k,i,j), k=1,NEQS), time
      enddo
    enddo
    close(10)
    return
  end subroutine output_dat

  subroutine output_vtk(h,u,v)
  implicit none
  real(dp), dimension(0:NX+1,0:NY+1), intent(in) :: h,u,v!,eta,w,p
  real(dp) :: x_0,y_0
  integer :: i,j,unit

  x_0 = 0.5d0*(XL + XR)
  y_0 = 0.5d0*(YL + YR)

  unit = 10
  open(unit,file=FILE_VTK)
    write(unit,'(''# vtk DataFile Version 2.0'')')
    write(unit,'(''Rectilinear 3D Dataset'')')
    write(unit,'(''ASCII'')')
    write(unit,'('' '')')
    write(unit,'(''DATASET STRUCTURED_POINTS'')')
    write(unit,FMT='(''DIMENSIONS'',I8,I8,I8)') NX/VTK_STEP+1, NY/VTK_STEP+1, 1
    write(unit,FMT='(''ORIGIN '',3(E11.4,1x))') x_0, y_0, 0.d0
    write(unit,FMT='(''SPACING'',3(E11.4,1x))') DFLOAT(VTK_STEP)*DX,&
                                                DFLOAT(VTK_STEP)*DY,0.0001d0
    write(unit,*) ' '
    write(unit,FMT='(''CELL_DATA '',I9)') NX/VTK_STEP*NY/VTK_STEP*1
    write(unit,*) ' '

    write(unit,FMT='(''SCALARS '',A6, '' float 1'')') 'depth'
    write(unit,'(''LOOKUP_TABLE default'')')
    do j=1,NY,VTK_STEP
      do i=1,NX,VTK_STEP
        write(unit,'(G11.4)') h(i,j)
      enddo
    enddo

    write(unit,FMT='(''VECTORS '',A13, '' float'')') 'velocity(m/s)'
    do j=1,NY,VTK_STEP
      do i=1,NX,VTK_STEP
        write(unit,'(3(E11.4,1x))') u(i,j), v(i,j), 0.d0
      enddo
    enddo
  close(unit)

  return
end

  subroutine output_single_prim_matrix(filename,array)
    implicit none
    character(len=NAMELEN), intent(in) :: filename
    real(dp), intent(in) :: array(0:NX+1,0:NY+1)
    integer :: i, j

    open(unit=10,file=filename)
    do j=0, NY+1
      write(10,*) (array(i,j), i=0,NX+1)
    enddo
    close(10)
    return
  end subroutine output_single_prim_matrix

  subroutine print_output_message(it,time,t1,t2)
    implicit none
    integer, intent(in) :: it
    real, intent(in) :: t1,t2
    real(dp), intent(in) :: time

    print*, ''
    write(*,*) 'NX', NX
    write(*,*) 'NY', NY
    write(*,*) 'time', time
    write(*,*) 'it', it
    write(*,*) 'DX', DX
    write(*,*) 'dt', dt
    write(*,*) 'cpu_time', t2-t1
    print*, ''
    return
  end subroutine print_output_message

end module aux
