module read_db

  integer*4 :: nx, ny, ny2
  real*8, allocatable :: u(:,:), v(:,:), tke(:,:), x(:),y(:)
  real*4, allocatable :: Ctlist(:), TIlist(:)
  integer*4 :: nct, nti
  CHARACTER(LEN=512) :: datapath

contains 

!********************************************************************************************

subroutine alloc_arrays()

  implicit none
  integer :: n, ierr
  logical :: iexist
  real*4 :: rdum
  character*200 :: filename
  character*200 :: cwd

  ny2=ny/2  !half height (size of axi data)

  allocate (U(nx,ny))
  allocate (V(nx,ny))
  allocate (tke(nx,ny))
  allocate (x(nx))
  allocate (y(ny))

  ! count the number of lines in Ct_set.txt
  cwd = TRIM(datapath)
  filename = trim(cwd)//'Ct_set.txt'
  inquire(file=trim(filename),EXIST=iexist)
  if ( .not. iexist ) then
    write(*,*) 'ERROR: '//trim(filename)//' not found'
    stop
  endif

  open(unit=12345,file=trim(filename))
  nct=0
  do
    read(12345,*,iostat=ierr) rdum
    ! break this loop when no valid real is read
    if ( ierr /= 0 ) then
      exit
    endif
    ! if no read error, increment line count
    nct=nct+1
  enddo
  allocate(CTlist(nct))
  rewind(12345)
  do n = 1,nct
    read(12345,*) CTlist(n)
    ! write(*,*) 'here 3',CTlist(n)
  enddo
  close(12345)

  ! count the number of lines in TI_set.txt
  cwd = TRIM(datapath)
  filename = trim(cwd)//'TI_set.txt'
  inquire(file=trim(filename),EXIST=iexist)
  if ( .not. iexist ) then
    write(*,*) 'ERROR: '//trim(filename)//' not found'
    stop
  endif

  open(unit=12345,file=trim(filename))
  nti=0
  do
    read(12345,*,iostat=ierr) rdum
    ! break this loop when no valid real is read
    if ( ierr /= 0 ) then
      exit
    endif
    ! if no read error, increment line count
    nti=nti+1
  enddo
  allocate(TIlist(nti))
  rewind(12345)
  do n = 1,nti
    read(12345,*) TIlist(n)
  enddo
  close(12345)

end subroutine alloc_arrays

!********************************************************************************************

subroutine dealloc_arrays()

  deallocate (U,V,tke,x,y)
  deallocate (Ctlist, TIlist)

end subroutine dealloc_arrays

!********************************************************************************************

subroutine read_u_v_tke(CT,TI)

  character*200 :: prefix, infile
  character*200 :: cwd

  real*4 :: CT, TI
  real :: tmpCT, tmpTI
  character*80 :: tagCT,tagTI

  integer :: i,j,k, n, nn

  integer*4 :: idum1, idum2

  ! bin ct an TI between list values
  if ( ( ct > CTlist(nct) ) .or.  ( ct < CTlist(1) ) ) then
    write(*,*) 'Ct out of bounds', ct, CTlist(1), CTlist(nct)
    stop
  endif
  if ( ( ti > TIlist(nti) ) .or.  ( ti < TIlist(1) ) ) then
    write(*,*) 'TI out of bounds', ti, TIlist(1), TIlist(nti)
    stop
  endif

  ! Read requested dataset
  write(tagCT,'(f3.1)') CT
  write(tagTI,'(f5.3)') TI
  ! write(*,*) 'Using: Ct_'//trim(tagCT)//'_TI_'//trim(tagTI)

  ! read in the data file corresponding to this point
  cwd = TRIM(datapath)
  prefix = trim(cwd)//'Ct_'//trim(tagCT)//'_TI_'//trim(tagTI)
  infile = trim(prefix)//'.xyz'

  open (unit=59,file=trim(infile),form='unformatted')
  !read(59) nx,ny
  read(59) idum1,idum2   ! nx, ny (already known)
  read(59) x(1:nx)
  !read(59) y(1:ny2)
  read(59) y(ny2+1:ny) ! read top half of y
  close(59)

  infile = trim(prefix)//'.U'
  open (unit=59,file=trim(infile),form='unformatted')
  read(59) u(1:nx,ny2+1:ny)  ! axi-symmetric data is the "top" half of y plane
  close(59)

  infile = trim(prefix)//'.V'
  open (unit=59,file=trim(infile),form='unformatted')
  read(59) v(1:nx,ny2+1:ny)  ! axi-symmetric data is the "top" half of y plane
  close(59)

  infile = trim(prefix)//'.K'
  open (unit=59,file=trim(infile),form='unformatted')
  read(59) tke(1:nx,ny2+1:ny)  ! axi-symmetric data is the "top" half of y plane
  close(59)

  ! flip the axi data to "bottom" half of y plane
  !y(ny2+1:ny) = y(1:ny2)
  !y(1:ny2)=-y(1:ny2)
  do j=1,ny2
    y(j)=-y(ny-j+1)
    do i=1,nx
      u(i,j)=u(i,ny-j+1)
      v(i,j)=v(i,ny-j+1)
      tke(i,j)=tke(i,ny-j+1)
    enddo
  enddo

end subroutine read_u_v_tke

end module read_db