module arrays
    use iso_fortran_env, only : dp => real64, i4 => int32
    use parameters
    implicit none

    real(dp), allocatable :: phi(:,:), alfa(:),etac(:)
    integer(i4),dimension(:),allocatable :: ipx,imx,ipt,imt

contains

  subroutine init_vecs()
  integer(i4) :: i
  real(dp) :: x1,x2
  allocate(ipx(Lx),imx(Lx))
  allocate(alfa(Lt))
  allocate(etac(Lt))
  do i=1,Lx-1
    ipx(i)=i+1
  end do
  ipx(Lx)=1
  do i=2,Lx
    imx(i)=i-1
  end do
  imx(1)=Lx
  x1=(1._dp/alfaf-1._dp/alfai)/real(Lt-1,dp)
  x2=1._dp/alfai -x1  
  do i=1,Lt 
    alfa(i)=-1._dp/(H0*(x1*real(i,dp)+x2 ))
    etac(i)=x1*real(i,dp)+x2
  end do
  end subroutine init_vecs

end module arrays
