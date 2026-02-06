module measurements
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  use statistics
  implicit none

contains

  subroutine thermalize(m0)
  real(dp), intent(in) :: m0
  real(dp), allocatable :: phi(:,:)
  real(dp) :: dphi
  integer(i4) :: i
  open(10, file = 'data/history.dat', status = 'replace')
  allocate(phi(Lt+1,Lx))
  !call hot_start(phi,hotphi)
  call cold_start(phi)
  do i=1,thermalization
    !50 sweeps for L=8, 500 sweeps for L=64
    call metropolis(m0,dphi,phi)
    if( mod(i,250)==0) then
      write(10,*) i,",", S(m0,phi)/real(Lt*Lx,dp)
    end if
  end do
  close(10)
  deallocate(phi)
  end subroutine thermalize

  subroutine correlation(phi,corr1,corr2)
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    real(dp), dimension(Lt), intent(inout) :: corr1
    real(dp), dimension(Lt,Lt), intent(inout) :: corr2
    real(dp), dimension(Lt) :: varphi
    real(dp) :: xx
    integer(i4) :: i1,i2
    varphi=0._dp
    do i1=1,Lt
      do i2=1,Lx
        varphi(i1)=varphi(i1)+phi(i1,i2)
      end do
    end do
    varphi(:)=varphi(:)/real(Lx,dp)
    xx=abs(meanphi(phi))/real(Lx*Lt,dp)
    do i1=1,Lt
      corr1(i1)=corr1(i1)+xx
      do i2=1,Lx
        corr2(i1,i2)=corr2(i1,i2)+(varphi(i1)*varphi(i2))
        !corr2(i1,i2)=corr2(i1,i2)+(phi(i1,1)*phi(i2,1))
      end do
    end do
  end subroutine correlation
  
  subroutine correlation2(phi,corr1,corr2,t)
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    integer(i4), intent(in) :: t
    real(dp), dimension(Lt), intent(inout) :: corr1
    real(dp), dimension(Lt,Lt), intent(inout) :: corr2
    real(dp) :: xx
    integer(i4) :: i1,i2
    !xx=abs(mean(phi))/real(N**2,dp)
    xx=abs(meanphi_t(phi,t))/real(Lx,dp)
    do i1=1,Lx
      corr1(i1)=corr1(i1)+xx
      do i2=1,Lx
        corr2(i1,i2)=corr2(i1,i2)+(phi(t,i1)*phi(t,i2))
      end do
    end do
  end subroutine correlation2

  subroutine histogram(x,A1,A2)
  real(dp), dimension(:,:), intent(in) :: x
  integer(i4), dimension(bins,Lt), intent(inout) :: A1
  real(dp), dimension(bins), intent(in) :: A2
  integer(i4) :: i,j,k
  do i=1,bins
    do j=1,Lt
      do k=1,size(x,dim=2)
        if(x(j,k) .le. real(A2(i),dp)+binwidth/2._dp .and. x(j,j)>real(A2(i),dp)-binwidth/2._dp ) then
          A1(i,j)=A1(i,j)+1
          cycle
        end if
      end do
    end do
  end do
  end subroutine histogram

  subroutine make_histogram(m0)
    real(dp), intent(in) :: m0
    real(dp) :: dphi,norm(Lt),AR(Lt),ARp(Lt,Nmsrs2),AR_ave(Lt),AR_err(Lt)
    integer(i4) :: i,j,k
    real(dp), allocatable :: phi(:,:), A2(:)
    integer(i4), allocatable :: A1(:,:)
    open(80, file = 'data/histogram.dat', status = 'replace')
    allocate(phi(Lt+1,Lx))
    allocate(A1(bins,Lt))
    allocate(A2(bins))
    call cold_start(phi)
    ARp=0._dp
    do i=1,bins
      A2(i)=minn+binwidth/2._dp+real(i-1,dp)*binwidth
    end do
    A1=0
    dphi=0.5_dp+m0/20._dp
    do i=1,thermalization
      call montecarlo(m0,dphi,phi,AR)
    end do

    do i=1,Nmsrs2
      do j=1,Nmsrs
        call flip_sign(phi)
        do k=1,eachsweep
          call montecarlo(m0,dphi,phi,AR)
        end do
        ARp(:,i)=ARp(:,i)+AR(:)
        call histogram(phi,A1,A2)
      end do
      ARp(:,i)=ARp(:,i)/real(Nmsrs,dp)
    end do

    do i=1,Lt
      call mean_scalar(ARp(i,:),AR_ave(i),AR_err(i))
      write(*,*) i, AR_ave(i), AR_err(i)
    end do

    norm(:)=0._dp
    do i=1,bins
      norm(:)=norm(:)+A1(i,:)
    end do
    norm(:)=norm(:)*(real(maxx-minn,dp) )/real(bins,dp)

    do i=1,bins
      do j=1,Lt
        write(80,*) j, A2(i), A1(i,j), sqrt( real(A1(i,j),dp) )
      end do
    end do
    deallocate(A1,A2)
    close(80)
  end subroutine make_histogram



end module measurements
