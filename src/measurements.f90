module measurements
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  use statistics
  implicit none

contains

  subroutine correlation(phi,corr1,corr2,dphi)
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    real(dp), intent(in) :: dphi
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
    xx=abs(varphi(1))
    do i1=1,Lt
      corr1(i1)=corr1(i1)+xx
      do i2=1,Lt
        corr2(i1,i2)=corr2(i1,i2)+(varphi(i1)*varphi(i2))
        !corr2(i1,i2)=corr2(i1,i2)+(phi(i1,1)*phi(i2,1))
      end do
    end do
  end subroutine correlation
  
  subroutine correlation2(phi,corr1,corr2,t,dphi)
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    integer(i4), intent(in) :: t
    real(dp), intent(in) :: dphi
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
  
  subroutine correlate(mi,mf,Nts)
  real(dp), intent(in) :: mi,mf
  !subroutine correlate(m0,lambi,lambf,Nts)
  !real(dp), intent(in) :: m0,lambi,lambf
  real(dp) :: x0,dphi
  integer(i4) :: Nts
  real(dp), allocatable :: phi(:,:),corr1(:),corr2(:,:),CF(:,:),CF_ave(:,:),CF_delta(:,:)
  integer(i4) :: i,j,k,i2
  open(60, file = 'data/corrfunc.dat', status = 'replace')
  allocate(phi(Lt+1,Lx))
  allocate(corr1(Lt))
  allocate(corr2(Lt,Lt))
  allocate(CF(Lt,Nmsrs2))
  allocate(CF_ave(Lt,Nts))
  allocate(CF_delta(Lt,Nts))

  !call cold_start(phi)
  do k=1,Nts
    CF(:,:)=0._dp
    x0=mi+(mf-mi)*real(k-1,dp)/real(Nts-1,dp)
    dphi=0.45_dp +x0/30._dp
    !x0=lambi+(lambf-lambi)*real(k-1,dp)/real(Nts-1,dp)
    write(*,*) x0
    call cold_start(phi)
    do j=1,2*thermalization
      !call cycles(x0,lamb0,phi,4)
      call metropolis(x0,dphi,phi)
    end do
    do j=1,Nmsrs2
      corr1(:)=0._dp
      corr2(:,:)=0._dp
      do i=1,Nmsrs
        call flip_sign(phi)
        do i2=1,eachsweep
          !call cycles(x0,lamb0,phi,4)
          call metropolis(x0,dphi,phi)
        end do
        call correlation(phi,corr1,corr2,dphi)
      end do
      corr1(:)=corr1(:)/real(Nmsrs,dp)
      corr2(:,:)=corr2(:,:)/real(Nmsrs,dp)
      !write(*,*) corr2(N/2,1), corr1(1)**2, corr2(N/2,1)-corr1(1)**2
      do i=1,Lt
        CF(i,j)=corr2(i,1) !-(corr1(1)**2)
      end do
    end do
    do j=1,Lt
      call mean_scalar(CF(j,:),CF_ave(j,k),CF_delta(j,k))
    end do
  end do
  
  do k=1,Lt+1
    write(60,*) abs(k-1), CF_ave(ivx(k),:), CF_delta(ivx(k),:)
  end do

  deallocate(corr1,corr2,CF,CF_ave,CF_delta)
  deallocate(phi)
  close(60)
  end subroutine correlate
  


end module measurements
