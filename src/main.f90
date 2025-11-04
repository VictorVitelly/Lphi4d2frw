program main

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  use measurements
  implicit none

  call cpu_time(starting)

  !Thermalization history and autocorretion functions
  !call thermalize(-1._dp)
  !call test(-1.0_dp)

  !Measure action, magnetization, susceptibility and heat cap.
  call vary_m0(-0.5_dp,-2.5_dp,51)

  !call test2(0._dp,-3.0_dp,11)
  
  !call correlate(-3.0_dp,0._dp,11)

  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"

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
  
subroutine vary_m0(mi,mf,Nps)
  real(dp),intent(in) :: mi,mf
  integer(i4),intent(in) :: Nps
  real(dp) :: m0,dphi,AR
  real(dp) :: M(Nmsrs2), M_ave, M_err
  real(dp), allocatable :: phi(:,:)
  real(dp), allocatable :: phi_xave(:),phi_ave(:,:),phi_res(:,:),phi_err(:,:)
  real(dp), allocatable :: phi_yave(:)
  integer(i4) :: i,j,k1,k2
  open(10, file = 'data/magnet.dat', status = 'replace')
  open(20, file = 'data/susceptibility.dat', status = 'replace')
  allocate(phi_res(Nps,Lt))
  allocate(phi_err(Nps,Lt))

  do k2=1,Nps
    M=0._dp
    write(*,*) k2
    !Initializations
    m0=mi+(mf-mi)*real(k2-1,dp)/real(Nps-1,dp)
    dphi=0.4_dp +m0/30._dp
    allocate(phi(Lt+1,Lx))
    allocate(phi_xave(Lt))
    allocate(phi_yave(Lt))
    allocate(phi_ave(Nmsrs2,Lt))
    !call hot_start(phi,hotphi)
    call cold_start(phi)
    
    !Thermalization
    do i=1,thermalization
      call metropolis(m0,dphi,phi)
      !call flip_sign(phi,j)
    end do
    
    !Measurements
    do i=1,Nmsrs2
      phi_xave=0._dp
      phi_yave=0._dp
      do j=1,Nmsrs
        call flip_sign(phi)
        do k1=1,eachsweep
          call metropolis(m0,dphi,phi)
        end do
        M(i)=M(i)+abs(meanphi(phi))
        do k1=1,Lt
          phi_xave(k1)=phi_xave(k1)+meanphi_t(phi,k1)**2
          phi_yave(k1)=phi_yave(k1)+abs(meanphi_t(phi,k1))
        end do
        !call flip_sign(phi,j)
      end do
      M(i)=M(i)/real(Nmsrs,dp)
      phi_ave(i,:)=(phi_xave(:)-(phi_yave(:)**2)/real(Nmsrs,dp))/real(Nmsrs,dp)
    end do
    
    !Averages
    call mean_scalar(M,M_ave,M_err)
    do i=1,Lt
      call mean_scalar(phi_ave(:,i),phi_res(k2,i),phi_err(k2,i))
    end do
    phi_res(k2,:)=phi_res(k2,:)/real(Lx,dp)
    phi_err(k2,:)=phi_err(k2,:)/real(Lx,dp)
    
    write(10,*) m0,M_ave/real(Lx*Lt,dp),M_err/real(Lx*Lt,dp)
    write(20,*) m0, phi_res(k2,1), phi_err(k2,1), phi_res(k2,2), phi_err(k2,2), &
                &phi_res(k2,3), phi_err(k2,3), phi_res(k2,4), phi_err(k2,5), &
                &phi_res(k2,5), phi_err(k2,5), phi_res(k2,6), phi_err(k2,6), &
                &phi_res(k2,7), phi_err(k2,7), phi_res(k2,8), phi_err(k2,8), &
                &phi_res(k2,9), phi_err(k2,9), phi_res(k2,10), phi_err(k2,10)
    deallocate(phi,phi_xave,phi_yave,phi_ave)
  end do
  close(10)
  close(20)
  deallocate(phi_res,phi_err)
  end subroutine vary_m0

  
  subroutine test(m0)
  real(dp),intent(in) :: m0
  real(dp), allocatable :: phi(:,:)
  real(dp), allocatable :: phi_xave(:),phi_ave(:,:),phi_res(:),phi_err(:)
  integer(i4) :: i,j,k1
  real(dp) :: dphi,AR,ARp(Nmsrs2),AR_ave,AR_err
  open(10, file = 'data/history.dat', status = 'replace')
  open(20, file = 'data/test.dat', status = 'replace')
  allocate(phi(Lt+1,Lx))
  allocate(phi_xave(Lt))
  allocate(phi_ave(Nmsrs2,Lt))
  allocate(phi_res(Lt))
  allocate(phi_err(Lt))
  !call hot_start(phi,hotphi)
  call cold_start(phi)
  dphi=0.65_dp!+0.1_dp*(m0)

  do i=1,thermalization
    !call metropolis(m0,phi)
    !call montecarlo(m0,dphi,phi,AR)
    call montecarlopbc(m0,dphi,phi,AR)
    if( mod(i,100)==0) then
      write(10,*) i, S(m0,phi)/real(Lt*Lx,dp)
    end if
    !call flip_sign(phi,j)
  end do
  ARp=0._dp
  do i=1,Nmsrs2
    write(*,*) i
    phi_xave=0._dp
    do j=1,Nmsrs
      do k1=1,eachsweep
        !call montecarlo(m0,dphi,phi,AR)
        call montecarlopbc(m0,dphi,phi,AR)
      end do
      ARp(i)=ARp(i)+AR
      do k1=1,Lt
        phi_xave(k1)=phi_xave(k1)+meanphi_t(phi,k1)
      end do
      !call flip_sign(phi,j)
    end do
    phi_ave(i,:)=phi_xave(:)/real(Nmsrs,dp)
    ARp(i)=ARp(i)/real(Nmsrs,dp)
  end do
  call mean_scalar(ARp,AR_ave,AR_err)
  write(*,*) 'Acc. Rate=', AR_ave, AR_err
  do i=1,Lt
    call mean_scalar(phi_ave(:,i),phi_res(i),phi_err(i))
    write(20,*) at*real(i,dp),phi_res(i),phi_err(i)
  end do
  close(10)
  close(20)
  deallocate(phi,phi_xave,phi_ave,phi_res,phi_err)
  end subroutine test
  
  subroutine test2(mi,mf,Nps)
  real(dp),intent(in) :: mi,mf
  integer(i4),intent(in) :: Nps
  real(dp) :: m0,AR,dphi,ARp(Nmsrs2),AR_ave,AR_err
  real(dp) :: M(Nmsrs2), M_ave,M_err
  real(dp), allocatable :: phi(:,:)
  real(dp), allocatable :: phi_xave(:),phi_ave(:,:),phi_res(:,:),phi_err(:,:)
  integer(i4) :: i,j,k1,k2
  open(10, file = 'data/magnet.dat', status = 'replace')
  open(20, file = 'data/tvsphi.dat', status = 'replace')
  open(30, file = 'data/m02vsphi.dat', status = 'replace')
  allocate(phi_res(Nps,Lt))
  allocate(phi_err(Nps,Lt))

  do k2=1,Nps
    ARp=0._dp
    M=0._dp
    write(*,*) k2
    !Initializations
    m0=mi+(mf-mi)*real(k2-1,dp)/real(Nps-1,dp)
    dphi=0.45_dp +m0/30._dp
    allocate(phi(Lt+1,Lx))
    allocate(phi_xave(Lt))
    allocate(phi_ave(Nmsrs2,Lt))
    !call hot_start(phi,hotphi)
    call cold_start(phi)
    
    !Thermalization
    do i=1,thermalization
      call montecarlo(m0,dphi,phi,AR)
      !call montecarlopbc(m0,dphi,phi,AR)
      !call metropolis(m0,phi)
      !call flip_sign(phi,j)
    end do
    
    !Measurements
    do i=1,Nmsrs2
      phi_xave=0._dp
      do j=1,Nmsrs
        call flip_sign(phi)
        do k1=1,eachsweep
          call montecarlo(m0,dphi,phi,AR)
          !call montecarlopbc(m0,dphi,phi,AR)
        end do
        M(i)=M(i)+abs(meanphi(phi))
        ARp(i)=ARp(i)+AR
        do k1=1,Lt
          phi_xave(k1)=phi_xave(k1)+abs(meanphi_t(phi,k1))
        end do
        !call flip_sign(phi,j)
      end do
      M(i)=M(i)/real(Nmsrs,dp)
      ARp(i)=ARp(i)/real(Nmsrs,dp)
      phi_ave(i,:)=phi_xave(:)/real(Nmsrs,dp)
    end do
    
    !Averages
    call mean_scalar(M,M_ave,M_err)
    call mean_scalar(ARp,AR_ave,AR_err)
    do i=1,Lt
      call mean_scalar(phi_ave(:,i),phi_res(k2,i),phi_err(k2,i))
    end do
    phi_res(k2,:)=phi_res(k2,:)/real(Lx,dp)
    phi_err(k2,:)=phi_err(k2,:)/real(Lx,dp)
    
    write(10,*) m0,M_ave/real(Lx*Lt,dp),M_err/real(Lx*Lt,dp)
    write(30,*) m0, phi_res(k2,1), phi_err(k2,1), phi_res(k2,2), phi_err(k2,2), &
                &phi_res(k2,3), phi_err(k2,3), phi_res(k2,4), phi_err(k2,5), &
                &phi_res(k2,5), phi_err(k2,5), phi_res(k2,6), phi_err(k2,6), &
                &phi_res(k2,7), phi_err(k2,7), phi_res(k2,8), phi_err(k2,8), &
                &phi_res(k2,9), phi_err(k2,9), phi_res(k2,10), phi_err(k2,10)
    write(*,*) 'The acc. rate for m02=', m0, 'is', AR_ave, '+-', AR_err
    deallocate(phi,phi_xave,phi_ave)
  end do

  do i=1,Lt
    write(20,*) at*real(i-1,dp), phi_res(:,i), phi_err(:,i)
  end do
  close(10)
  close(20)
  close(30)
  deallocate(phi_res,phi_err)
  end subroutine test2
  

end program main
