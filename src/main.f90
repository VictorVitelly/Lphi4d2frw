program main

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  implicit none

  call cpu_time(starting)

  !Thermalization history and autocorretion functions
  !call thermalize(-1._dp)
  !call test(-1.0_dp)
  
  call test2(0._dp,-2._dp,11)

  !Measure action, magnetization, susceptibility and heat cap.
  !call vary_m0(-1.3_dp,-1.1_dp,50)

  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"

contains

  subroutine thermalize(m0)
  real(dp), intent(in) :: m0
  real(dp), allocatable :: phi(:,:)
  integer(i4) :: i
  open(10, file = 'data/history.dat', status = 'replace')
  allocate(phi(Lt+1,Lx))
  !call hot_start(phi,hotphi)
  call cold_start(phi)
  do i=1,thermalization
    !50 sweeps for L=8, 500 sweeps for L=64
    call metropolis(m0,phi)
    if( mod(i,250)==0) then
      write(10,*) i,",", S(m0,phi)/real(Lt*Lx,dp)
    end if
  end do
  close(10)
  deallocate(phi)
  end subroutine thermalize
  
 subroutine vary_m0(mi,mf,Nts)
  real(dp), intent(in) :: mi,mf
  integer(i4), intent(in) :: Nts
  real(dp), dimension(Lt+1,Lx) :: phi
  integer(i4) :: i,j,k
  real(dp), dimension(Nmsrs2) :: E
  real(dp) :: m0,vol,norm,EE,E_ave,E_delta
  open(10, file = 'data/action.dat', status = 'replace')
  norm=real(Nmsrs,dp)
  vol=real(Lt*Lx,dp)
  do k=1,Nts
  write(*,*) k
  call cold_start(phi)
    m0=mi+(mf-mi)*real(k-1,dp)/real(Nts-1)
    E(:)=0._dp
    do j=1,Nmsrs2
      do i=1,sweeps
        if(i>thermalization .and. mod(i,eachsweep)==0 ) then
          EE=S(m0,phi)
          E(j)=E(j)+EE
        end if
        call metropolis(m0,phi)
      end do
      E(j)=E(j)/norm
    end do
    call mean_scalar(E,E_ave,E_delta)
    write(10,*) m0, E_ave/vol, E_delta/vol
  end do
  close(10)
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
  dphi=0.5_dp+0.1_dp*(m0)
  do i=1,thermalization
    !50 sweeps for L=8, 500 sweeps for L=64
    !call metropolis(m0,phi)
    call montecarlo(m0,dphi,phi,AR)
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
        call montecarlo(m0,dphi,phi,AR)
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
    write(20,*) at*real(i-1,dp),phi_res(i),phi_err(i)
  end do
  close(10)
  close(20)
  deallocate(phi,phi_xave,phi_ave,phi_res,phi_err)
  end subroutine test
  
  subroutine test2(mi,mf,Nps)
  real(dp),intent(in) :: mi,mf
  integer(i4),intent(in) :: Nps
  real(dp) :: m0,AR,dphi,ARp(Nmsrs2),AR_ave,AR_err
  real(dp), allocatable :: phi(:,:)
  real(dp), allocatable :: phi_xave(:),phi_ave(:,:),phi_res(:,:),phi_err(:,:)
  integer(i4) :: i,j,k1,k2
  open(20, file = 'data/tvsphi.dat', status = 'replace')
  open(30, file = 'data/m02vsphi.dat', status = 'replace')
  allocate(phi_res(Nps,Lt))
  allocate(phi_err(Nps,Lt))

  do k2=1,Nps
    ARp=0._dp
    write(*,*) k2
    !Initializations
    m0=mi+(mf-mi)*real(k2-1,dp)/real(Nps-1,dp)
    dphi=0.5_dp+0.1*m0
    allocate(phi(Lt+1,Lx))
    allocate(phi_xave(Lt))
    allocate(phi_ave(Nmsrs2,Lt))
    !call hot_start(phi,hotphi)
    call cold_start(phi)
    
    !Thermalization
    do i=1,thermalization
      call montecarlo(m0,dphi,phi,AR)
      !call metropolis(m0,phi)
      !call flip_sign(phi,j)
    end do
    
    !Measurements
    do i=1,Nmsrs2
      phi_xave=0._dp
      do j=1,Nmsrs
        do k1=1,eachsweep
          call montecarlo(m0,dphi,phi,AR)
        end do
        ARp(i)=ARp(i)+AR
        do k1=1,Lt
          phi_xave(k1)=phi_xave(k1)+meanphi_t(phi,k1)
        end do
        !call flip_sign(phi,j)
      end do
      ARp(i)=ARp(i)/real(Nmsrs,dp)
      phi_ave(i,:)=phi_xave(:)/real(Nmsrs,dp)
    end do
    
    !Averages
    call mean_scalar(ARp,AR_ave,AR_err)
    do i=1,Lt
      call mean_scalar(phi_ave(:,i),phi_res(k2,i),phi_err(k2,i))
    end do
    write(30,*) m0, phi_res(k2,1), phi_err(k2,1), phi_res(k2,4), phi_err(k2,4), &
                &phi_res(k2,8), phi_err(k2,8), phi_res(k2,12), phi_err(k2,12), &
                &phi_res(k2,16), phi_err(k2,16), phi_res(k2,20), phi_err(k2,20), &
                &phi_res(k2,24), phi_err(k2,24), phi_res(k2,28), phi_err(k2,28), &
                &phi_res(k2,32), phi_err(k2,32)
    write(*,*) 'The acc. rate for m02=', m0, 'is', AR_ave, '+-', AR_err
    deallocate(phi,phi_xave,phi_ave)
  end do

  do i=1,Lt
    write(20,*) at*real(i-1,dp), phi_res(:,i), phi_err(:,i)
  end do
  close(10)
  close(30)
  deallocate(phi_res,phi_err)
  end subroutine test2
  

end program main
