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
  !call vary_m0(0._dp,-3.0_dp,11)
  !call vary_m0(-0.0_dp,-3.0_dp,11)
  !call vary_m0(0._dp,-3._dp,21)

  call make_histogram(-0.128_dp)

  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"

contains
  
  subroutine vary_m0(mi,mf,Nps)
  real(dp),intent(in) :: mi,mf
  integer(i4),intent(in) :: Nps
  real(dp) :: m0,AR(Lt),dphi,ARp(Lt,Nmsrs2),AR_ave(Lt),AR_err(Lt)
  real(dp) :: M(Nmsrs2), M_ave, M_err, E(Nmsrs2), E_ave, E_err
  real(dp), allocatable :: phi(:,:)
  real(dp), allocatable :: phi_xave(:),phi_ave(:,:),phi_res(:,:),phi_err(:,:)
  real(dp), allocatable :: Sphi_xave(:),Sphi_ave(:,:),Sphi_res(:,:),Sphi_err(:,:)
  real(dp), allocatable :: corr1(:),corr2(:,:),CF(:,:),CF_ave(:,:),CF_delta(:,:)
  integer(i4) :: i,j,k1,k2
  open(10, file = 'data/magnet.dat', status = 'replace')
  open(20, file = 'data/tvsphi.dat', status = 'replace')
  open(30, file = 'data/m02vsphi.dat', status = 'replace')
  open(40, file = 'data/action.dat', status = 'replace')  
  open(50, file = 'data/tvssphi.dat', status = 'replace')  
  open(60, file = 'data/m02vssphi.dat', status = 'replace')
  open(70, file = 'data/corrfunc.dat', status = 'replace')
  allocate(phi_res(Nps,Lt))
  allocate(phi_err(Nps,Lt))
  allocate(Sphi_res(Nps,Lt))
  allocate(Sphi_err(Nps,Lt))
  allocate(corr1(Lt))
  allocate(corr2(Lt,Lt))
  allocate(CF(Lt,Nmsrs2))
  allocate(CF_ave(Lt,Nps))
  allocate(CF_delta(Lt,Nps))

  do k2=1,Nps
    ARp=0._dp
    M=0._dp
    E=0._dp
    CF(:,:)=0._dp
    write(*,*) k2
    !Initializations
    m0=mi+(mf-mi)*real(k2-1,dp)/real(Nps-1,dp)
    dphi=0.5_dp +m0/20._dp
    allocate(phi(Lt+1,Lx))
    allocate(phi_xave(Lt))
    allocate(phi_ave(Nmsrs2,Lt))
    allocate(Sphi_xave(Lt))
    allocate(Sphi_ave(Nmsrs2,Lt))
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
      Sphi_xave=0._dp
      corr1(:)=0._dp
      corr2(:,:)=0._dp
      do j=1,Nmsrs
        call flip_sign(phi)
        do k1=1,eachsweep
          call montecarlo(m0,dphi,phi,AR)
          !call montecarlopbc(m0,dphi,phi,AR)
        end do
        M(i)=M(i)+abs(meanphi(phi))
        E(i)=E(i)+S(m0,phi)
        ARp(:,i)=ARp(:,i)+AR(:)
        call correlation(phi,corr1,corr2)
        do k1=1,Lt
          phi_xave(k1)=phi_xave(k1)+abs(meanphi_t(phi,k1))
          Sphi_xave(k1)=Sphi_xave(k1)+meanS_t(m0,phi,k1)
        end do
        !call flip_sign(phi,j)
      end do
      M(i)=M(i)/real(Nmsrs,dp)
      E(i)=E(i)/real(Nmsrs,dp)
      ARp(:,i)=ARp(:,i)/real(Nmsrs,dp)
      phi_ave(i,:)=phi_xave(:)/real(Nmsrs,dp)
      Sphi_ave(i,:)=Sphi_xave(:)/real(Nmsrs,dp)
      corr1(:)=corr1(:)/real(Nmsrs,dp)
      corr2(:,:)=corr2(:,:)/real(Nmsrs,dp)
      do j=1,Lt
        CF(j,i)=corr2(j,1) !-(corr1(1)**2)
      end do
    end do
    
    !Averages
    call mean_scalar(M,M_ave,M_err)
    call mean_scalar(E,E_ave,E_err)
    do i=1,Lt
      call mean_scalar(ARp(i,:),AR_ave(i),AR_err(i))
      call mean_scalar(phi_ave(:,i),phi_res(k2,i),phi_err(k2,i))
      call mean_scalar(Sphi_ave(:,i),Sphi_res(k2,i),Sphi_err(k2,i))
      call mean_scalar(CF(i,:),CF_ave(i,k2),CF_delta(i,k2))
    end do
    phi_res(k2,:)=phi_res(k2,:)/real(Lx,dp)
    phi_err(k2,:)=phi_err(k2,:)/real(Lx,dp)
    Sphi_res(k2,:)=Sphi_res(k2,:)/real(Lx,dp)
    Sphi_err(k2,:)=Sphi_err(k2,:)/real(Lx,dp)
    
    write(10,*) m0,M_ave/real(Lx*Lt,dp),M_err/real(Lx*Lt,dp)
    write(40,*) m0,E_ave/real(Lx*(Lt-1),dp),E_err/real(Lx*(Lt-1),dp)
    write(30,*) m0, phi_res(k2,1), phi_err(k2,1), phi_res(k2,2), phi_err(k2,2), &
                &phi_res(k2,4), phi_err(k2,4), phi_res(k2,6), phi_err(k2,6), &
                &phi_res(k2,8), phi_err(k2,8), phi_res(k2,10), phi_err(k2,10), &
                &phi_res(k2,12), phi_err(k2,12), phi_res(k2,14), phi_err(k2,14), &
                &phi_res(k2,15), phi_err(k2,15), phi_res(k2,16), phi_err(k2,16)
    write(60,*) m0, Sphi_res(k2,1), Sphi_err(k2,1), Sphi_res(k2,2), Sphi_err(k2,2), &
                &Sphi_res(k2,4), Sphi_err(k2,4), Sphi_res(k2,6), Sphi_err(k2,6), &
                &Sphi_res(k2,8), Sphi_err(k2,8), Sphi_res(k2,10), Sphi_err(k2,10), &
                &Sphi_res(k2,12), Sphi_err(k2,12), Sphi_res(k2,14), Sphi_err(k2,14), &
                &Sphi_res(k2,15), Sphi_err(k2,15), Sphi_res(k2,16), Sphi_err(k2,16)
    write(*,*) 'The acc. rate for m02=',m0
    do i=1,Lt
      write(*,*) i, AR_ave(i), AR_err(i)
    end do
    deallocate(phi,phi_xave,phi_ave)
    deallocate(Sphi_xave,Sphi_ave)
  end do

  do i=1,Lt
    write(20,*) at*real(i-1,dp), phi_res(:,i), phi_err(:,i)
    write(50,*) at*real(i-1,dp), Sphi_res(:,i), Sphi_err(:,i)
  end do
  do i=1,Lt+1
    write(70,*) abs(i-1), CF_ave(ivx(i),:), CF_delta(ivx(i),:)
  end do
  close(10)
  close(20)
  close(30)
  close(40)
  close(50)
  close(60)
  close(70)
  deallocate(phi_res,phi_err)
  deallocate(corr1,corr2,CF,CF_ave,CF_delta)
  end subroutine vary_m0
  
end program main
