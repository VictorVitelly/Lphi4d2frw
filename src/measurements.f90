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
    real(dp), dimension(Lt,Lx), intent(inout) :: corr2
    real(dp), dimension(Lt) :: varphi
    integer(i4) :: i1,i2
    varphi=0._dp
    do i1=1,Lt
      do i2=1,Lx
        varphi(i1)=varphi(i1)+phi(i1,i2)
      end do
    end do
    varphi(:)=varphi(:)/real(Lx,dp)
    do i1=1,Lt
      corr1(i1)=corr1(i1)+abs(varphi(i1))
      do i2=1,Lx
        corr2(i1,i2)=corr2(i1,i2)+phi(i1,i2)*phi(i1,1)
      end do
    end do
  end subroutine correlation
  
  subroutine secondmomentum(CF,xi2_ave,xi2_err)
    real(dp),dimension(Lt,Lx,Nmsrs2),intent(in) :: CF
    real(dp),intent(out),dimension(Lt) :: xi2_ave,xi2_err 
    integer(i4) :: i0,i1,i2
    real(dp),allocatable :: xi2(:,:),F1(:,:),F2(:,:),F12(:,:),F12_ave(:),F12_err(:)
    allocate(xi2(Lt,Nmsrs2),F1(Lt,Nmsrs2),F2(Lt,Nmsrs2),F12(Lt,Nmsrs2),F12_ave(Lt),F12_err(Lt))
    F1=0._dp
    F2=0._dp
    do i0=1,Lt
      do i1=1,Nmsrs2
        do i2=1,Lx
          F1(i0,i1)=F1(i0,i1)+CF(i0,i2,i1)
          F2(i0,i1)=F2(i0,i1)+CF(i0,i2,i1)*COS(real(i2-1,dp)*2._dp*PI/real(Lx,dp))
        end do
      end do
    end do
    do i0=1,Lt
      do i1=1,Nmsrs2
        F12(i0,i1)=F1(i0,i1)/F2(i0,i1)
      end do
    end do
    do i0=1,Lt
      call mean_scalar(F12(i0,:),F12_ave(i0),F12_err(i0))
      write(*,*) F12_ave(i0)
      xi2_ave(i0)=sqrt( (F12_ave(i0) -1._dp))/(2._dp*abs(SIN(PI/real(Lx,dp))) ) 
      xi2_err(i0)=F12_err(i0)/(4._dp*sqrt(F12_ave(i0)-1._dp)*abs(SIN(PI/real(Lx,dp))) )
    end do
    deallocate(xi2,F1,F2,F12,F12_ave,F12_err)
  end subroutine secondmomentum

  subroutine histogram(x,A1,A2)
  real(dp), dimension(:,:), intent(in) :: x
  integer(i4), dimension(bins,Lt), intent(inout) :: A1
  real(dp), dimension(bins), intent(in) :: A2
  integer(i4) :: i,j,k
  do i=1,bins
    do j=1,Lt
      do k=1,size(x,dim=2)
        if(x(j,k) .le. real(A2(i),dp)+binwidth/2._dp .and. x(j,k)>real(A2(i),dp)-binwidth/2._dp ) then
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
    real(dp), allocatable :: phi(:,:), A2(:), aux(:,:), aux_err(:,:)
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
    dphi=0.1_dp
    do i=1,thermalization
      !call montecarlo(m0,dphi,phi,AR)
      call metropolis(m0,dphi,phi)
    end do

    do i=1,Nmsrs2
      do j=1,Nmsrs
        call flip_sign(phi)
        do k=1,eachsweep
          call metropolis(m0,dphi,phi)
          !call montecarlo(m0,dphi,phi,AR)
        end do
        call montecarlo(m0,dphi,phi,AR)
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
      norm(:)=norm(:)+real(A1(i,:),dp)
    end do
    norm(:)=norm(:)*(real(maxx-minn,dp) )/real(bins,dp)
    allocate(aux(bins,Lt))
    allocate(aux_err(bins,Lt))
    do i=1,bins
      do j=1,Lt
      aux(i,j)=A1(i,j)/norm(j)
      aux_err(i,j)=sqrt( real(A1(i,j),dp) )/norm(j)
      end do
    end do
    
    do i=1,bins
      do j=1,Lt
        write(80,*) at*real(j,dp), A2(i), aux(i,j), aux_err(i,j)
      end do
    end do
    
    !do i=1,bins
    !  write(80,*) A2(i), aux(i,:), aux_err(i,:)
    !end do
    deallocate(aux,aux_err,A1,A2)
    close(80)
  end subroutine make_histogram

  subroutine fixed_m0(m0)
  real(dp),intent(in) :: m0
  real(dp) :: AR(Lt),dphi,ARp(Lt,Nmsrs2),AR_ave(Lt),AR_err(Lt),auxphi,auxsphi,tot,vol
  real(dp), allocatable :: phi(:,:),M2(:),M4(:),suscep(:,:),U4(:,:),S2(:),heat(:,:)
  real(dp), allocatable :: sus_ave(:),sus_err(:),U4_ave(:),U4_err(:),hea_ave(:),hea_err(:)
  real(dp), allocatable :: phi_xave(:),phi_ave(:,:),phi_res(:),phi_err(:)
  real(dp), allocatable :: Sphi_xave(:),Sphi_ave(:,:),Sphi_res(:),Sphi_err(:)
  real(dp), allocatable :: corr1(:),corr2(:,:),CF(:,:,:),CF_ave(:,:),CF_err(:,:)
  !real(dp), allocatable :: xi2_ave(:),xi2_err(:) 
  real(dp), allocatable :: norm(:), A2(:), aux(:,:), aux_err(:,:)
  integer(i4), allocatable :: A1(:,:)
  integer(i4) :: i,j,k1
  open(10, file = 'data/mag.dat', status = 'replace')
  open(20, file = 'data/ene.dat', status = 'replace')  
  open(30, file = 'data/cfs.dat', status = 'replace')
  open(40, file = 'data/sus.dat', status = 'replace')
  open(50, file = 'data/bin.dat', status = 'replace')
  !open(60, file = 'data/xi2.dat', status = 'replace')
  open(70, file = 'data/hgm.dat', status = 'replace')
  open(80, file = 'data/hea.dat', status = 'replace')
  allocate(phi(Lt+1,Lx))
  allocate(phi_res(Lt),phi_err(Lt),Sphi_res(Lt),Sphi_err(Lt))
  allocate(corr1(Lt),corr2(Lt,Lx),CF(Lt,Lx,Nmsrs2),CF_ave(Lt,Lx),CF_err(Lt,Lx))
  allocate(phi_xave(Lt),phi_ave(Nmsrs2,Lt),Sphi_xave(Lt),Sphi_ave(Nmsrs2,Lt))
  allocate(M2(Lt),M4(Lt),suscep(Nmsrs2,Lt),U4(Nmsrs2,Lt),S2(Lt),heat(Nmsrs2,Lt))
  allocate(sus_ave(Lt),sus_err(Lt),U4_ave(Lt),U4_err(Lt),hea_ave(Lt),hea_err(Lt))
  !allocate(xi2_ave(Lt),xi2_err(Lt))
  allocate(A1(bins,Lt),A2(bins),norm(Lt))
  dphi=0.1_dp
  ARp=0._dp
  CF=0._dp
  do i=1,bins
      A2(i)=minn+binwidth/2._dp+real(i-1,dp)*binwidth
  end do
  A1=0
  tot=real(Nmsrs,dp)
  vol=real(Lx,dp)
  !call hot_start(phi,hotphi)
  call cold_start(phi)
  !Thermalization
  do i=1,thermalization
    call metropolis(m0,dphi,phi)
  end do
    
  !Measurements
  do i=1,Nmsrs2
    phi_xave=0._dp
    Sphi_xave=0._dp
    M2=0._dp
    M4=0._dp
    S2=0._dp
    corr1(:)=0._dp
    corr2(:,:)=0._dp
    do j=1,Nmsrs
      call flip_sign(phi)
      do k1=1,eachsweep-1
        call metropolis(m0,dphi,phi)
      end do
      call montecarlo(m0,dphi,phi,AR)
      ARp(:,i)=ARp(:,i)+AR(:)
      do k1=1,Lt
        auxphi=meanphi_t(phi,k1)
        auxsphi=meanS_t(m0,phi,k1)
        phi_xave(k1)=phi_xave(k1)+abs(auxphi)
        Sphi_xave(k1)=Sphi_xave(k1)+auxsphi
        M2(k1)=M2(k1)+auxphi**2
        M4(k1)=M4(k1)+auxphi**4
        S2(k1)=S2(K1)+auxsphi**2
      end do
      call correlation(phi,corr1,corr2)
      call histogram(phi,A1,A2)
      !call flip_sign(phi,j)
    end do
    ARp(:,i)=ARp(:,i)/tot
    phi_ave(i,:)=phi_xave(:)/tot
    Sphi_ave(i,:)=Sphi_xave(:)/tot
    M2(:)=M2(:)/tot 
    M4(:)=M4(:)/tot
    S2(:)=S2(:)/tot
    corr1(:)=corr1(:)/tot
    corr2(:,:)=corr2(:,:)/tot
    do j=1,Lt
      do k1=1,Lx
        CF(j,k1,i)=corr2(j,k1) !-(phi_res(j)**2)
      end do
      suscep(i,j)=M2(j)-phi_ave(i,j)**2
      U4(i,j)=1._dp -M4(j)/(3._dp*M2(j)**2)
      heat(i,j)=S2(j)-Sphi_ave(i,j)**2
    end do
  end do
  
  
  !Averages
  do i=1,Lt
    call mean_scalar(ARp(i,:),AR_ave(i),AR_err(i))
    call mean_scalar(phi_ave(:,i),phi_res(i),phi_err(i))
    call mean_scalar(Sphi_ave(:,i),Sphi_res(i),Sphi_err(i))
    do j=1,Lx
      call mean_scalar(CF(i,j,:),CF_ave(i,j),CF_err(i,j))
    end do
    call mean_scalar(suscep(:,i),sus_ave(i),sus_err(i))
    call mean_scalar(U4(:,i),U4_ave(i),U4_err(i))
    call mean_scalar(heat(:,i),hea_ave(i),hea_err(i))
  end do
  phi_res(:)=phi_res(:)/vol
  phi_err(:)=phi_err(:)/vol
  Sphi_res(:)=Sphi_res(:)/vol
  Sphi_err(:)=Sphi_err(:)/vol
  sus_ave(:)=sus_ave(:)/vol
  sus_err(:)=sus_err(:)/vol
  hea_ave(:)=hea_ave(:)/vol
  hea_err(:)=hea_err(:)/vol
  
  !Histograms
  norm(:)=0._dp
  do i=1,bins
    norm(:)=norm(:)+real(A1(i,:),dp)
  end do
  norm(:)=norm(:)*(real(maxx-minn,dp) )/real(bins,dp)
  allocate(aux(bins,Lt))
  allocate(aux_err(bins,Lt))
  do i=1,bins
    do j=1,Lt
    aux(i,j)=A1(i,j)/norm(j)
    aux_err(i,j)=sqrt( real(A1(i,j),dp) )/norm(j)
    end do
  end do
  
  !Second momentum correlation length
  !do i=1,Nmsrs2
  !  do j=1,Lt 
  !    do k1=1,Lx
  !      CF(j,k1,i)=corr2(j,k1) -(phi_res(j)**2)
  !    end do
  !  end do
  !end do
  !call secondmomentum(CF,xi2_ave,xi2_err)
  
    
  write(*,*) 'The acc. rate for m02=',m0
  do i=1,Lt
    write(*,*) i, AR_ave(i), AR_err(i)
  end do

  
  do i=1,Lt
    write(10,*) etac(i), phi_res(i), phi_err(i)
    write(20,*) etac(i), Sphi_res(i), Sphi_err(i)
    write(40,*) etac(i), sus_ave(i), sus_err(i)
    write(50,*) etac(i), U4_ave(i), U4_err(i)
    write(80,*) etac(i), hea_ave(i), hea_err(i)
    !write(60,*) etac(i), xi2_ave(i), xi2_err(i)
  end do
  do i=1,Lx+1
    write(30,*) abs(i-1), CF_ave(:,ivx(i)), CF_err(:,ivx(i))
  end do
  do i=1,bins
    do j=1,Lt
      write(70,*) etac(j), A2(i), aux(i,j), aux_err(i,j)
    end do
  end do
  close(10)
  close(20)
  close(30)
  close(40)
  close(50)
  !close(60)
  close(70)
  close(80)
  deallocate(phi)
  deallocate(phi_xave,phi_ave,Sphi_xave,Sphi_ave)
  deallocate(phi_res,phi_err)
  deallocate(M2,M4,suscep,U4,sus_ave,sus_err,U4_ave,U4_err)
  deallocate(corr1,corr2,CF,CF_ave,CF_err,hea_ave,hea_err)
  deallocate(aux,aux_err,A1,A2)
  end subroutine fixed_m0


end module measurements
