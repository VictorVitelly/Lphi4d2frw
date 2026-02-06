module statistics
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  implicit none

contains

  subroutine random_phi(x,bound)
    implicit none
    real(dp),intent(out) :: x
    real(dp), intent(in) :: bound
    real(dp) :: y
    call random_number(y)
    x = 2._dp*bound*y -bound
  end subroutine random_phi

  subroutine cold_start(phi)
    real(dp), dimension(:,:), intent(out) :: phi
    phi=0.0_dp
  end subroutine cold_start

  subroutine hot_start(phi,hotphi)
    real(dp), dimension(Lt+1,Lx), intent(out) :: phi
    real(dp), intent(in) :: hotphi
    integer(i4) :: i1,i2
    do i1=1,Lt+1
      do i2=1,Lx
        call random_phi(phi(i1,i2),hotphi)
      end do
    end do
  end subroutine hot_start

  subroutine metropolis(m0,dphi,phi)
    real(dp), intent(in) :: m0,dphi
    real(dp), dimension(Lt+1,Lx), intent(inout) :: phi
    real(dp) :: deltaphi,phi2,DS,r,p
    integer(i4) :: i1,i2
    do i1=1,Lt
      do i2=1,Lx
        call random_phi(deltaphi,dphi)
        phi2=phi(i1,i2)+deltaphi
        !DS=DeltaSdbc(m0,phi,i1,i2,phi2)
        !DS=DeltaSfbc(m0,phi,i1,i2,phi2)
        DS=DeltaSgbc(m0,phi,i1,i2,phi2)
        !DS=DeltaSrbc(m0,phi,i1,i2,phi2)
        !DS=DeltaSpbc(m0,phi,i1,i2,phi2)
        if(DS .le. 0._dp) then
          phi(i1,i2)=phi2
        else
          call random_number(r)
          p=Exp(-DS)
          if(r < p ) then
            phi(i1,i2)=phi2
          end if
        end if
      end do
    end do
  end subroutine metropolis
  
  subroutine montecarlo(m0,dphi,phi,AR)
    real(dp), intent(in) :: m0,dphi
    real(dp), dimension(Lt+1,Lx), intent(inout) :: phi
    real(dp), dimension(Lt), intent(out) :: AR
    real(dp) :: deltaphi,phi2,DS,r,p
    integer(i4) :: i1,i2
    AR=0._dp
    do i1=1,Lt
      do i2=1,Lx
        call random_phi(deltaphi,dphi-0.08_dp*real(i1,dp)/real(Lt,dp))
        phi2=phi(i1,i2)+deltaphi
        !DS=DeltaSdbc(m0,phi,i1,i2,phi2)
        !DS=DeltaSfbc(m0,phi,i1,i2,phi2)
        DS=DeltaSgbc(m0,phi,i1,i2,phi2)
        !DS=DeltaSrbc(m0,phi,i1,i2,phi2)
        !DS=DeltaSpbc(m0,phi,i1,i2,phi2)
        if(DS .le. 0._dp) then
          phi(i1,i2)=phi2
          AR(i1)=AR(i1)+1._dp
        else
          call random_number(r)
          p=Exp(-DS)
          AR(i1)=AR(i1)+p
          if(r < p ) then
            phi(i1,i2)=phi2
          end if
        end if
      end do
    end do
    AR=AR/real(Lx,dp)
  end subroutine montecarlo

  subroutine montecarlopbc(m0,dphi,phi,AR)
    real(dp), intent(in) :: m0,dphi
    real(dp), dimension(Lt+1,Lx), intent(inout) :: phi
    real(dp),intent(out) :: AR
    real(dp) :: deltaphi,phi2,DS,r,p
    integer(i4) :: i1,i2
    AR=0._dp
    do i1=1,Lt
      do i2=1,Lx
        call random_phi(deltaphi,dphi)
        phi2=phi(i1,i2)+deltaphi
        DS=DeltaSpbc(m0,phi,i1,i2,phi2)
        if(DS .le. 0._dp) then
          phi(i1,i2)=phi2
          AR=AR+1._dp
        else
          call random_number(r)
          p=Exp(-DS)
          AR=AR+p
          if(r < p ) then
            phi(i1,i2)=phi2
          end if
        end if
      end do
    end do
    AR=AR/real(Lx*Lt,dp)
  end subroutine montecarlopbc
  
  subroutine flip_sign(phi)
    real(dp), dimension(Lt+1,Lx), intent(inout) :: phi
    integer(i4) :: i1,i2
      do i1=1,Lt
        do i2=1,Lx
          phi(i1,i2)=-phi(i1,i2)
        end do
      end do
  end subroutine flip_sign
  
  !Statistics for measurements
  !
  !
  
  subroutine jackknife(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: jackk
    real(dp), allocatable :: xmean(:), delta_y(:)
    integer(i4) :: k,Narr,i,j
      Narr=size(x)
      allocate(delta_y(size(Mbin)))
      do j=1,size(Mbin)
        allocate(xmean(Mbin(j)))
        jackk=0._dp
        xmean=0._dp
        do i=1,Mbin(j)
          do k=1,Narr
            if(k .le. (i-1)*Narr/Mbin(j)) then
              xmean(i)=xmean(i)+x(k)
            else if(k > i*Narr/Mbin(j)) then
              xmean(i)=xmean(i)+x(k)
            end if
          end do
          xmean(i)=xmean(i)/(real(Narr,dp) -real(Narr/Mbin(j),dp))
        end do
        do k=1,Mbin(j)
          jackk=jackk+(xmean(k)-y )**2
        end do
        delta_y(j)=Sqrt(real(Mbin(j)-1,dp)*jackk/real(Mbin(j),dp))
        deallocate(xmean)
      end do
      deltay=maxval(delta_y)
  end subroutine jackknife

  subroutine standard_error(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: variance
    integer(i4) :: k,Narr
    Narr=size(x)
    deltay=0._dp
    variance=0._dp
    do k=1,Narr
      variance=variance+(x(k) -y)**2
    end do
    variance=variance/real(Narr-1,dp)
    deltay=Sqrt(variance/real(Narr,dp))
  end subroutine standard_error

  subroutine mean_0(x,y)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y
    integer(i4) :: k,Narr
    Narr=size(x)
    y=0._dp
    do k=1,Narr
      y=y+x(k)
    end do
    y=y/real(Narr,dp)
  end subroutine mean_0

  subroutine mean_scalar(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y,deltay
    call mean_0(x,y)
    !call standard_error(x,y,deltay)
    call jackknife(x,y,deltay)
  end subroutine mean_scalar

end module statistics
