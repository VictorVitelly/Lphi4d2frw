module functions
    use iso_fortran_env, only : dp => real64, i4 => int32
    use parameters
    implicit none

contains

  function ivx(i)
    integer(i4), intent(in) :: i
    integer(i4) :: ivx
    if(i==Lx+1) then
      ivx=1
    else if(i==0) then
      ivx=Lx
    else
      ivx=i
    end if
  end function ivx

  function ivt(i)
    integer(i4), intent(in) :: i
    integer(i4) :: ivt
    if(i==0) then
      ivt=Lt+1
    else
      ivt=i
    end if
  end function ivt
  
  function alfa(t)
    integer(i4),intent(in) :: t
    real(dp) :: alfa
    !alfa=1._dp
    alfa=exp(H0*at*real(t-1,dp))
    !alfa=sin(at*real(t-1,dp))
  end function alfa

  function lagrangian(m02,phi,i1,i2)
    real(dp), intent(in) :: m02
    real(dp), dimension(:,:), intent(in) :: phi
    integer(i4), intent(in) :: i1,i2
    real(dp) :: lagrangian
    real(dp) :: laga,lagb,lagc
    laga=(phi(i1+1,i2)-phi(i1,i2) )**2 /at**2
    lagb=(phi(i1,ivx(i2+1))-phi(i1,i2) )**2 /ax**2
    lagc=alfa(i1)**2 *(m02*phi(i1,i2)**2+0.5_dp*lambda0*phi(i1,i2)**4)
    lagrangian=(laga+lagb+lagc)/2._dp
  end function lagrangian

  function S(m02,phi)
    real(dp), intent(in) :: m02
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    real(dp) :: S
    integer(i4) :: i1,i2
    S=0._dp
    do i1=1,Lt-1
      do i2=1,Lx
        S=S+lagrangian(m02,phi,i1,i2)
      end do
    end do
    S=at*ax*S
  end function S

  function DeltaSdbc(m02,phi,i1,i2,phi2)
    real(dp), intent(in) :: m02
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    integer(i4), intent(in) :: i1,i2
    real(dp), intent(in) :: phi2
    real(dp) :: DeltaSdbc
    real(dp) :: DSa,DSb,DSc,DSd
      DSa=(phi2**2-phi(i1,i2)**2)*(1._dp/at**2+1._dp/ax**2)
      DSb=-(phi2-phi(i1,i2))*(phi(ivt(i1+1),i2) +phi(ivt(i1-1),i2))/at**2 
      DSc=-(phi2-phi(i1,i2))*(phi(i1,ivx(i2+1)) +phi(i1,ivx(i2-1)))/ax**2
      DSd=alfa(i1)**2 *(m02*(phi2**2-phi(i1,i2)**2)&
          &+lambda0*(phi2**4-phi(i1,i2)**4)/2._dp)/2._dp
      DeltaSdbc=at*ax*(DSa+DSb+DSc+DSd)
  end function DeltaSdbc
  
  function DeltaSgbc(m02,phi,i1,i2,phi2)
    real(dp), intent(in) :: m02
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    integer(i4), intent(in) :: i1,i2
    real(dp), intent(in) :: phi2
    real(dp) :: DeltaSgbc
    real(dp) :: DSa,DSb,DSc,DSd
    If(i1==Lt) then 
      DSa=(phi2**2-phi(i1,i2)**2)*(1.0_dp/at**2+1._dp/ax**2)
      DSb=-2._dp*(phi2-phi(i1,i2))*(phi(i1-1,i2))/at**2 
      DSc=-(phi2-phi(i1,i2))*(phi(i1,ivx(i2+1)) +phi(i1,ivx(i2-1)))/ax**2
      DSd=alfa(i1)**2 *(m02*(phi2**2-phi(i1,i2)**2)&
          &+lambda0*(phi2**4-phi(i1,i2)**4)/2._dp)/2._dp
      DeltaSgbc=at*ax*(DSa+DSb+DSc+DSd)
    else if (i1==1) then 
      DSa=(phi2**2-phi(i1,i2)**2)*(1.0_dp/at**2+1._dp/ax**2)
      DSb=-2._dp*(phi2-phi(i1,i2))*(phi(i1+1,i2) )/at**2 
      DSc=-(phi2-phi(i1,i2))*(phi(i1,ivx(i2+1)) +phi(i1,ivx(i2-1)))/ax**2
      DSd=alfa(i1)**2 *(m02*(phi2**2-phi(i1,i2)**2)&
          &+lambda0*(phi2**4-phi(i1,i2)**4)/2._dp)/2._dp
      DeltaSgbc=at*ax*(DSa+DSb+DSc+DSd)
    else 
      DSa=(phi2**2-phi(i1,i2)**2)*(1._dp/at**2+1._dp/ax**2)
      DSb=-(phi2-phi(i1,i2))*(phi(i1+1,i2) +phi(i1-1,i2))/at**2 
      DSc=-(phi2-phi(i1,i2))*(phi(i1,ivx(i2+1)) +phi(i1,ivx(i2-1)))/ax**2
      DSd=alfa(i1)**2 *(m02*(phi2**2-phi(i1,i2)**2)&
          &+lambda0*(phi2**4-phi(i1,i2)**4)/2._dp)/2._dp
      DeltaSgbc=at*ax*(DSa+DSb+DSc+DSd)
    end if
  end function DeltaSgbc
  
  function DeltaSfbc(m02,phi,i1,i2,phi2)
    real(dp), intent(in) :: m02
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    integer(i4), intent(in) :: i1,i2
    real(dp), intent(in) :: phi2
    real(dp) :: DeltaSfbc
    real(dp) :: DSa,DSb,DSc,DSd
    If(i1==Lt) then 
      DSa=(phi2**2-phi(i1,i2)**2)*(0.5_dp/at**2+1._dp/ax**2)
      DSb=-(phi2-phi(i1,i2))*(phi(i1-1,i2))/at**2 
      DSc=-(phi2-phi(i1,i2))*(phi(i1,ivx(i2+1)) +phi(i1,ivx(i2-1)))/ax**2
      DSd=alfa(i1)**2 *(m02*(phi2**2-phi(i1,i2)**2)&
          &+lambda0*(phi2**4-phi(i1,i2)**4)/2._dp)/2._dp
      DeltaSfbc=at*ax*(DSa+DSb+DSc+DSd)
    else if (i1==1) then 
      DSa=(phi2**2-phi(i1,i2)**2)*(0.5_dp/at**2+1._dp/ax**2)
      DSb=-(phi2-phi(i1,i2))*(phi(i1+1,i2) )/at**2 
      DSc=-(phi2-phi(i1,i2))*(phi(i1,ivx(i2+1)) +phi(i1,ivx(i2-1)))/ax**2
      DSd=alfa(i1)**2 *(m02*(phi2**2-phi(i1,i2)**2)&
          &+lambda0*(phi2**4-phi(i1,i2)**4)/2._dp)/2._dp
      DeltaSfbc=at*ax*(DSa+DSb+DSc+DSd)
    else 
      DSa=(phi2**2-phi(i1,i2)**2)*(1._dp/at**2+1._dp/ax**2)
      DSb=-(phi2-phi(i1,i2))*(phi(i1+1,i2) +phi(i1-1,i2))/at**2 
      DSc=-(phi2-phi(i1,i2))*(phi(i1,ivx(i2+1)) +phi(i1,ivx(i2-1)))/ax**2
      DSd=alfa(i1)**2 *(m02*(phi2**2-phi(i1,i2)**2)&
          &+lambda0*(phi2**4-phi(i1,i2)**4)/2._dp)/2._dp
      DeltaSfbc=at*ax*(DSa+DSb+DSc+DSd)
    end if
  end function DeltaSfbc
  
  function DeltaSrbc(m02,phi,i1,i2,phi2)
    real(dp), intent(in) :: m02
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    integer(i4), intent(in) :: i1,i2
    real(dp), intent(in) :: phi2
    real(dp) :: DeltaSrbc
    real(dp) :: DSa,DSb,DSc,DSd
    If(i1==Lt) then 
      DSa=(phi(Lt,i2)-phi2)**2 /(2._dp*at**2)
      DSb=(phi2-phi(Lt-1,i2))**2 /(2._dp*at**2)
      DSc=-(phi(Lt,i2)-phi(Lt-1,i2))**2 /(2._dp*at**2)
      DSd=alfa(i1)**2 *(m02*(phi2**2-phi(i1,i2)**2)&
          &+lambda0*(phi2**4-phi(i1,i2)**4)/2._dp)/2._dp
      DeltaSrbc=at*ax*(DSa+DSb+DSc+DSd)
    else if (i1==1) then 
      DSa=(phi(2,i2)-phi2)**2 /(2._dp*at**2)
      DSb=(phi2-phi(1,i2))**2 /(2._dp*at**2)
      DSc=-(phi(2,i2)-phi(1,i2))**2 /(2._dp*at**2)
      DSd=alfa(i1)**2 *(m02*(phi2**2-phi(i1,i2)**2)&
          &+lambda0*(phi2**4-phi(i1,i2)**4)/2._dp)/2._dp
      DeltaSrbc=at*ax*(DSa+DSb+DSc+DSd)
    else 
      DSa=(phi2**2-phi(i1,i2)**2)*(1._dp/at**2+1._dp/ax**2)
      DSb=-(phi2-phi(i1,i2))*(phi(i1+1,i2) +phi(i1-1,i2))/at**2 
      DSc=-(phi2-phi(i1,i2))*(phi(i1,ivx(i2+1)) +phi(i1,ivx(i2-1)))/ax**2
      DSd=alfa(i1)**2 *(m02*(phi2**2-phi(i1,i2)**2)&
          &+lambda0*(phi2**4-phi(i1,i2)**4)/2._dp)/2._dp
      DeltaSrbc=at*ax*(DSa+DSb+DSc+DSd)
    end if
  end function DeltaSrbc

  function DeltaSpbc(m02,phi,i1,i2,phi2)
    real(dp), intent(in) :: m02
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    integer(i4), intent(in) :: i1,i2
    real(dp), intent(in) :: phi2
    real(dp) :: DeltaSpbc
    real(dp) :: DSa,DSb,DSc,DSd
      DSa=(phi2**2-phi(i1,i2)**2)*(1._dp/at**2 +1._dp/ax**2)
      DSb=-(phi2-phi(i1,i2))*(phi(ivx(i1+1),i2) +phi(ivx(i1-1),i2))/at**2
      DSc=-(phi2-phi(i1,i2))*(phi(i1,ivx(i2+1)) +phi(i1,ivx(i2-1)))/ax**2
      DSd=alfa(i1)**2 *( m02*(phi2**2-phi(i1,i2)**2) &
          &+lambda0*(phi2**4-phi(i1,i2)**4)/2._dp)/2._dp
      DeltaSpbc=at*ax*(DSa+DSb+DSc+DSd)
  end function DeltaSpbc

  function meanphi(phi)
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    integer(i4):: i1,i2
    real(dp) :: meanphi
    meanphi=0._dp
    do i1=1,Lt
      do i2=1,Lx
        meanphi=meanphi+phi(i1,i2)
      end do
    end do
  end function meanphi
  
  function meanphi_t(phi,t)
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    integer(i4), intent(in) :: t
    real(dp) :: meanphi_t
    integer(i4) :: i
    meanphi_t=0._dp
    do i=1,Lx
      meanphi_t=meanphi_t+phi(t,i)
    end do
  end function meanphi_t
  
  function meandphi_t(phi,t)
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    integer(i4), intent(in) :: t
    real(dp) :: meandphi_t
    integer(i4) :: i
    meandphi_t=0._dp
    do i=1,Lx
      meandphi_t=meandphi_t+phi(t+1,i)-phi(t,i)
    end do
  end function meandphi_t
  
  function meanS_t(m0,phi,t)
    real(dp), dimension(Lt+1,Lx), intent(in) :: phi
    real(dp), intent(in) :: m0
    integer(i4), intent(in) :: t
    real(dp) :: meanS_t
    integer(i4) :: i
    meanS_t=0._dp
    do i=1,Lx
      meanS_t=meanS_t+lagrangian(m0,phi,t,i)
    end do
  end function meanS_t
  
end module functions
