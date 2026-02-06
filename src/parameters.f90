module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    integer(i4), parameter :: Lt=16,Lx=16,Ttot=16,Xtot=16
    real(dp), parameter :: PI=4._dp*DATAN(1._dp)
    real(dp) :: lambda0=0.08_dp,H0=0.1_dp
    real(dp) :: at=real(Lt,dp)/real(Ttot,dp),ax=real(Lx,dp)/real(Xtot,dp)
    
    integer(i4), parameter :: thermalization=5000,Nmsrs=100,eachsweep=300,Nmsrs2=120
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    integer(i4), parameter :: Mbins=10,bins=101,Nauto=15000,Mbin(4)=(/5,10,15,20/)
    real(dp), parameter :: dphi_m=0.5_dp,hotphi=1._dp!, dphi=0.5_dp, hotphi=2._dp*dphi
    !real(dp) :: dphi=0.5_dp, hotphi=1._dp
    
    real(dp), parameter :: maxx=1.5_dp, minn=-1.5_dp
    real(dp) :: binwidth=(maxx-minn)/real(bins,dp)

    real :: starting,ending

end module parameters
