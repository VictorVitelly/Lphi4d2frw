module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    integer(i4), parameter :: Lt=32,Lx=32,tau=10
    real(dp) :: lambda0=1._dp,at=real(tau,dp)/real(Lt,dp),ax=real(tau,dp)/real(Lt,dp)
    
    integer(i4), parameter :: thermalization=20000,Nmsrs=100,eachsweep=300,Nmsrs2=120
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    integer(i4), parameter :: Mbins=10,bins=101,Nauto=15000,Mbin(4)=(/5,10,15,20/)
    real(dp), parameter :: dphi_m=0.5_dp,hotphi=1._dp!, dphi=0.5_dp, hotphi=2._dp*dphi
    !real(dp) :: dphi=0.5_dp, hotphi=1._dp
    
    real(dp), parameter :: maxx=1.5_dp, minn=-1.5_dp
    real(dp) :: binwidth=(maxx-minn)/real(bins,dp)

    real :: starting,ending

end module parameters
