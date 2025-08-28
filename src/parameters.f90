module parameters

    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    integer(i4), parameter :: Lt=10,Lx=10
    real(dp), parameter :: PI=4._dp*DATAN(1._dp),tau=2*PI
    real(dp) :: lambda0=1._dp,at=tau/real(Lt,dp),ax=tau/real(Lt,dp)
    
    integer(i4), parameter :: thermalization=10000,Nmsrs=200,eachsweep=500,Nmsrs2=120
    integer(i4) :: sweeps=thermalization+eachsweep*Nmsrs
    integer(i4), parameter :: Mbins=10,bins=101,Nauto=15000,Mbin(4)=(/5,10,15,20/)
    real(dp), parameter :: dphi_m=0.5_dp,hotphi=1._dp!, dphi=0.5_dp, hotphi=2._dp*dphi
    !real(dp) :: dphi=0.5_dp, hotphi=1._dp
    
    real(dp), parameter :: maxx=1.5_dp, minn=-1.5_dp
    real(dp) :: binwidth=(maxx-minn)/real(bins,dp)

    real :: starting,ending

end module parameters
