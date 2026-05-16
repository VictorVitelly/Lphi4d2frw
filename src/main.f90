program main

  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use arrays
  use functions
  use statistics
  use measurements
  implicit none

  call init_vecs()
  call cpu_time(starting)

  !Thermalization history and autocorretion functions
  !call thermalize(-1._dp)
  !call test(-1.0_dp)

  !Measure action, magnetization, susceptibility and heat cap.
  !call vary_m0(0._dp,-3.0_dp,11)
  !call vary_m0(0.0_dp,-3.0_dp,11)
  !call vary_m0(0._dp,1.5_dp,11)

  call fixed_m0(-1._dp)
  !call make_histogram(-1.0_dp)

  call cpu_time(ending)
  write(*,*) "Elapsed time: ", (ending-starting), " s"
  
end program main
