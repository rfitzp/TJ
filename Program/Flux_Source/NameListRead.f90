! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read FLUX namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (NPSI, PACK, NTHETA, H0, ACC)&
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: NPSI
  real    (kind = c_double), intent (inout) :: PACK
  integer (kind = c_int),    intent (inout) :: NTHETA
  real    (kind = c_double), intent (inout) :: H0
  real    (kind = c_double), intent (inout) :: ACC

  namelist /FLUX_CONTROL/ NPSI, PACK, NTHETA, H0, ACC
  open  (unit = 100, file = 'Inputs/Flux.nml', status = 'old')
  read  (unit = 100, nml  = FLUX_CONTROL)
  close (unit = 100)

endsubroutine namelistRead
