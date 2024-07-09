! NameListEquilibrium.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
! Function to read EQUILIBRIUM namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

subroutine NameListEquilibrium (QC, NU, PC, MU, EPSA,&
     EPS, NS, NR, NF, NW,&
     ACC, H0, HMIN, HMAX)&
     bind (c, name = 'NameListEquilibrium')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  real    (kind = c_double), intent (inout) :: QC
  real    (kind = c_double), intent (inout) :: NU
  real    (kind = c_double), intent (inout) :: PC
  real    (kind = c_double), intent (inout) :: MU
  real    (kind = c_double), intent (inout) :: EPSA

  real    (kind = c_double), intent (inout) :: EPS
  integer (kind = c_int),    intent (inout) :: NS
  integer (kind = c_int),    intent (inout) :: NR
  integer (kind = c_int),    intent (inout) :: NF
  integer (kind = c_int),    intent (inout) :: NW

  real    (kind = c_double), intent (inout) :: ACC
  real    (kind = c_double), intent (inout) :: H0
  real    (kind = c_double), intent (inout) :: HMIN
  real    (kind = c_double), intent (inout) :: HMAX
 
  namelist /EQUILIBRIUM_CONTROL/ QC, NU, PC, MU, EPSA,&
       EPS, NS, NR, NF, NW,&
       ACC, H0, HMIN, HMAX
       
  open  (unit = 100, file = 'Inputs/Namelist.nml', status = 'old')
  read  (unit = 100, nml  = EQUILIBRIUM_CONTROL)
  close (unit = 100)

endsubroutine NameListEquilibrium
