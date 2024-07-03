! NameListTJ.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read TJ namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListTJ (NTOR, MMIN, MMAX,&
     EPS, DEL, NFIX, NDIAG, NULC, ITERMAX, FREE, TVAC,&
     ACC, H0, HMIN, HMAX, EPSF, POWR) &
     bind (c, name = 'NameListTJ')
  
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  integer (kind = c_int),    intent (inout) :: NTOR
  integer (kind = c_int),    intent (inout) :: MMIN
  integer (kind = c_int),    intent (inout) :: MMAX

  real    (kind = c_double), intent (inout) :: EPS
  real    (kind = c_double), intent (inout) :: DEL
  integer (kind = c_int),    intent (inout) :: NFIX
  integer (kind = c_int),    intent (inout) :: NDIAG
  real    (kind = c_double), intent (inout) :: NULC
  integer (kind = c_int),    intent (inout) :: ITERMAX
  integer (kind = c_int),    intent (inout) :: FREE
  integer (kind = c_int),    intent (inout) :: TVAC

  real    (kind = c_double), intent (inout) :: ACC
  real    (kind = c_double), intent (inout) :: H0
  real    (kind = c_double), intent (inout) :: HMIN
  real    (kind = c_double), intent (inout) :: HMAX
  real    (kind = c_double), intent (inout) :: EPSF
  real    (kind = c_double), intent (inout) :: POWR
 
  namelist /TJ_CONTROL/ NTOR, MMIN, MMAX,&
       EPS, DEL, NFIX, NDIAG, NULC, ITERMAX, FREE, TVAC,&
       ACC, H0, HMIN, HMAX, EPSF, POWR
       
  open  (unit = 100, file = 'Inputs/Namelist.nml', status = 'old')
  read  (unit = 100, nml  = TJ_CONTROL)
  close (unit = 100)

endsubroutine NameListTJ
