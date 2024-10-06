! WriteEFIT.f

! Program to input data from Outputs/Equilibrium/EFIT.nc and write EFIT file Outputs/WriteEFIT/EFIT.txt      

program Write_EFIT

use netcdf
implicit none

integer :: ncid, err, varid, file_unit, i, j, k
integer :: ivals(4), values(8), i3
integer :: NRBOX, NZBOX, NPBOUND, NLIMITER

double precision :: rvals(15), zero
double precision :: RBOXLEN, ZBOXLEN, RBOXLFT, ZOFF, R0EXP, B0EXP
double precision :: RAXIS, ZAXIS, PSIAXIS, PSIBOUND, CURRENT
double precision, dimension (:),    allocatable :: T, P, TTP, PP, Q, RBOUND, ZBOUND, RLIMITER, ZLIMITER
double precision, dimension (:, :), allocatable :: PSI

character(len=*), parameter :: filename = '../Outputs/Equilibrium/EFIT.nc'
character(len=*), parameter :: comment  = ' Output from TJ code: '
character(len=8)  :: date_str
character(len=10) :: time_str
character(len=48) :: string
character(len=60) :: efit_name

i3   = 3
zero = 0.

! -----------------------------
! Get the current date and time
! -----------------------------
call date_and_time (date_str, time_str, values=values)
string = comment // date_str

! ----------------
! Open netcdf file
! ----------------
err = nf90_open (filename, NF90_NOWRITE, ncid)
if (err /= nf90_noerr) then
   print *, 'Error opening file: ', trim (filename)
   stop
end if

! -----------------------
! Read integer parameters
! -----------------------
err = nf90_inq_varid(ncid, 'IntegerParameters', varid)
if (err /= nf90_noerr) then
   print *, 'Error getting variable ID for IntegerParameters'
   err =  nf90_close(ncid)
   stop
end if

err = nf90_get_var(ncid, varid, ivals)
if (err /= nf90_noerr) then
   print *, 'Error reading IntegerParameters data'
   err =  nf90_close(ncid)
   stop
end if

NRBOX    = ivals(1)
NZBOX    = ivals(2)
NPBOUND  = ivals(3)
NLIMITER = ivals(4)

allocate (T        (NRBOX))
allocate (P        (NRBOX))
allocate (TTp      (NRBOX))
allocate (Pp       (NRBOX))
allocate (PSI      (NRBOX, NZBOX))
allocate (Q        (NRBOX))
allocate (RBOUND   (NPBOUND))
allocate (ZBOUND   (NPBOUND))
allocate (RLIMITER (NLIMITER))
allocate (ZLIMITER (NLIMITER))

! -------------------------
! Read equilibrium profiles
! -------------------------
err = nf90_inq_varid (ncid, 'T', varid)
if (err /= nf90_noerr) then
   print *, 'Error getting variable ID for T'
   err = nf90_close (ncid)
   stop
end if

err = nf90_get_var (ncid, varid, T)
if (err /= nf90_noerr) then
   print *, 'Error reading T data'
   err = nf90_close (ncid)
   stop
end if

err = nf90_inq_varid (ncid, 'P', varid)
if (err /= nf90_noerr) then
   print *, 'Error getting variable ID for P'
   err = nf90_close (ncid)
   stop
end if

err = nf90_get_var (ncid, varid, P)
if (err /= nf90_noerr) then
   print *, 'Error reading P data'
   err = nf90_close (ncid)
   stop
end if

err = nf90_inq_varid (ncid, 'TTp', varid)
if (err /= nf90_noerr) then
   print *, 'Error getting variable ID for TTP'
   err = nf90_close (ncid)
   stop
end if

err = nf90_get_var (ncid, varid, TTP)
if (err /= nf90_noerr) then
   print *, 'Error reading TTP data'
   err = nf90_close (ncid)
   stop
end if

err = nf90_inq_varid (ncid, 'Pp', varid)
if (err /= nf90_noerr) then
   print *, 'Error getting variable ID for PP'
   err =  nf90_close (ncid)
   stop
end if

err = nf90_get_var (ncid, varid, PP)
if (err /= nf90_noerr) then
   print *, 'Error reading PP data'
   err = nf90_close (ncid)
   stop
end if

open (unit = 101, file = '../Outputs/Equilibrium/PsiSequential.txt')
do i = 1, NRBOX
   do j = 1, NZBOX
      read (101, '(i4,1x,i4,1x,e17.9)') k, k, PSI (i, j)
   end do
end do
close (unit = 101)

err = nf90_inq_varid (ncid, 'Q', varid)
if (err /= nf90_noerr) then
   print *, 'Error getting variable ID for Q'
   err = nf90_close (ncid)
   stop
end if

err = nf90_get_var (ncid, varid, Q)
if (err /= nf90_noerr) then
   print *, 'Error reading Q data'
   err = nf90_close (ncid)
   stop
end if

err = nf90_inq_varid (ncid, 'RBOUND', varid)
if (err /= nf90_noerr) then
   print *, 'Error getting variable ID for RBOUND'
   err = nf90_close (ncid)
   stop
end if

err = nf90_get_var (ncid, varid, RBOUND)
if (err /= nf90_noerr) then
   print *, 'Error reading RBOUND data'
   err = nf90_close (ncid)
   stop
end if

err = nf90_inq_varid (ncid, 'ZBOUND', varid)
if (err /= nf90_noerr) then
   print *, 'Error getting variable ID for ZBOUND'
   err = nf90_close (ncid)
   stop
end if

err = nf90_get_var (ncid, varid, ZBOUND)
if (err /= nf90_noerr) then
   print *, 'Error reading ZBOUND data'
   err = nf90_close (ncid)
   stop
end if

err = nf90_inq_varid (ncid, 'RLIMITER', varid)
if (err /= nf90_noerr) then
   print *, 'Error getting variable ID for RLIMITER'
   err = nf90_close (ncid)
   stop
end if

err = nf90_get_var (ncid, varid, RLIMITER)
if (err /= nf90_noerr) then
   print *, 'Error reading RLIMITER data'
   err =  nf90_close (ncid)
   stop
end if

err = nf90_inq_varid (ncid, 'ZLIMITER', varid)
if (err /= nf90_noerr) then
   print *, 'Error getting variable ID for ZLIMITER'
   err =  nf90_close (ncid)
   stop
end if

err = nf90_get_var (ncid, varid, ZLIMITER)
if (err /= nf90_noerr) then
   print *, 'Error reading ZLIMITER data'
   err =  nf90_close (ncid)
   stop
end if

! --------------------
! Read real parameters
! --------------------
err = nf90_inq_varid(ncid, 'RealParameters', varid)
if (err /= nf90_noerr) then
   print *, 'Error getting variable ID for RealParameters'
   err =  nf90_close(ncid)
   stop
end if

err = nf90_get_var(ncid, varid, rvals)
if (err /= nf90_noerr) then
   print *, 'Error reading RealParameters data'
   err =  nf90_close(ncid)
   stop
end if

RBOXLEN  = rvals(1)
ZBOXLEN  = rvals(2)
RBOXLFT  = rvals(3)
ZOFF     = rvals(4)
R0EXP    = rvals(5)
B0EXP    = rvals(6)
RAXIS    = rvals(7)
ZAXIS    = rvals(8)
PSIAXIS  = rvals(9)
PSIBOUND = rvals(10)
CURRENT  = rvals(11)

! -----------------
! Close netcdf file
! -----------------
err = nf90_close(ncid)
if (err /= nf90_noerr) then
   print *, 'Error closing file'
   stop
end if

! ---------------
! Write EFIT file
! ---------------
efit_name = "../Outputs/WriteEFIT/EFIT.txt"
  
open (unit = 10, file = efit_name, status = "replace", action = "write", iostat = file_unit)
    
! Check for errors in opening the file
if (file_unit /= 0) then
   print *, "Error opening Outputs/WriteEFIT/EFIT.txt"
   stop
end if

write (10, '(a48, 3i4)') string,  i3,      NRBOX,   NZBOX
write (10, '(5e16.9)'  ) RBOXLEN, ZBOXLEN, R0EXP,   RBOXLFT,  ZOFF
write (10, '(5e16.9)'  ) RAXIS,   ZAXIS,   PSIAXIS, PSIBOUND, B0EXP
write (10, '(5e16.9)'  ) CURRENT, PSIAXIS, zero,    RAXIS,    zero
write (10, '(5e16.9)'  ) ZAXIS,   zero,    zero,    zero,     zero

write (10, '(5e16.9)') (T   (i),    i = 1, NRBOX)
write (10, '(5e16.9)') (P   (i),    i = 1, NRBOX)
write (10, '(5e16.9)') (TTp (i),    i = 1, NRBOX)
write (10, '(5e16.9)') (PP  (i),    i = 1, NRBOX)
write (10, '(5e16.9)') ((PSI(i, j), i = 1, NRBOX), j = 1, NZBOX)
write (10, '(5e16.9)') (Q   (i),    i = 1, NRBOX)

write (10, '(2i5)')     NPBOUND, NLIMITER
write (10, '(5e16.9)') (RBOUND  (i), ZBOUND  (i), i = 1, NPBOUND)
write (10, '(5e16.9)') (RLIMITER(i), ZLIMITER(i), i = 1, NLIMITER)

close (10)

! --------
! Clean up
! --------
deallocate (T)
deallocate (P)
deallocate (TTp)
deallocate (Pp)
deallocate (PSI)
deallocate (Q)
deallocate (RBOUND)
deallocate (ZBOUND)
deallocate (RLIMITER)
deallocate (ZLIMITER)

stop
      
end program Write_EFIT

 
