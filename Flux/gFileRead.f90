! gFileRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read EFIT file and output data to Stage1 files
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine gFileRead () bind (c, name = 'gFileRead')

  use netcdf
  implicit none
  
  character (len = 100) :: string
  integer               :: i,       j,        i3,      NRBOX,   NZBOX, NBOUND, NLIM
  double precision      :: RBOXLEN, ZBOXLEN,  R0,      RBOXLFT, zero,  ZOFF,   RMAX
  double precision      :: RAXIS,   ZAXIS,    B0,      MU0,     ZMAX
  double precision      :: RLEFT,   RRIGHT,   ZLOW,    ZHIGH
  double precision      :: PSIAXIS, PSIBOUND, CURRENT

  double precision, dimension (:),    allocatable :: PSIN, T,     P,       TTp,    Pp,   Q
  double precision, dimension (:),    allocatable :: R,    Z,     RBOUND,  ZBOUND, RLIM, ZLIM  
  double precision, dimension (:, :), allocatable :: PSI, PSIT
  double precision, dimension (:),    allocatable :: xt, xp, xttp, xpp, xq;

  character (len = *), parameter :: file = "../Outputs/Flux/Stage1.nc"

  integer          :: err = 0, file_id
  integer          :: para_d_id, para_id, NPARA = 10, bound_d_id, rbound_id, zbound_id
  integer          :: lim_d_id, rlim_id, zlim_id, R_d_id, R_id, Z_d_id, Z_id
  integer          :: PS_d_id (2), PS_id, PN_id, T_id, P_id, TTp_id, Pp_id, Q_id
  integer          :: xt_id, xp_id, xttp_id, xpp_id, xq_id
  double precision :: para (10)

  MU0 = 16. * atan(1.0) * 1.e-7

  ! --------------
  ! Read EFIT file
  ! --------------
  open (unit = 100, file = '../Outputs/WriteEFIT/EFIT.txt', status = 'old')
  
  read (100, '(a48, 3i4)') string,  i3,      NRBOX,   NZBOX
  read (100, '(5e16.9)'  ) RBOXLEN, ZBOXLEN, R0,      RBOXLFT,  ZOFF
  read (100, '(5e16.9)'  ) RAXIS,   ZAXIS,   PSIAXIS, PSIBOUND, B0
  read (100, '(5e16.9)'  ) CURRENT, zero,    zero,    zero,     zero
  read (100, '(5e16.9)'  ) zero,    zero,    zero,    zero,     zero

  ! Use Raxis as scale major radius to ensure that normalized on-axis major radius is unity
  R0 = RAXIS

  RLEFT  =  RBOXLFT               /R0
  RRIGHT = (RBOXLFT + RBOXLEN)    /R0
  ZLOW   = (- ZBOXLEN /2. + ZOFF) /R0
  ZHIGH  = (+ ZBOXLEN /2. + ZOFF) /R0
  RAXIS  = RAXIS                  /R0
  ZAXIS  = ZAXIS                  /R0
   
  allocate (R (NRBOX))
  allocate (Z (NZBOX))

  RMAX = dble (NRBOX-1)
  do i = 1, NRBOX
     R (i) =  RLEFT + (RRIGHT - RLEFT) * dble (i-1) /RMAX
  end do

  ZMAX = dble (NZBOX-1)
  do j = 1, NZBOX
     Z (j) = ZLOW + (ZHIGH - ZLOW) * dble (j-1) /ZMAX
  end do

  allocate (PSIN  (NRBOX))
  allocate (T     (NRBOX))
  allocate (P     (NRBOX))
  allocate (TTp   (NRBOX))
  allocate (Pp    (NRBOX))
  allocate (Q     (NRBOX))
  allocate (PSI   (NRBOX, NZBOX))
  allocate (PSIT  (NZBOX, NRBOX))
  allocate (xt    (NRBOX))
  allocate (xp    (NRBOX))
  allocate (xttp  (NRBOX))
  allocate (xpp   (NRBOX))
  allocate (xq    (NRBOX))

  read (100, '(5e16.9)') ( T  (i),    i = 1, NRBOX)
  read (100, '(5e16.9)') ( P  (i),    i = 1, NRBOX)
  read (100, '(5e16.9)') ( TTp(i),    i = 1, NRBOX)
  read (100, '(5e16.9)') ( Pp (i),    i = 1, NRBOX)
  read (100, '(5e16.9)') ((PSI(i, j), i = 1, NRBOX), j = 1, NZBOX)
  read (100, '(5e16.9)') ( Q  (i),    i = 1, NRBOX)

  do i = 1, NRBOX
     PSIN (i) = dble (i-1) /RMAX
  end do

  read (100, '(2i5)') NBOUND, NLIM

  allocate (RBOUND (NBOUND))
  allocate (ZBOUND (NBOUND))
  allocate (RLIM   (NLIM))
  allocate (ZLIM   (NLIM))

  read (100, '(5e16.9)') (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  read (100, '(5e16.9)') (RLIM   (i), ZLIM   (i), i = 1, NLIM)

  do i = 1, NBOUND
     RBOUND (i) = RBOUND (i) /R0
     ZBOUND (i) = ZBOUND (i) /R0
  end do

  do i = 1, NLIM
     RLIM (i) = RLIM (i) /R0
     ZLIM (i) = ZLIM (i) /R0
  end do   

  close (unit = 100)
  
  para (1)  = R0
  para (2)  = dabs(B0)
  para (3)  = RLEFT
  para (4)  = RRIGHT
  para (5)  = ZLOW  
  para (6)  = ZHIGH 
  para (7)  = RAXIS
  para (8)  = ZAXIS 
  para (9)  = PSIAXIS  /R0/R0/dabs(B0)
  para (10) = PSIBOUND /R0/R0/dabs(B0)

  ! --------------------------------------
  ! Output equilibrium data to ascii files
  ! --------------------------------------  
  open  (unit = 101, file = '../Outputs/Flux/R0B0.txt')
  write (101, '(2e17.9)') R0, dabs(B0)
  close (unit = 101)

  open  (unit = 101, file = '../Outputs/Flux/Box.txt')
  write (101, '(4e17.9)') RLEFT, ZLOW, RRIGHT, ZHIGH
  close (unit = 101)

  open  (unit = 101, file = '../Outputs/Flux/Axis.txt')
  write (101, '(2e17.9)') RAXIS, ZAXIS
  close (unit = 101)
 
  open  (unit = 101, file = '../Outputs/Flux/R.txt')
  do i = 1, NRBOX
     write (101, '(1e17.9)') R (i)
  end do
  close (unit = 101)

  open  (unit = 101, file = '../Outputs/Flux/Z.txt')
  do j = 1, NZBOX
     write (101, '(1e17.9)') Z (j)
  end do
  close (unit = 101)
  
  open (unit = 101, file = '../Outputs/Flux/PsiSequential.txt')
  do i = 1, NRBOX
     do j = 1, NZBOX
         write (101, '(2i4,e17.9)') i, j, PSI (i, j)
     end do
  end do
  close (unit = 101)
  
  do i = 1, NRBOX
     do j = 1, NZBOX
        PSI  (i, j) = PSI (i, j) /R0/R0/dabs(B0)
        PSIT (j, i) = PSI (i, j)
      end do
  end do

  do i = 1, NRBOX
     xt   (i) = T   (i) /R0/B0
     xp   (i) = P   (i) *MU0/B0/B0
     xttp (i) = TTp (i) /dabs(B0)
     xpp  (i) = Pp  (i) *MU0*R0*R0/dabs(B0)
     xq   (i) = Q   (i)
  end do

  open  (unit = 101, file = '../Outputs/Flux/Profiles.txt')
  write (101, '(6e17.9)') (PSIN (i), xt (i), xp (i), xttp (i), xpp (i), xq (i), i = 1, NRBOX)
  close (unit = 101)
  
  open  (unit = 101, file = '../Outputs/Flux/Points.txt')
  write (101, '(4i5)') NRBOX, NZBOX, NBOUND, NLIM
  close (unit = 101)
 
  open  (unit = 101, file = '../Outputs/Flux/Boundary.txt')
  write (101, '(2e17.9)') (RBOUND (i), ZBOUND (i), i = 1, NBOUND)
  close (unit = 101)

  open  (unit = 101, file = '../Outputs/Flux/Limiter.txt')
  write (101, '(2e17.9)') (RLIM (i), ZLIM (i), i = 1, NLIM)
  close (unit = 101)
  
  ! -------------------------------------
  ! Output equilibrium data to NETCF file
  ! -------------------------------------
  err = nf90_create (file, NF90_CLOBBER, file_id)

  err = err + nf90_def_dim (file_id, "n_para",  NPARA,  para_d_id)
  err = err + nf90_def_dim (file_id, "i",       NRBOX,  R_d_id)
  err = err + nf90_def_dim (file_id, "j",       NZBOX,  Z_d_id)
  err = err + nf90_def_dim (file_id, "i_bound", NBOUND, bound_d_id)
  err = err + nf90_def_dim (file_id, "i_lim",   NLIM,   lim_d_id)

  PS_d_id = (/R_d_id, Z_d_id/) 
  
  err = err + nf90_def_var (file_id, "Parameters", NF90_DOUBLE, para_d_id,  para_id)
  err = err + nf90_def_var (file_id, "R",          NF90_DOUBLE, R_d_id,     R_id)
  err = err + nf90_def_var (file_id, "Z",          NF90_DOUBLE, Z_d_id,     Z_id)
  err = err + nf90_def_var (file_id, "PSI",        NF90_DOUBLE, PS_d_id,    PS_id)
  err = err + nf90_def_var (file_id, "PSI_N",      NF90_DOUBLE, R_d_id,     PN_id)
  err = err + nf90_def_var (file_id, "T",          NF90_DOUBLE, R_d_id,     T_id)
  err = err + nf90_def_var (file_id, "P",          NF90_DOUBLE, R_d_id,     P_id)
  err = err + nf90_def_var (file_id, "TTp",        NF90_DOUBLE, R_d_id,     TTp_id)
  err = err + nf90_def_var (file_id, "Pp",         NF90_DOUBLE, R_d_id,     Pp_id)
  err = err + nf90_def_var (file_id, "Q",          NF90_DOUBLE, R_d_id,     Q_id)
  err = err + nf90_def_var (file_id, "t",          NF90_DOUBLE, R_d_id,     xt_id)
  err = err + nf90_def_var (file_id, "p",          NF90_DOUBLE, R_d_id,     xp_id)
  err = err + nf90_def_var (file_id, "ttp",        NF90_DOUBLE, R_d_id,     xttp_id)
  err = err + nf90_def_var (file_id, "pp",         NF90_DOUBLE, R_d_id,     xpp_id)
  err = err + nf90_def_var (file_id, "q",          NF90_DOUBLE, R_d_id,     xq_id)
  err = err + nf90_def_var (file_id, "RBOUND",     NF90_DOUBLE, bound_d_id, rbound_id)
  err = err + nf90_def_var (file_id, "ZBOUND",     NF90_DOUBLE, bound_d_id, zbound_id)
  err = err + nf90_def_var (file_id, "RLIM",       NF90_DOUBLE, lim_d_id,   rlim_id)
  err = err + nf90_def_var (file_id, "ZLIM",       NF90_DOUBLE, lim_d_id,   zlim_id)

  err = err + nf90_enddef  (file_id)

  err = err + nf90_put_var (file_id, para_id,   para)
  err = err + nf90_put_var (file_id, R_id,      R)
  err = err + nf90_put_var (file_id, Z_id,      Z)
  err = err + nf90_put_var (file_id, PS_id,     PSIT)
  err = err + nf90_put_var (file_id, PN_id,     PSIN)
  err = err + nf90_put_var (file_id, T_id,      T)
  err = err + nf90_put_var (file_id, P_id,      P)
  err = err + nf90_put_var (file_id, TTp_id,    TTp)
  err = err + nf90_put_var (file_id, Pp_id,     Pp)
  err = err + nf90_put_var (file_id, Q_id,      Q)
  err = err + nf90_put_var (file_id, xt_id,     xt)
  err = err + nf90_put_var (file_id, xp_id,     xp)
  err = err + nf90_put_var (file_id, xttp_id,   xttp)
  err = err + nf90_put_var (file_id, xpp_id,    xpp)
  err = err + nf90_put_var (file_id, xq_id,     xq)
  err = err + nf90_put_var (file_id, rbound_id, RBOUND)
  err = err + nf90_put_var (file_id, zbound_id, ZBOUND)
  err = err + nf90_put_var (file_id, rlim_id,   RLIM)
  err = err + nf90_put_var (file_id, zlim_id,   ZLIM)
  
  err = err + nf90_close (file_id)

  if (err .NE. 0) then
     print*, "Error writing Outputs/Flux/Stage1.nc"
     stop
  end if
  
  ! --------
  ! Clean up
  ! --------
  deallocate (PSIN)
  deallocate (T)
  deallocate (P)
  deallocate (TTp)
  deallocate (Pp)
  deallocate (Q)
  deallocate (R)
  deallocate (Z)
  deallocate (RBOUND)
  deallocate (ZBOUND)
  deallocate (RLIM)
  deallocate (ZLIM)
  deallocate (PSI)
  deallocate (PSIT)
  deallocate (xt)
  deallocate (xp)
  deallocate (xttp)
  deallocate (xpp)
  deallocate (xq)

endsubroutine gFileRead
