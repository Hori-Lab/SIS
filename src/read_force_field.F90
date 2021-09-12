subroutine read_force_field(cfilepath)
      
   use const
   use const_phys
   use var_potential
   use var_io, only : iopen_hdl
  
   implicit none

   character(len=CHAR_FILE_PATH), intent(in) :: cfilepath
  
   character(len=CHAR_FILE_LINE) :: cline
   integer :: istat

   integer :: hdl

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   bond_k  = INVALID_VALUE
   bond_r0 = INVALID_VALUE
   angl_k  = INVALID_VALUE
   angl_t0 = INVALID_VALUE

   Ubp_cutoff = INVALID_VALUE
   Ubp0 = INVALID_VALUE
   Ubp_bond_k = INVALID_VALUE
   Ubp_bond_r = INVALID_VALUE
   Ubp_angl_k = INVALID_VALUE
   Ubp_angl_theta1 = INVALID_VALUE
   Ubp_angl_theta2 = INVALID_VALUE
   Ubp_dihd_k = INVALID_VALUE
   Ubp_dihd_phi1 = INVALID_VALUE
   Ubp_dihd_phi2 = INVALID_VALUE
   Ubp_min_loop = -1

   open(hdl, file=cfilepath, status='old', action='read', iostat=istat)

   do
      read (hdl, '(a)', iostat = istat) cline
      if(istat < 0) then
         exit
      else if(istat > 0) then
         write(*,*) 'Error: cannot read force field file'
         stop (2)
      end if
 
      cline = adjustl(cline)

      if (cline(1:6) == 'bond_k') read(cline(7:), *) bond_k
      if (cline(1:7) == 'bond_r0') read(cline(8:), *) bond_r0

      if (cline(1:6) == 'angl_k') read(cline(7:), *) angl_k
      if (cline(1:7) == 'angl_t0') read(cline(8:), *) angl_t0

      if (cline(1:10) == 'Ubp_cutoff') read(cline(11:), *) Ubp_cutoff
      if (cline(1:4) == 'Ubp0') read(cline(5:), *) Ubp0
      if (cline(1:10) == 'Ubp_bond_k') read(cline(11:), *) Ubp_bond_k
      if (cline(1:10) == 'Ubp_bond_r') read(cline(11:), *) Ubp_bond_r
      if (cline(1:10) == 'Ubp_angl_k') read(cline(11:), *) Ubp_angl_k
      if (cline(1:15) == 'Ubp_angl_theta1') read(cline(16:), *) Ubp_angl_theta1
      if (cline(1:15) == 'Ubp_angl_theta2') read(cline(16:), *) Ubp_angl_theta2
      if (cline(1:10) == 'Ubp_dihd_k') read(cline(11:), *) Ubp_dihd_k
      if (cline(1:13) == 'Ubp_dihd_phi1') read(cline(14:), *) Ubp_dihd_phi1
      if (cline(1:13) == 'Ubp_dihd_phi2') read(cline(14:), *) Ubp_dihd_phi2
      if (cline(1:12) == 'Ubp_min_loop') read(cline(13:), *) Ubp_min_loop

   enddo

   close(hdl)
   iopen_hdl = iopen_hdl - 1

   ! Check
   if (bond_k  > INVALID_JUDGE) then
      write(*,*) "INVALID bond_k in the force field file"
   else
      write(*,*) "bond_k", bond_k
   endif

   if (bond_r0 > INVALID_JUDGE) then
      write(*,*) "INVALID bond_r0 in the force field file"
   else
      write(*,*) "bond_r0", bond_r0
   endif

   if (angl_k  > INVALID_JUDGE) then
      write(*,*) "INVALID angl_k in the force field file"
   else
      write(*,*) "angl_k", angl_k
   endif

   if (angl_t0 > INVALID_JUDGE) then
      write(*,*) "INVALID angl_t0 in the force field file"
   else
      write(*,*) "angl_t0", angl_t0 
   endif

   if (Ubp_cutoff > INVALID_JUDGE) then
      write(*,*) "INVALID Ubp_cutoff in the force field file"
   else
      write(*,*) "Ubp_cutoff", Ubp_cutoff
   endif

   if (Ubp0 > INVALID_JUDGE) then
      write(*,*) "INVALID Ubp0 in the force field file"
   else
      write(*,*) "Ubp0", Ubp0
   endif

   if (Ubp_bond_k > INVALID_JUDGE) then
      write(*,*) "INVALID Ubp_bond_k in the force field file"
   else
      write(*,*) "Ubp_bond_k", Ubp_bond_k
   endif

   if (Ubp_bond_r > INVALID_JUDGE) then
      write(*,*) "INVALID Ubp_bond_r in the force field file"
   else
      write(*,*) "Ubp_bond_r", Ubp_bond_r
   endif

   if (Ubp_angl_k > INVALID_JUDGE) then
      write(*,*) "INVALID Ubp_angl_k in the force field file"
   else
      write(*,*) "Ubp_angl_k", Ubp_angl_k
   endif

   if (Ubp_angl_theta1 > INVALID_JUDGE) then
      write(*,*) "INVALID Ubp_angl_theta1 in the force field file"
   else
      write(*,*) "Ubp_angl_theta1", Ubp_angl_theta1
   endif

   if (Ubp_angl_theta2 > INVALID_JUDGE) then
      write(*,*) "INVALID Ubp_angl_theta2 in the force field file"
   else
      write(*,*) "Ubp_angl_theta2", Ubp_angl_theta2
   endif

   if (Ubp_dihd_k > INVALID_JUDGE) then
      write(*,*) "INVALID Ubp_dihd_k in the force field file"
   else
      write(*,*) "Ubp_dihd_k", Ubp_dihd_k
   endif

   if (Ubp_dihd_phi1 > INVALID_JUDGE) then
      write(*,*) "INVALID Ubp_dihd_phi1 in the force field file"
   else
      write(*,*) "Ubp_dihd_phi1", Ubp_dihd_phi1
   endif

   if (Ubp_dihd_phi2 > INVALID_JUDGE) then
      write(*,*) "INVALID Ubp_dihd_phi2 in the force field file"
   else
      write(*,*) "Ubp_dihd_phi2", Ubp_dihd_phi2
   endif

   if (Ubp_min_loop < 0) then
      write(*,*) "INVALID Ubp_min_loop in the force field file"
   else
      write(*,*) "Ubp_min_loop", Ubp_dihd_phi2
   endif


end subroutine read_force_field
