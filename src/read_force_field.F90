subroutine read_force_field(stat)

   use tomlf
      
   use const_phys, only : INVALID_VALUE, INVALID_JUDGE
   use var_potential
   use var_io, only : iopen_hdl, cfile_ff
  
   implicit none

   logical, intent(out) :: stat
  
   integer :: istat
   integer :: hdl

   character(len=:), allocatable :: cline

   !======= TOML
   type(toml_table), allocatable :: table
   type(toml_table), pointer :: group, node
   !type(toml_array), pointer :: array
   !======= 

   stat = .False.

   write(*,*) 'Reading force-field file: '//trim(cfile_ff)

   call set_invalid()

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   open(hdl, file=cfile_ff, status='old', action='read', iostat=istat)

   if (istat /= 0) then
      write(*,*) 'Error: failed to open the force-field file. '//trim(cfile_ff)
      stop (2)
   endif

   call toml_parse(table, hdl)

   close(hdl)
   iopen_hdl = iopen_hdl - 1

   call get_value(table, "title", cline)
   write(*,*) '# title: ', trim(cline)

   call get_value(table, "potential", group)
   if (.not. associated(group)) then
      write(*,*) 'Error: [potential] required in FF file'
      return
   endif

   call get_value(group, "bond", node)
   if (associated(node)) then
      call get_value(node, "k", bond_k, stat=istat)
      call get_value(node, "r0", bond_r0)

   else
      write(*,*) 'Error: [bond] parameters required in FF file'
      return
   endif

   call get_value(group, "angle", node)
   if (associated(node)) then
      call get_value(node, "k", angl_k)
      call get_value(node, "a0", angl_t0)

   else
      write(*,*) 'Error: [angle] parameters required in FF file'
      return
   endif

   call get_value(group, "basepair", node)
   if (associated(node)) then
      call get_value(node, "min_loop", bp_min_loop)
      call get_value(node, "cutoff", bp_cutoff)
      call get_value(node, "bond_k", bp_bond_k)
      call get_value(node, "bond_r", bp_bond_r)
      call get_value(node, "angl_k", bp_angl_k)
      call get_value(node, "angl_theta1", bp_angl_theta1)
      call get_value(node, "angl_theta2", bp_angl_theta2)
      call get_value(node, "dihd_k", bp_dihd_k)
      call get_value(node, "dihd_phi1", bp_dihd_phi1)
      call get_value(node, "dihd_phi2", bp_dihd_phi2)

      bp_seqdep = 0
      ! = 0 (Default): No sequence dependence. Only U0_GC, U0_AU, U0_GU are required.
      ! = 1: Sequence dependent parameters. All possible combinations of trinucleotide-dimer are required.
      call get_value(node, "seqdep", bp_seqdep)

      if (bp_seqdep == 0) then
         call get_value(node, "U0_GC", bp_U0_GC)
         call get_value(node, "U0_AU", bp_U0_AU)
         call get_value(node, "U0_GU", bp_U0_GU)
      endif

   else
      write(*,*) 'Error: [basepair] parameters required in FF file'
      return
   endif

   call get_value(group, "wca", node)
   if (associated(node)) then
      call get_value(node, "sigma", wca_sigma)
      call get_value(node, "epsilon", wca_eps)

   else
      write(*,*) 'Error: [wca] parameters required in FF file'
      return
   endif

   call check()

   write(*,*) 'Done: reading force-field file'
   write(*,*)

   stat = .True.

contains

   subroutine set_invalid()
      bond_k  = INVALID_VALUE
      bond_r0 = INVALID_VALUE

      angl_k  = INVALID_VALUE
      angl_t0 = INVALID_VALUE
   
      bp_min_loop = -1
      bp_cutoff = INVALID_VALUE
      bp_U0_GC = INVALID_VALUE
      bp_U0_AU = INVALID_VALUE
      bp_U0_GU = INVALID_VALUE
      bp_bond_k = INVALID_VALUE
      bp_bond_r = INVALID_VALUE
      bp_angl_k = INVALID_VALUE
      bp_angl_theta1 = INVALID_VALUE
      bp_angl_theta2 = INVALID_VALUE
      bp_dihd_k = INVALID_VALUE
      bp_dihd_phi1 = INVALID_VALUE
      bp_dihd_phi2 = INVALID_VALUE
   
      wca_sigma = INVALID_VALUE
      wca_eps = INVALID_VALUE
   endsubroutine set_invalid

   subroutine check
      ! Check
      if (bond_k  > INVALID_JUDGE) then
         write(*,*) "INVALID bond_k in the force field file"
         return
      else
         write(*,*) "# bond_k: ", bond_k
      endif

      if (bond_r0 > INVALID_JUDGE) then
         write(*,*) "INVALID bond_r0 in the force field file"
         return
      else
         write(*,*) "# bond_r0: ", bond_r0
      endif

      if (angl_k  > INVALID_JUDGE) then
         write(*,*) "INVALID angl_k in the force field file"
         return
      else
         write(*,*) "# angl_k: ", angl_k
      endif

      if (angl_t0 > INVALID_JUDGE) then
         write(*,*) "INVALID angl_t0 in the force field file"
         return
      else
         write(*,*) "# angl_t0: ", angl_t0 
      endif

      if (bp_cutoff > INVALID_JUDGE) then
         write(*,*) "INVALID bp_cutoff in the force field file"
         return
      else
         write(*,*) "# bp_cutoff: ", bp_cutoff
      endif

      if (bp_seqdep == 0) then
         if (bp_U0_GC > INVALID_JUDGE) then
            write(*,*) "INVALID bp_U0_GC in the force field file"
            return
         else
            write(*,*) "# bp_U0_GC: ", bp_U0_GC
         endif

         if (bp_U0_AU > INVALID_JUDGE) then
            write(*,*) "INVALID bp_U0_AU in the force field file"
            return
         else
            write(*,*) "# bp_U0_AU: ", bp_U0_AU
         endif

         if (bp_U0_GU > INVALID_JUDGE) then
            write(*,*) "INVALID bp_U0_GU in the force field file"
            return
         else
            write(*,*) "# bp_U0_GU: ", bp_U0_GU
         endif

      else if (bp_seqdep == 1) then
         continue

      else
         write(*,*) "INVALID bp_seqdep in the force field file"
         return
      endif

      if (bp_bond_k > INVALID_JUDGE) then
         write(*,*) "INVALID bp_bond_k in the force field file"
         return
      else
         write(*,*) "# bp_bond_k: ", bp_bond_k
      endif

      if (bp_bond_r > INVALID_JUDGE) then
         write(*,*) "INVALID bp_bond_r in the force field file"
         return
      else
         write(*,*) "# bp_bond_r: ", bp_bond_r
      endif

      if (bp_angl_k > INVALID_JUDGE) then
         write(*,*) "INVALID bp_angl_k in the force field file"
         return
      else
         write(*,*) "# bp_angl_k: ", bp_angl_k
      endif

      if (bp_angl_theta1 > INVALID_JUDGE) then
         write(*,*) "INVALID bp_angl_theta1 in the force field file"
         return
      else
         write(*,*) "# bp_angl_theta1: ", bp_angl_theta1
      endif

      if (bp_angl_theta2 > INVALID_JUDGE) then
         write(*,*) "INVALID bp_angl_theta2 in the force field file"
         return
      else
         write(*,*) "# bp_angl_theta2: ", bp_angl_theta2
      endif

      if (bp_dihd_k > INVALID_JUDGE) then
         write(*,*) "INVALID bp_dihd_k in the force field file"
         return
      else
         write(*,*) "# bp_dihd_k: ", bp_dihd_k
      endif

      if (bp_dihd_phi1 > INVALID_JUDGE) then
         write(*,*) "INVALID bp_dihd_phi1 in the force field file"
         return
      else
         write(*,*) "# bp_dihd_phi1: ", bp_dihd_phi1
      endif

      if (bp_dihd_phi2 > INVALID_JUDGE) then
         write(*,*) "INVALID bp_dihd_phi2 in the force field file"
         return
      else
         write(*,*) "# bp_dihd_phi2: ", bp_dihd_phi2
      endif

      if (bp_min_loop < 0) then
         write(*,*) "INVALID bp_min_loop in the force field file"
         return
      else
         write(*,*) "# bp_min_loop: ", bp_min_loop
      endif

      if (wca_sigma > INVALID_JUDGE) then
         write(*,*) "INVALID wca_sigma in the force field file"
         return
      else
         write(*,*) "# wca_sigma: ", wca_sigma
      endif

      if (wca_eps > INVALID_JUDGE) then
         write(*,*) "INVALID wca_eps in the force field file"
         return
      else
         write(*,*) "# wca_eps: ", wca_eps
      endif
   endsubroutine check

end subroutine read_force_field
