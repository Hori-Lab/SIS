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
   type(toml_table), pointer :: group, node, subnode
   !type(toml_array), pointer :: array
   !======= 

   stat = .False.

   print '(2a)', 'Reading force-field file: ', trim(cfile_ff)
   flush(6)

   call set_invalid()

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   open(hdl, file=cfile_ff, status='old', action='read', iostat=istat)

   if (istat /= 0) then
      print '(2a)', 'Error: failed to open the force-field file. ', trim(cfile_ff)
      flush(6)
      error stop
   endif

   call toml_parse(table, hdl)

   close(hdl)
   iopen_hdl = iopen_hdl - 1

   call get_value(table, "title", cline)
   print '(2a)', '# title: ', trim(cline)
   flush(6)

   call get_value(table, "potential", group)
   if (.not. associated(group)) then
      print '(a)', 'Error: [potential] required in FF file'
      return
   endif

   call get_value(group, "bond", node)
   if (associated(node)) then
      call get_value(node, "k", bond_k, stat=istat)
      call get_value(node, "r0", bond_r0)

   else
      print '(a)', 'Error: [bond] parameters required in FF file'
      return
   endif

   call get_value(group, "angle", node)
   if (associated(node)) then
      call get_value(node, "k", angl_k)
      call get_value(node, "a0", angl_t0)
      
   call get_value(group, "BAdihedral", node)
   if (associated(node)) then
      call get_value(node, "k", angl_kphi)
      call get_value(node, "phi0", angl_phi0)

   else
      print '(a)', 'Error: [angle] parameters required in FF file'
      return
   endif

   call get_value(group, "basepair", node)
   bp_cutoff_energy = 0.01_PREC  ! Default /kcal/mol
   if (associated(node)) then
      !call get_value(node, "min_loop", bp_min_loop)
      !call get_value(node, "cutoff", bp_cutoff_dist)
      call get_value(node, "cutoff_energy", bp_cutoff_energy)

      ! older format
      !call get_value(node, "bond_k", bp_bond_k)
      !call get_value(node, "bond_r", bp_bond_r)
      !call get_value(node, "angl_k", bp_angl_k)
      !call get_value(node, "angl_theta1", bp_angl_theta1)
      !call get_value(node, "angl_theta2", bp_angl_theta2)
      !call get_value(node, "dihd_k", bp_dihd_k)
      !call get_value(node, "dihd_phi1", bp_dihd_phi1)
      !call get_value(node, "dihd_phi2", bp_dihd_phi2)

      !bp_seqdep = 0
      !! = 0 (Default): No sequence dependence. Only U0_GC, U0_AU, U0_GU are required.
      !! = 1: Sequence dependent parameters. All possible combinations of trinucleotide-dimer are required.
      !call get_value(node, "seqdep", bp_seqdep)

      !if (bp_seqdep == 0) then
      !   call get_value(node, "U0_GC", bp_U0_GC)
      !   call get_value(node, "U0_AU", bp_U0_AU)
      !   call get_value(node, "U0_GU", bp_U0_GU)
      !endif

      call get_value(node, "GC", subnode)
      if (associated(subnode)) then
         call get_value(subnode, "U0", bp_paras(BPT%GC)%U0)
         call get_value(subnode, "bond_k", bp_paras(BPT%GC)%bond_k)
         call get_value(subnode, "bond_r", bp_paras(BPT%GC)%bond_r)
         call get_value(subnode, "angl_k1", bp_paras(BPT%GC)%angl_k1)
         call get_value(subnode, "angl_k2", bp_paras(BPT%GC)%angl_k2)
         call get_value(subnode, "angl_theta1", bp_paras(BPT%GC)%angl_theta1)
         call get_value(subnode, "angl_theta2", bp_paras(BPT%GC)%angl_theta2)
         call get_value(subnode, "dihd_k1", bp_paras(BPT%GC)%dihd_k1)
         call get_value(subnode, "dihd_k2", bp_paras(BPT%GC)%dihd_k2)
         call get_value(subnode, "dihd_phi1", bp_paras(BPT%GC)%dihd_phi1)
         call get_value(subnode, "dihd_phi2", bp_paras(BPT%GC)%dihd_phi2)
      endif

      call get_value(node, "AU", subnode)
      if (associated(subnode)) then
         call get_value(subnode, "U0", bp_paras(BPT%AU)%U0)
         call get_value(subnode, "bond_k", bp_paras(BPT%AU)%bond_k)
         call get_value(subnode, "bond_r", bp_paras(BPT%AU)%bond_r)
         call get_value(subnode, "angl_k1", bp_paras(BPT%AU)%angl_k1)
         call get_value(subnode, "angl_k2", bp_paras(BPT%AU)%angl_k2)
         call get_value(subnode, "angl_theta1", bp_paras(BPT%AU)%angl_theta1)
         call get_value(subnode, "angl_theta2", bp_paras(BPT%AU)%angl_theta2)
         call get_value(subnode, "dihd_k1", bp_paras(BPT%AU)%dihd_k1)
         call get_value(subnode, "dihd_k2", bp_paras(BPT%AU)%dihd_k2)
         call get_value(subnode, "dihd_phi1", bp_paras(BPT%AU)%dihd_phi1)
         call get_value(subnode, "dihd_phi2", bp_paras(BPT%AU)%dihd_phi2)
      endif

      call get_value(node, "GU", subnode)
      if (associated(subnode)) then
         call get_value(subnode, "U0", bp_paras(BPT%GU)%U0)
         call get_value(subnode, "bond_k", bp_paras(BPT%GU)%bond_k)
         call get_value(subnode, "bond_r", bp_paras(BPT%GU)%bond_r)
         call get_value(subnode, "angl_k1", bp_paras(BPT%GU)%angl_k1)
         call get_value(subnode, "angl_k2", bp_paras(BPT%GU)%angl_k2)
         call get_value(subnode, "angl_theta1", bp_paras(BPT%GU)%angl_theta1)
         call get_value(subnode, "angl_theta2", bp_paras(BPT%GU)%angl_theta2)
         call get_value(subnode, "dihd_k1", bp_paras(BPT%GU)%dihd_k1)
         call get_value(subnode, "dihd_k2", bp_paras(BPT%GU)%dihd_k2)
         call get_value(subnode, "dihd_phi1", bp_paras(BPT%GU)%dihd_phi1)
         call get_value(subnode, "dihd_phi2", bp_paras(BPT%GU)%dihd_phi2)
      endif

   else
      print '(a)', 'Error: [basepair] parameters required in FF file'
      return
   endif

   call get_value(group, "wca", node)
   if (associated(node)) then
      call get_value(node, "sigma", wca_sigma)
      call get_value(node, "epsilon", wca_eps)

   else
      print '(a)', 'Error: [wca] parameters required in FF file'
      return
   endif

   call check(stat)  ! stat will be .True. if everything is ok

   if (stat) then
      print '(a)', 'Done: reading force-field file'
      print *
   endif

contains

   subroutine set_invalid()
      integer :: bptype

      bond_k  = INVALID_VALUE
      bond_r0 = INVALID_VALUE

      angl_k  = INVALID_VALUE
      angl_t0 = INVALID_VALUE
   
      !bp_min_loop = -1
      !bp_cutoff_dist = INVALID_VALUE
      !bp_cutoff_energy = INVALID_VALUE
      do bptype = 1, BPT%MAX
         bp_paras(bptype)%U0 = INVALID_VALUE
         bp_paras(bptype)%bond_k = INVALID_VALUE
         bp_paras(bptype)%bond_r = INVALID_VALUE
         bp_paras(bptype)%angl_k1 = INVALID_VALUE
         bp_paras(bptype)%angl_k2 = INVALID_VALUE
         bp_paras(bptype)%angl_theta1 = INVALID_VALUE
         bp_paras(bptype)%angl_theta2 = INVALID_VALUE
         bp_paras(bptype)%dihd_k1 = INVALID_VALUE
         bp_paras(bptype)%dihd_k2 = INVALID_VALUE
         bp_paras(bptype)%dihd_phi1 = INVALID_VALUE
         bp_paras(bptype)%dihd_phi2 = INVALID_VALUE
      enddo
   
      wca_sigma = INVALID_VALUE
      wca_eps = INVALID_VALUE
   endsubroutine set_invalid

   subroutine check_bp_paras(bptype)

      use const_idx, only : BPTYPE_CHAR
      integer, intent(in) :: bptype

      if (bp_paras(bptype)%U0 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp U0 for ", BPTYPE_CHAR(bptype)
         stat = .False.
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " U0: ", bp_paras(bptype)%U0
      endif

      if (bp_paras(bptype)%bond_k > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp bond_k for ", BPTYPE_CHAR(bptype)
         stat = .False.
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " bond_k: ", bp_paras(bptype)%bond_k
      endif

      if (bp_paras(bptype)%bond_r > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp bond_r for ", BPTYPE_CHAR(bptype)
         stat = .False.
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " bond_r: ", bp_paras(bptype)%bond_r
      endif

      if (bp_paras(bptype)%angl_k1 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp angl_k1 for ", BPTYPE_CHAR(bptype)
         stat = .False.
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " angl_k1: ", bp_paras(bptype)%angl_k1
      endif

      if (bp_paras(bptype)%angl_k2 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp angl_k2 for ", BPTYPE_CHAR(bptype)
         stat = .False.
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " angl_k2: ", bp_paras(bptype)%angl_k2
      endif

      if (bp_paras(bptype)%angl_theta1 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp angl_theta1 for ", BPTYPE_CHAR(bptype)
         stat = .False.
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " angl_theta1: ", bp_paras(bptype)%angl_theta1
      endif

      if (bp_paras(bptype)%angl_theta2 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp angl_theta2 for ", BPTYPE_CHAR(bptype)
         stat = .False.
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " angl_theta2: ", bp_paras(bptype)%angl_theta2
      endif

      if (bp_paras(bptype)%dihd_k1 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp dihd_k1 for ", BPTYPE_CHAR(bptype)
         stat = .False.
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " dihd_k1: ", bp_paras(bptype)%dihd_k1
      endif

      if (bp_paras(bptype)%dihd_k2 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp dihd_k2 for ", BPTYPE_CHAR(bptype)
         stat = .False.
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " dihd_k2: ", bp_paras(bptype)%dihd_k2
      endif

      if (bp_paras(bptype)%dihd_phi1 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp dihd_phi1 for ", BPTYPE_CHAR(bptype)
         stat = .False.
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " dihd_phi1: ", bp_paras(bptype)%dihd_phi1
      endif

      if (bp_paras(bptype)%dihd_phi2 > INVALID_JUDGE) then
         print '(2a)', "INVALID value in the force field file. bp dihd_phi2 for ", BPTYPE_CHAR(bptype)
         stat = .False.
      else
         print '(3a,g15.8)', "# ", BPTYPE_CHAR(bptype), " dihd_phi2: ", bp_paras(bptype)%dihd_phi2
      endif

   endsubroutine check_bp_paras

   subroutine check(stat)

      logical, intent(inout) :: stat

      stat = .True.

      ! Check
      if (bond_k  > INVALID_JUDGE) then
         print '(a)', "INVALID bond_k in the force field file"
         stat = .False.
      else
         print '(a,g15.8)', "# bond_k: ", bond_k
      endif

      if (bond_r0 > INVALID_JUDGE) then
         print '(a)', "INVALID bond_r0 in the force field file"
         stat = .False.
      else
         print '(a,g15.8)', "# bond_r0: ", bond_r0
      endif

      if (angl_k  > INVALID_JUDGE) then
         print '(a)', "INVALID angl_k in the force field file"
         stat = .False.
      else
         print '(a,g15.8)', "# angl_k: ", angl_k
      endif

      if (angl_t0 > INVALID_JUDGE) then
         print '(a)', "INVALID angl_t0 in the force field file"
         stat = .False.
      else
         print '(a,g15.8)', "# angl_t0: ", angl_t0 
      endif

      !if (bp_cutoff_dist > INVALID_JUDGE) then
      !   print '(a)', "INVALID bp_cutoff in the force field file"
      !   stat = .False.
      !else
      !   print '(a)', "# bp_cutoff: ", bp_cutoff_dist
      !endif
      print '(a,g15.8)', '# bp_cutoff_energy: ', bp_cutoff_energy

!      if (bp_seqdep == 0) then
!         if (bp_U0_GC > INVALID_JUDGE) then
!            print '(a)', "INVALID bp_U0_GC in the force field file"
!            stat = .False.
!         else
!            print '(a)', "# bp_U0_GC: ", bp_U0_GC
!         endif
!
!         if (bp_U0_AU > INVALID_JUDGE) then
!            print '(a)', "INVALID bp_U0_AU in the force field file"
!            stat = .False.
!         else
!            print '(a)', "# bp_U0_AU: ", bp_U0_AU
!         endif
!
!         if (bp_U0_GU > INVALID_JUDGE) then
!            print '(a)', "INVALID bp_U0_GU in the force field file"
!            stat = .False.
!         else
!            print '(a)', "# bp_U0_GU: ", bp_U0_GU
!         endif
!
!      else if (bp_seqdep == 1) then
!         continue
!
!      else
!         print '(a)', "INVALID bp_seqdep in the force field file"
!         stat = .False.
!      endif

      call check_bp_paras(BPT%GC)
      call check_bp_paras(BPT%AU)
      call check_bp_paras(BPT%GU)

      !if (bp_min_loop < 0) then
      !   print '(a)', "INVALID bp_min_loop in the force field file"
      !   stat = .False.
      !else
      !   print '(a)', "# bp_min_loop: ", bp_min_loop
      !endif

      if (wca_sigma > INVALID_JUDGE) then
         print '(a)', "INVALID wca_sigma in the force field file"
         stat = .False.
      else
         print '(a,g15.8)', "# wca_sigma: ", wca_sigma
      endif

      if (wca_eps > INVALID_JUDGE) then
         print '(a)', "INVALID wca_eps in the force field file"
         stat = .False.
      else
         print '(a,g15.8)', "# wca_eps: ", wca_eps
      endif

   endsubroutine check

end subroutine read_force_field
