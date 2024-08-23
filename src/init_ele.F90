subroutine init_ele()
   
   use, intrinsic :: iso_fortran_env, only : output_unit
   use const, only : PREC
   use const_idx, only : REPT, SEQT
   use var_top, only : nmp, has_charge, inp_no_charge, dummy_has_charge, &
                       nchains, imp_chain, nmp_chain, seq
   use var_state, only : tempK, lambdaD, diele, ionic_strength, &
                         temp_independent, diele_dTcoef
   use var_potential, only : ele_coef, ele_cutoff_type, ele_cutoff_inp, ele_cutoff
   use var_replica, only : flg_replica, flg_repvar, rep2val, nrep_proc, irep2grep

   implicit none

   integer :: i, ichain, imp, irep, grep
   real(PREC) :: tK, lb, Zp, ionstr

   interface
   subroutine set_ele(irep, tempk, ionstr, out_lb, out_Zp)
      use const, only : PREC
      integer, intent(in) :: irep
      real(PREC), intent(in) :: tempk
      real(PREC), intent(in) :: ionstr
      real(PREC), intent(out), optional :: out_lb
      real(PREC), intent(out), optional :: out_Zp
   endsubroutine set_ele
   endinterface

   !! Allocation
   allocate(has_charge(nmp))
   allocate(lambdaD(nrep_proc))
   allocate(diele(nrep_proc))
   allocate(ele_coef(nrep_proc))
   !allocate(charge(nmp, nrep_proc))
   allocate(ele_cutoff(nrep_proc))

   if (temp_independent > 0) then
      allocate(diele_dTcoef(nrep_proc))
      diele_dTcoef(:) = 0.0_PREC
   endif

   lambdaD(:) = 0.0_PREC
   diele(:) = 0.0_PREC
   ele_coef(:) = 0.0_PREC
   !charge(:,:) = 0.0_PREC
   ele_cutoff(:) = 0.0_PREC

   tK = tempK
   ionstr = ionic_strength

   has_charge(:) = .True.

   ! "dummy_has_charge" option (default = False)
   if (.not. dummy_has_charge) then
      do ichain =  1, nchains
         do i = 1, nmp_chain(ichain)
            if (seq(i, ichain) == SEQT%D) then
               imp = imp_chain(i, ichain)
               has_charge(imp) = .False.
            endif
         enddo
      enddo
   endif

   ! "no_charge" option
   if (allocated(inp_no_charge)) then
      do i = 1, size(inp_no_charge)
         has_charge(inp_no_charge(i)) = .False.
      enddo
   endif

   ! Set cutoff
   if (ele_cutoff_type == 1) then
      ele_cutoff(:) = ele_cutoff_inp
   else if (ele_cutoff_type == 2) then
      continue   ! This will be calcualted in the next loop.
   else
      error stop "Error: Invalid ele_cutoff_type value in set_ele."
   endif


   print '(a)', 'Set electrostaic parameters'

   do irep = 1, nrep_proc

      if (flg_replica) then
         grep = irep2grep(irep)

         if (flg_repvar(REPT%TEMP)) then
            tK = rep2val(grep, REPT%TEMP)
         endif

         if (flg_repvar(REPT%ION)) then
            ionstr = rep2val(grep, REPT%ION)
         endif
      endif

      call set_ele(irep, tK, ionstr, lb, Zp)

      if (ele_cutoff_type == 2) then
         ele_cutoff(irep) = ele_cutoff_inp * lambdaD(irep)
      endif

      print '(a,i8)', '# Replica: ', irep
      print '(a,g15.8)', '#    Ionic strength: ', ionstr
      print '(a,g15.8)', '#    Dielectric constant (H2O): ', diele(irep)
      print '(a,g15.8)', '#    coef: ', ele_coef(irep) / (Zp**2)
      print '(a,g15.8)', '#    Bjerrum length: ', lb
      print '(a,g15.8)', '#    Debye length (lambdaD): ', lambdaD(irep)
      print '(a,g15.8)', '#    Reduced charge on each phosphate: ', -Zp
      print '(a,g15.8)', '#    ele_cutoff: ', ele_cutoff(irep)

   enddo

   print *
   flush(output_unit)

endsubroutine init_ele
