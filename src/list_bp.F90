subroutine list_bp()

   use const, only : PREC
   use var_top, only : nmp
   use var_state, only : bp_status, ene_bp, for_bp, nt_bp_excess
   use var_potential, only : bp_model, nbp, bp_mp, bp_min_loop, bp_coef, &
                             bp3_map, bp3_dH, bp3_dS, bp_map
   use var_replica, only : nrep_proc

   implicit none
  
   integer :: n, imp, jmp
   integer :: ibp

   n = count(bp_map /= 0)

   allocate(nbp(nrep_proc))
   allocate(bp_mp(3, n, nrep_proc))
   allocate(bp_status(n, nrep_proc))
   allocate(ene_bp(n, nrep_proc))
   allocate(for_bp(3, 6, n))
   allocate(nt_bp_excess(nmp))

   nbp(:) = n
   bp_mp(:,:,:) = 0
   bp_status(:,:) = .False.
   ene_bp(:,:) = 0.0_PREC
   for_bp(:,:,:) = 0.0_PREC
   nt_bp_excess(:) = 0

   if (bp_model == 5) then
      allocate(bp_coef(2, n, nrep_proc))
      bp_coef(:,:,:) = 0.0_PREC
   endif

   do imp = 1, nmp
      do jmp = imp + bp_min_loop + 1, nmp

         if (bp_map(imp, jmp) == 0) cycle

         ibp = ibp + 1

         bp_mp(1, ibp, 1:nrep_proc) = imp
         bp_mp(2, ibp, 1:nrep_proc) = jmp
         bp_mp(3, ibp, 1:nrep_proc) = bp_map(imp, jmp)

         if (bp_model == 5) then
            print *, imp, jmp, bp_map(imp, jmp), bp3_map(imp, jmp)
            flush(6)
            bp_coef(1, ibp, 1:nrep_proc) = bp3_dH(bp3_map(imp, jmp))  ! dH
            bp_coef(2, ibp, 1:nrep_proc) = bp3_dS(bp3_map(imp, jmp))  ! dS
                        ! (0.001 already multiplied to dS so that the unit is kcal/mol/K)
         endif
      enddo
   enddo

   print '(a,i8)', '# nbp: ', nbp(1)

end subroutine list_bp
