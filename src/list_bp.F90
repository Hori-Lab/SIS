subroutine list_bp()

   use const, only : PREC
   use const_idx, only : MOLT
   use var_top, only : nmp, lmp_mp, moltypes, ichain_mp, nmp_chain
   use var_state, only : bp_status, bp_status_MC, ene_bp, for_bp, nt_bp_excess
   use var_potential, only : bp_model, nbp, bp_mp, bp_min_loop, bp_coef, &
                             bp3_map, bp3_dH, bp3_dS, bp_map
   use var_replica, only : nrep_proc

   implicit none
  
   integer :: n, imp, jmp
   integer :: ichain, jchain
   logical :: i_circ, j_circ
   integer :: ibp

   n = count(bp3_map /= 0)

   allocate(nbp(nrep_proc))
   allocate(bp_mp(7, n, nrep_proc))
   allocate(bp_status(n, nrep_proc))
   allocate(ene_bp(n, nrep_proc))
   allocate(for_bp(3, 6, n))
   allocate(nt_bp_excess(nmp))
   allocate(bp_status_MC(n, nrep_proc))

   nbp(:) = n
   bp_mp(:,:,:) = 0
   bp_status(:,:) = .False.
   ene_bp(:,:) = 0.0_PREC
   for_bp(:,:,:) = 0.0_PREC
   nt_bp_excess(:) = 0
   bp_status_MC(:,:) = .False.

   if (bp_model == 5) then
      allocate(bp_coef(2, n, nrep_proc))
      bp_coef(:,:,:) = 0.0_PREC
   endif

   ibp = 0
   do imp = 1, nmp
      do jmp = imp + bp_min_loop + 1, nmp

         if (bp3_map(imp, jmp) == 0) cycle

         ibp = ibp + 1

         bp_mp(1, ibp, 1:nrep_proc) = imp
         bp_mp(2, ibp, 1:nrep_proc) = jmp

         ichain = ichain_mp(imp)
         jchain = ichain_mp(jmp)
         i_circ = .False.
         j_circ = .False.
         if (moltypes(ichain) == MOLT%CIRCRNA) i_circ = .True.
         if (moltypes(jchain) == MOLT%CIRCRNA) j_circ = .True.

         ! If circular and imp is the first nucleotide,
         ! imp-1 should be the last nucleotide of the chain.
         if (i_circ .and. lmp_mp(imp) == 1) then
            bp_mp(3, ibp, 1:nrep_proc) = imp + nmp_chain(ichain) - 1
         else
            bp_mp(3, ibp, 1:nrep_proc) = imp - 1
         endif
         if (j_circ .and. lmp_mp(jmp) == 1) then
            bp_mp(4, ibp, 1:nrep_proc) = jmp + nmp_chain(jchain) - 1
         else
            bp_mp(4, ibp, 1:nrep_proc) = jmp - 1
         endif

         ! If circular and imp is the last nucleotide of the chain,
         ! imp+1 should be the first nucleotide.
         if (i_circ .and. lmp_mp(imp) == nmp_chain(ichain)) then
            bp_mp(5, ibp, 1:nrep_proc) = imp - nmp_chain(ichain) + 1
         else
            bp_mp(5, ibp, 1:nrep_proc) = imp + 1
         endif
         if (j_circ .and. lmp_mp(jmp) == nmp_chain(jchain)) then
            bp_mp(6, ibp, 1:nrep_proc) = jmp - nmp_chain(jchain) + 1
         else
            bp_mp(6, ibp, 1:nrep_proc) = jmp + 1
         endif

         bp_mp(7, ibp, 1:nrep_proc) = bp_map(imp, jmp)

         if (bp_model == 5) then
            bp_coef(1, ibp, 1:nrep_proc) = bp3_dH(bp3_map(imp, jmp))  ! dH
            bp_coef(2, ibp, 1:nrep_proc) = bp3_dS(bp3_map(imp, jmp))  ! dS
                        ! (0.001 already multiplied to dS so that the unit is kcal/mol/K)
         endif
      enddo
   enddo

   print '(a,i8)', '# nbp: ', nbp(1)

end subroutine list_bp
