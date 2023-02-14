subroutine check_int_stage()
   use const, only : PREC
   use const_phys, only : INVALID_JUDGE
   use var_potential, only : flg_stage, stage_sigma
   use var_top, only : nchains, nmp_chain, imp_chain
   use var_state, only : xyz

   implicit none
   integer :: imp, ichain, iinchain
   real(PREC) :: z, minz, stage_sigma3

   !! Testing molecule positions if stage energy is on
   if (flg_stage) then
      stage_sigma3 = 3*stage_sigma
      do ichain = 1, nchains
         minz = INVALID_JUDGE
         do iinchain = 1, nmp_chain(ichain)
            imp = imp_chain(iinchain, ichain)
            z = xyz(3,imp)
            if (z <= 0) then
               print '(a,i4,a,i7,a)', 'Error in initial structures. Move chain ', ichain, ' above the stage (resid ', iinchain, 's z is bellow 0)'
               flush(6)
               stop (2)
            elseif (z < minz) then
               minz = z
            endif
         enddo
         if (minz > stage_sigma3) then
            print '(a,i4,a)', 'Warning in initial structures. Chain ', ichain, ' is too far away from the stage.'
            flush(6)
         endif
      enddo
   endif
   !! Done testing

end subroutine check_int_stage
