subroutine set_bp_map()

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const
   use const_idx, only : SEQT, seqt2char, seqt2nnt
   use var_io, only : iopen_hdl, cfile_prefix
   use var_top, only : nmp, seq, lmp_mp, ichain_mp, nmp_chain
   use var_potential, only : bp_model, bp_map, bp_map_0, bp_map_dG, &
                             NN_dH, NN_dS, dH0, dS0, coef_dG
   use var_state, only : tempK, temp_independent

   implicit none

   integer :: imp, jmp
   integer :: i, j, ichain, jchain
   integer :: hdl, istat
   real(PREC) :: dG, dH, dS

   bp_map(:,:) = bp_map_0(:,:)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! All pairwise
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (bp_model == 5) then

      iopen_hdl = iopen_hdl + 1
      hdl = iopen_hdl
      open(hdl, file=trim(cfile_prefix)//'.bpcoef', status='unknown', action='write', position='append', iostat=istat)
      if (istat /= 0) then
         print '(2a)', 'Error: failed to open bpcoef file. ', trim(cfile_prefix)//'.bpcoef'
         flush(output_unit)
         error stop
      endif

      if (temp_independent == 0) then
         write(hdl, '(a, f8.3)') '########## tempK = ', tempK
      else
         write(hdl, '(a)') '########## temperature independent'
      endif

      bp_map_dG(:,:) = 0.0_PREC

      do imp = 1, nmp-1
         i = lmp_mp(imp)
         ichain = ichain_mp(imp)

         ! Either 5' or 3' end
         if (i == 1 .or. i == nmp_chain(ichain)) cycle

         do jmp = imp+1, nmp

            if (bp_map_0(imp, jmp) == 0) cycle

            j = lmp_mp(jmp)
            jchain = ichain_mp(jmp)

            dH = 0.0_PREC
            dS = 0.0_PREC

            if (bp_map_0(imp-1, jmp+1) > 0) then
               dH = 0.5_PREC * (NN_dH(seqt2nnt(seq(i-1, ichain), seq(i, ichain), seq(j+1, jchain), seq(j, jchain))) - dH0)
               dS = 0.5_PREC * (NN_dS(seqt2nnt(seq(i-1, ichain), seq(i, ichain), seq(j+1, jchain), seq(j, jchain))) - dS0)
            endif

            if (bp_map_0(imp+1, jmp-1) > 0) then
               dH = dH + 0.5_PREC * (NN_dH(seqt2nnt(seq(i, ichain), seq(i+1, ichain), seq(j, jchain), seq(j-1, jchain))) - dH0) 
               dS = dS + 0.5_PREC * (NN_dS(seqt2nnt(seq(i, ichain), seq(i+1, ichain), seq(j, jchain), seq(j-1, jchain))) - dS0)
            endif


            ! Default
            if (temp_independent == 0) then
               dG = coef_dG * (dH - tempK * 1.0e-3_PREC * dS)
            else
               dG = coef_dG * dH
            endif


            if (dG < 0.0_PREC) then
               bp_map_dG(imp, jmp) = dG
               bp_map_dG(jmp, imp) = dG
               write(hdl, '(i5,1x,i5,3x,7a1,3x,f8.3)') imp, jmp, &
                         seqt2char(seq(i-1,ichain)), seqt2char(seq(i,ichain)), seqt2char(seq(i+1,ichain)), '/', &
                         seqt2char(seq(j+1,jchain)), seqt2char(seq(j,jchain)), seqt2char(seq(j-1,jchain)), dG

            else
               bp_map(imp, jmp) = 0
               bp_map(jmp, imp) = 0

            endif
         enddo
      enddo

      close(hdl)
      iopen_hdl = iopen_hdl - 1

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Specific pairs given in CT file
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   else if (bp_model == 2) then

      continue

   endif

endsubroutine set_bp_map
