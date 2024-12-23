subroutine write_bpcoef()

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const
   use const_idx, only : seqt2char, REPT, MOLT
   use var_io, only : flg_out_bpcoef, hdl_bpcoef
   use var_top, only : nmp, seq, lmp_mp, ichain_mp, moltypes, nmp_chain
   use var_replica, only : flg_replica, irep2grep, rep2val, nrep_proc
   use var_potential, only : bp_model, coef_dG, bp3_map, bp3_dH, bp3_dS
   use var_state, only : tempK, temp_independent

   implicit none

   integer :: imp, jmp, irep, grep
   integer :: i, j, ichain, jchain
   integer :: i_pre, i_nxt, j_pre, j_nxt
   real(PREC) :: dG, dH, dS, tK

   if (.not. flg_out_bpcoef) return

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! All pairwise
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (bp_model == 5) then

      !if (temp_independent == 0) then
      !   write(hdl, '(a, f8.3)') '########## tempK = ', tempK
      !else
      !   write(hdl, '(a)') '########## temperature independent'
      !endif

      do imp = 1, nmp-1
         i = lmp_mp(imp)
         ichain = ichain_mp(imp)

         i_pre = i - 1
         i_nxt = i + 1
         if (moltypes(ichain) == MOLT%CIRCRNA) then
            if (i == 1) i_pre = nmp_chain(ichain)
            if (i == nmp_chain(ichain)) i_nxt = 1
         endif

         do jmp = imp+1, nmp
            if (bp3_map(imp, jmp) == 0) cycle

            j = lmp_mp(jmp)
            jchain = ichain_mp(jmp)

            j_pre = j - 1
            j_nxt = j + 1
            if (moltypes(jchain) == MOLT%CIRCRNA) then
               if (j == 1) j_pre = nmp_chain(jchain)
               if (j == nmp_chain(jchain)) j_nxt = 1
            endif

            dH = bp3_dH(bp3_map(imp, jmp))  ! dH
            dS = bp3_dS(bp3_map(imp, jmp))  ! dS (0.001 already multiplied so that the unit is kcal/mol/K)

            do irep = 1, nrep_proc

               if (flg_replica) then
                  grep = irep2grep(irep)
                  tK = rep2val(grep, REPT%TEMP)
               else
                  grep = 1
                  tK = tempK
               endif

               if (temp_independent == 0) then
                  dG = coef_dG * (dH - tK * dS)
               else
                  dG = coef_dG * dH
               endif

               if (dG < 0.0_PREC) then
                  write(hdl_bpcoef(irep), '(i5,1xi5,1x,i5,3x,7a1,3x,f8.3)') grep, imp, jmp, &
                     seqt2char(seq(i_pre,ichain)), seqt2char(seq(i,ichain)), seqt2char(seq(i_nxt,ichain)), '/', &
                     seqt2char(seq(j_nxt,jchain)), seqt2char(seq(j,jchain)), seqt2char(seq(j_pre,jchain)), dG

               endif

            enddo
         enddo
      enddo

   endif

   do irep = 1, nrep_proc
      flush(hdl_bpcoef(irep))
   enddo
endsubroutine write_bpcoef
