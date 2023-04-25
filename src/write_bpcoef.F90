subroutine write_bpcoef()

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const
   use const_idx, only : seqt2char, REPT
   use var_io, only : flg_out_bpcoef, hdl_bpcoef
   use var_top, only : nmp, seq, lmp_mp, ichain_mp
   use var_replica, only : flg_replica, irep2grep, rep2val, nrep_proc
   use var_potential, only : bp_model, coef_dG, bp3_map, bp3_dH, bp3_dS
   use var_state, only : tempK, temp_independent

   implicit none

   integer :: imp, jmp, irep, grep
   integer :: i, j, ichain, jchain
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

         do jmp = imp+1, nmp
            if (bp3_map(imp, jmp) == 0) cycle

            j = lmp_mp(jmp)
            jchain = ichain_mp(jmp)

            dH = bp3_dH(bp3_map(imp, jmp))  ! dH
            dS = bp3_dS(bp3_map(imp, jmp))  ! dS (0.001 already multiplied so that the unit is kcal/mol/K)

            do irep = 1, nrep_proc

               if (flg_replica) then
                  grep = irep2grep(irep)
                  tK = rep2val(grep, REPT%TEMP)
               else
                  tK = tempK
               endif

               if (temp_independent == 0) then
                  dG = coef_dG * (dH - tK * dS)
               else
                  dG = coef_dG * dH
               endif

               if (dG < 0.0_PREC) then
                  write(hdl_bpcoef(irep), '(i5,1xi5,1x,i5,3x,7a1,3x,f8.3)') grep, imp, jmp, &
                            seqt2char(seq(i-1,ichain)), seqt2char(seq(i,ichain)), seqt2char(seq(i+1,ichain)), '/', &
                            seqt2char(seq(j+1,jchain)), seqt2char(seq(j,jchain)), seqt2char(seq(j-1,jchain)), dG

               endif

            enddo
         enddo
      enddo

   endif

endsubroutine write_bpcoef
