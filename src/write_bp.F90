subroutine write_bp(irep, tempK_in)

   use ieee_exceptions, only : IEEE_GET_HALTING_MODE, IEEE_SET_HALTING_MODE, IEEE_UNDERFLOW
   use const, only : PREC
   use const_phys, only : BOLTZ_KCAL_MOL
   use var_state, only : bp_status, ene_bp
   use var_potential, only : bp_mp, nbp
   use var_io, only : flg_out_bp, flg_out_bpall, flg_out_bpe, hdl_bp, hdl_bpall, hdl_bpe, KIND_OUT_BP, KIND_OUT_BPE

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(in) :: tempK_in

   integer :: ibp, imp, jmp
   real(PREC) :: kT
   logical :: halt_mode

   if (.not. (flg_out_bp .or. flg_out_bpall .or. flg_out_bpe)) then
      return
   endif

   kT = tempK_in * BOLTZ_KCAL_MOL

   if (flg_out_bp) then

      call ieee_get_halting_mode(IEEE_UNDERFLOW, halt_mode)
      call ieee_set_halting_mode(IEEE_UNDERFLOW, halting=.false. )

      do ibp = 1, nbp(irep)

         if (bp_status(ibp, irep)) then
            if (ene_bp(ibp, irep) <= -kT) then
               imp = bp_mp(1, ibp, irep)
               jmp = bp_mp(2, ibp, irep)
               write(hdl_bp(irep)) int(imp,kind=KIND_OUT_BP), int(jmp,kind=KIND_OUT_BP), &
                                   real(ene_bp(ibp, irep), kind=KIND_OUT_BPE)
            endif
         endif
      enddo

      write(hdl_bp(irep)) int(0,kind=KIND_OUT_BP), int(0,kind=KIND_OUT_BP), &
                          real(0.0, kind=KIND_OUT_BPE)

      call ieee_set_halting_mode(IEEE_UNDERFLOW, halting=halt_mode)
   endif

   if (flg_out_bpall) then

      call ieee_get_halting_mode(IEEE_UNDERFLOW, halt_mode)
      call ieee_set_halting_mode(IEEE_UNDERFLOW, halting=.false. )

      do ibp = 1, nbp(irep)

         if (bp_status(ibp, irep)) then
            imp = bp_mp(1, ibp, irep)
            jmp = bp_mp(2, ibp, irep)
            write(hdl_bpall(irep)) int(imp,kind=KIND_OUT_BP), int(jmp,kind=KIND_OUT_BP), &
                                   real(ene_bp(ibp, irep), kind=KIND_OUT_BPE)
         endif
      enddo

      write(hdl_bpall(irep)) int(0,kind=KIND_OUT_BP), int(0,kind=KIND_OUT_BP), &
                             real(0.0, kind=KIND_OUT_BPE)

      call ieee_set_halting_mode(IEEE_UNDERFLOW, halting=halt_mode)
   endif

   if (flg_out_bpe) then

      do ibp = 1, nbp(irep)

         if (bp_status(ibp, irep)) then
            if (ene_bp(ibp, irep) <= -kT) then
               imp = bp_mp(1, ibp, irep)
               jmp = bp_mp(2, ibp, irep)
               write(hdl_bpe(irep), '(1x,i5,1x,i5,1x,f6.2)', advance='no') imp, jmp, ene_bp(ibp, irep)
            endif
         endif
      enddo

      write(hdl_bpe(irep), '(a)') ''
   endif

endsubroutine write_bp
