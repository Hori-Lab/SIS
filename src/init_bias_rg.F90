subroutine init_bias_rg()
   
   use const_idx, only : REPT
   use var_potential, only : bias_rg_0_inp, bias_rg_0
   use var_replica, only : flg_repvar, nrep_proc, irep2grep, rep2val

   implicit none

   integer :: irep, grep

   allocate(bias_rg_0(nrep_proc))

   if (flg_repvar(REPT%RG)) then
      do irep = 1, nrep_proc
         grep = irep2grep(irep)
         bias_rg_0(irep) = rep2val(grep, REPT%RG)
      enddo

   else
      bias_rg_0(:) = bias_rg_0_inp

   endif

endsubroutine init_bias_rg
