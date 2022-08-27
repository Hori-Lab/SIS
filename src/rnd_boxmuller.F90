function rnd_boxmuller(irep)

   use const, only : PREC
   !use mt19937_64
   use mt_stream
   use var_state, only : mts

   implicit none

   integer, intent(in) :: irep
   real(PREC) :: rnd_boxmuller

   logical, save :: flg_first = .true.
   logical, save :: flg_stored
   real(PREC), save :: stored_value
   real(PREC) :: vx, vy, r2, rf
  
   ! --------------------------------------------------------------------
   if (flg_first) then
      flg_first = .false.
      flg_stored = .false.
      stored_value = 0.0
   endif


   if (.not. flg_stored) then
      flg_stored = .true.

      do
         !vx = 2.0e0_PREC * genrand64_real1() - 1.0e0_PREC
         !vy = 2.0e0_PREC * genrand64_real1() - 1.0e0_PREC
         vx = 2.0e0_PREC * genrand_double1(mts(irep)) - 1.0e0_PREC
         vy = 2.0e0_PREC * genrand_double1(mts(irep)) - 1.0e0_PREC

         r2 = vx*vx + vy*vy
         if(r2 < 1.0e0_PREC .and. r2 > 0.0e0_PREC) exit
      enddo

      rf = sqrt(-2.0e0_PREC * log(r2) / r2)
      stored_value = vx * rf
      rnd_boxmuller = vy * rf

   else
      flg_stored = .False.
      rnd_boxmuller = stored_value

   endif
  
endfunction rnd_boxmuller
