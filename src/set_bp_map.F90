subroutine set_bp_map()

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const
   use const_idx, only : SEQT, seqt2char, seqt2nnt
   use var_top, only : nmp, seq, lmp_mp, ichain_mp, nmp_chain
   use var_potential, only : bp_model, bp_map, bp_map_0, bp_min_loop, bp_map_dG, &
                             NN_dG, NN_dH, NN_dS, dH0, dS0, coef_dG
   use var_state, only : tempK, temp_independent

   implicit none

   integer :: imp, jmp
   integer :: i, j, ichain, jchain
   integer :: l, n, idummy
   integer :: istat, hdl
   real(PREC) :: dG, dH, dS
   character(len=1) :: nt

   bp_map(:,:) = bp_map_0(:,:)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! All pairwise
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (bp_model == 5) then

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
            if (i == 2 .or. j+1 == nmp_chain(jchain)) then
               continue
            else if (is_complement(seq(i-1, ichain), seq(j+1, jchain))) then
               dH = 0.5 * (NN_dH(seqt2nnt(seq(i-1, ichain), seq(i, ichain), seq(j+1, jchain), seq(j, jchain))) - dH0)
               dS = 0.5 * (NN_dS(seqt2nnt(seq(i-1, ichain), seq(i, ichain), seq(j+1, jchain), seq(j, jchain))) - dS0)
            endif

            if (i+1 == nmp_chain(ichain) .or. j == 2) then
               continue
            else if (is_complement(seq(i+1, ichain), seq(j-1, jchain))) then
               dH = dH + 0.5 * (NN_dH(seqt2nnt(seq(i, ichain), seq(i+1, ichain), seq(j, jchain), seq(j-1, jchain))) - dH0) 
               dS = dS + 0.5 * (NN_dS(seqt2nnt(seq(i, ichain), seq(i+1, ichain), seq(j, jchain), seq(j-1, jchain))) - dS0)
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
               !print '(i5,1x,i5,3x,7a1,3x,f8.3)', imp, jmp, &
               !          seqt2char(seq(i-1,ichain)), seqt2char(seq(i,ichain)), seqt2char(seq(i+1,ichain)), '/', &
               !          seqt2char(seq(j+1,jchain)), seqt2char(seq(j,jchain)), seqt2char(seq(j-1,jchain)), dG

            else
               bp_map(imp, jmp) = 0
               bp_map(jmp, imp) = 0

            endif
         enddo
      enddo

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Specific pairs given in CT file
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   else if (bp_model == 2) then

      continue

   endif

contains

   logical function is_complement(s1, s2)

      integer, intent(in) :: s1, s2

      is_complement = .False.

      if (s1 == SEQT%A) then
         if (s2 == SEQT%U) then
            is_complement = .True.
         endif

      else if (s1 == SEQT%U) then
         if (s2 == SEQT%A .or. s2 == SEQT%G) then
            is_complement = .True.
         endif

      else if (s1 == SEQT%G) then
         if (s2 == SEQT%C .or. s2 == SEQT%U) then
            is_complement = .True.
         endif

      else if (s1 == SEQT%C) then
         if (s2 == SEQT%G) then
            is_complement = .True.
         endif

      endif

   end function is_complement

endsubroutine set_bp_map
