subroutine energy_bp(Ebp)

   use const
   use const_idx, only : SEQT
   use var_state, only : xyz
   use var_top, only : nmp_chain, seq, imp_chain, nchains
   use var_potential

   implicit none
  
   real(PREC), intent(inout) :: Ebp

   integer :: ichain, jchain
   integer :: i, j, j_start, imp, jmp
   real(PREC) :: nhb, u
   real(PREC) :: d, theta, phi

   do ichain = 1, nchains
      do jchain = ichain, nchains

         do i = 2, nmp_chain(ichain)-1

            imp = imp_chain(i, ichain)

            if (ichain == jchain) then
               j_start = i + 5
            else
               j_start = 2
            endif

            do j = j_start, nmp_chain(jchain)-1
               
               jmp = imp_chain(j, jchain)

               if ((seq(i,ichain) == SEQT%G .and. seq(j, jchain) == SEQT%C) .or. &
                   (seq(i,ichain) == SEQT%C .and. seq(j, jchain) == SEQT%G) ) then
                  nhb = 3

               else if ((seq(i,ichain) == SEQT%A .and. seq(j, jchain) == SEQT%U) .or. &
                        (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%A) ) then
                  nhb = 2

               else
                  cycle

               endif
                  
               d = mp_distance(imp, jmp)

               if (d >= 18.0) cycle
               
               u = Ubp_bond_k * (d - Ubp_bond_r)**2

               theta = mp_angle(imp, jmp, jmp-1)
               u = u + Ubp_angl_k * (theta - Ubp_angl_theta1)**2

               theta = mp_angle(imp-1, imp, jmp)
               u = u + Ubp_angl_k * (theta - Ubp_angl_theta1)**2

               theta = mp_angle(imp, jmp, jmp+1)
               u = u + Ubp_angl_k * (theta - Ubp_angl_theta2)**2

               theta = mp_angle(imp+1, imp, jmp)
               u = u + Ubp_angl_k * (theta - Ubp_angl_theta2)**2

               phi = mp_dihedral(imp-1, imp, jmp, jmp-1)
               u = u + Ubp_dihd_k * (1.0 + cos(phi + Ubp_dihd_phi1))

               phi = mp_dihedral(imp+1, imp, jmp, jmp+1)
               u = u + Ubp_dihd_k * (1.0 + cos(phi + Ubp_dihd_phi2))

               Ebp = Ebp + nhb * Ubp0 * exp(-u)

            enddo
         enddo

      enddo
   enddo

contains

   function pbc_vec(v) result (new_vec)
      
      use var_top, only : pbc_box, pbc_box_half

      real(PREC) :: new_vec(3)
      real(PREC), intent(in) :: v(3)

      integer :: i

      do i = 1, 3
         if(v(i) > pbc_box_half(i)) then
            new_vec(i) = v(i) - pbc_box(i)

         else if(v(i) < -pbc_box_half(i)) then
            new_vec(i) = v(i) + pbc_box(i)

         else
            new_vec(i) = v(i)

         end if

      end do
    
   end function pbc_vec

   function mp_distance(imp1, imp2) result(d)

      real(PREC) :: d
      integer, intent(in) :: imp1, imp2

      d = norm2( pbc_vec(xyz(:,imp1) - xyz(:, imp2)) )

   endfunction mp_distance

   function mp_angle(imp1, imp2, imp3) result(theta)

      real(PREC) :: theta
      integer, intent(in) :: imp1, imp2, imp3
      real(PREC) :: v12(3), v32(3)
      real(PREC) :: co

      v12(:) = pbc_vec(xyz(:, imp1) - xyz(:, imp2))
      v32(:) = pbc_vec(xyz(:, imp3) - xyz(:, imp2))

      co = dot_product(v32, v12) / sqrt(dot_product(v12,v12) * dot_product(v32,v32))

      if(co > 1.0e0_PREC) then
         co = 1.0e0_PREC
      else if(co < -1.0e0_PREC) then
         co = -1.0e0_PREC
      end if

      theta = acos(co)

   endfunction mp_angle

   function mp_dihedral(imp1, imp2, imp3, imp4) result(phi)

      real(PREC) :: phi
      integer, intent(in) :: imp1, imp2, imp3, imp4
      real(PREC) :: v12(3), v32(3), v34(3), m(3), n(3)

      v12(:) =  pbc_vec(xyz(:, imp1) - xyz(:, imp2))
      v32(:) =  pbc_vec(xyz(:, imp3) - xyz(:, imp2))
      v34(:) =  pbc_vec(xyz(:, imp3) - xyz(:, imp4))

      m(1) = v12(2)*v32(3) - v12(3)*v32(2)
      m(2) = v12(3)*v32(1) - v12(1)*v32(3)
      m(3) = v12(1)*v32(2) - v12(2)*v32(1)
      n(1) = v32(2)*v34(3) - v32(3)*v34(2)
      n(2) = v32(3)*v34(1) - v32(1)*v34(3)
      n(3) = v32(1)*v34(2) - v32(2)*v34(1)

      phi = atan2(dot_product(v12,n)*norm2(v32), dot_product(m,n))

   endfunction mp_dihedral

end subroutine energy_bp
