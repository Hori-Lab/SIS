subroutine set_ele()

   use const
   use const_phys, only : PI, EPS0, BOLTZ_J, N_AVO, ELE, JOUL2KCAL_MOL
   use var_state, only : ionic_strength, lambdaD, diele, length_per_charge, &
                         tempK, temp_independent, diele_dTcoef, temp_ref
   use var_potential, only : ele_coef, ele_cutoff_type, ele_cutoff_inp, ele_cutoff
   use var_top, only : nmp, inp_no_charge, has_charge, charge

   implicit none

   integer :: imp, i
   real(PREC) :: Tc, lb, Zp
   real(PREC), parameter ::  MM_A=87.740e0_PREC, MM_B=-0.4008e0_PREC
   real(PREC), parameter ::  MM_C=9.398e-4_PREC, MM_D=-1.410e-6_PREC

   logical, save :: flg_first = .True.

#ifdef _HTN_CONSISTENT
   real(PREC) :: temp_kT, lb_kT, kappaD, rho
#endif

   ! -----------------------------------------------------------------------
   ! Dielectric constant
   ! Temperature dependent (Malmberg and Maryott, 1956)

#ifdef _HTN_CONSISTENT
   !To be consistent with Hung's code
   temp_kT = tempk * BOLTZ_J * JOUL2KCAL_MOL
   diele = 296.0736276 - 619.2813716 * temp_kT + 531.2826741 * temp_kT**2 - 180.0369914 * temp_kT**3

   ! Bjerrum length
   lb_kT = 332.0637 / diele
   lb = lb_kT / temp_kT
   Zp = length_per_charge / lb

   if (temp_independent /== 0) then
      print '(a)', 'Error: temp_independnet /= 0 not implemented in _HTN_CONSISTENT (set_ele.F90)'
      flush(6)
      error stop
   endif

#else
   ! Default
   if (temp_independent == 0) then
      Tc = tempk - 273.15_PREC
      diele =  MM_A + MM_B*Tc + MM_C*Tc*Tc + MM_D*Tc*Tc*Tc

   else
      Tc = temp_ref - 273.15_PREC
      diele =  MM_A + MM_B*Tc + MM_C*Tc*Tc + MM_D*Tc*Tc*Tc
      diele_dTcoef = 1.0_PREC + temp_ref / diele &
                    * (MM_B + 2.0_PREC*MM_C*Tc + 3.0_PREC*MM_D*Tc*Tc)
   endif

   ! Bjerrum length
   lb = ELE * ELE / (4.0e0_PREC * PI * EPS0 * diele * BOLTZ_J * tempk) * 1.0e10_PREC

   ! Calculate phosphate charge taking into account monovalent salt condensation
   !xi = lb / length_per_charge
   !theta = 1.0 - 1.0/xi
   Zp = length_per_charge / lb
#endif

   allocate(has_charge(nmp))
   allocate(charge(nmp))

   has_charge(:) = .True.
   if (allocated(inp_no_charge)) then
      do i = 1, size(inp_no_charge)
         has_charge(inp_no_charge(i)) = .False.
      enddo
   endif

   charge(:) = 0.0_PREC

   ! Reflect Zp
   do imp = 1, nmp
      if (.not. has_charge(imp)) cycle
      charge(imp) = -Zp
   enddo

   ! ----------------------------------------------------------------------
   ! coef: j_kcal * eq**2 / (4.0e0_PREC * PI * e0 * ek * rij)
   !   =  332.063713019 / ek
#ifdef _HTN_CONSISTENT
   ele_coef = lb_kT
#else
   ele_coef = JOUL2KCAL_MOL * 1.0e10_PREC * ELE**2 / (4.0e0_PREC * PI * EPS0 * diele)
#endif

   ! Currently, all charges are phosphate with the same value of charge (-Zp).
   ! Therefore it can be included in coef.
   ele_coef = ele_coef * Zp**2

   ! Kd: sqrt(e0 * ek * RT / 2 * NA**2 * eq**2 * I)
   !   = sqrt(1.0e-3 * e0 * ek * kb / 2 * NA * eq**2) * sqrt(T(K) / I(M))

#ifdef _HTN_CONSISTENT
   rho = 2 * ionic_strength * 6.022e-4_PREC
   kappaD = sqrt(4 * 3.14159 * lb_kT * rho / temp_kT)
   lambdaD = 1.0_PREC / kappaD
#else
   lambdaD = 1.0e10_PREC * sqrt( (1.0e-3_PREC * EPS0 * diele * BOLTZ_J) &
                                  / (2.0_PREC * N_AVO * ELE**2)  )     &
                                * sqrt(tempk / ionic_strength)
#endif

   ! Set cutoff
   if (ele_cutoff_type == 1) then
      ele_cutoff = ele_cutoff_inp
   else if (ele_cutoff_type == 2) then
      ele_cutoff = ele_cutoff_inp * lambdaD
   else
      error stop "Error: Invalid ele_cutoff_type value in set_ele."
   endif

   if (flg_first) then
      write(6, '(a)') 'Set electrostaic parameters'
      write(6, '(a,g15.8)') '# Dielectric constant (H2O): ', diele
      write(6, '(a,g15.8)') '# coef: ', ele_coef / (Zp**2)
      write(6, '(a,g15.8)') '# Bjerrum length: ', lb
      write(6, '(a,g15.8)') '# Debye length (lambdaD): ', lambdaD
      write(6, '(a,g15.8)') '# Reduced charge on each phosphate: ', -Zp
      write(6, *)

      flg_first = .False.
   endif

end subroutine set_ele
