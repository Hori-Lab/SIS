subroutine force_ele_DH(forces)

  use const, only : PREC
  use pbc, only : pbc_vec_d
  use var_top, only : nmp
  use var_state, only : xyz, lambdaD
  use var_potential, only : ele_coef, ele_cutoff, nele, ele_mp

  implicit none

  real(PREC), intent(inout) :: forces(3, nmp)

  integer :: iele, imp1, imp2
  real(PREC) :: dist1, dist2, rdist1
  real(PREC) :: dvdw_dr, rcdist, cutoff2
  real(PREC) :: v21(3), for(3)
  !character(ARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------

  cutoff2 = ele_cutoff ** 2
  rcdist = 1.0_PREC / lambdaD

  !$omp do private(imp1,imp2,v21,dist2,dist1,rdist1,dvdw_dr,for)
  do iele = 1, nele
     imp1 = ele_mp(1, iele)
     imp2 = ele_mp(2, iele)

     v21(:) = pbc_vec_d(xyz(:,imp2), xyz(:,imp1))

     dist2 = dot_product(v21,v21)

     if(dist2 > cutoff2) cycle

     ! -----------------------------------------------------------------
     dist1 = sqrt(dist2)
     rdist1 = 1.0 / dist1

     dvdw_dr = ele_coef * rdist1 * rdist1 &
              * (rdist1 + rcdist) * exp(-dist1 * rcdist)

     
     !if(dvdw_dr > DE_MAX) then
     !   write(error_message,*) 'force_ele_DH > DE_MAX', istep, grep, imp1, imp2, dist1, dvdw_dr, DE_MAX
     !   write(6, '(a)') trim(error_message)
     !   dvdw_dr = DE_MAX
     !end if
        
     for(1:3) = dvdw_dr * v21(1:3)
     forces(1:3, imp1) = forces(1:3, imp1) - for(1:3)
     forces(1:3, imp2) = forces(1:3, imp2) + for(1:3)
  end do
  !$omp end do nowait

end subroutine force_ele_DH
