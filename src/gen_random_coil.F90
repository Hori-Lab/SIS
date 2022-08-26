subroutine gen_random_coil(nmp, xyz, origin)

   use mt19937_64
   use const, only : PREC
   use const_phys, only : PI

   implicit none

   integer, intent(in) :: nmp
   real(PREC), intent(out) :: xyz(3, nmp)
   real(PREC), intent(out) :: origin(3)

   integer :: imp, jmp
   real(PREC) :: trial(3), v(3)
   real(PREC) :: dih, d2
   logical :: flg

   integer, parameter :: Nback = 10           ! Number of particles removed backward when new bead causes clash.
   real(PREC), parameter :: BL = 5.84_PREC    ! Bond length
   real(PREC), parameter :: BA = 2.643_PREC   ! Bond angle
   real(PREC), parameter :: EV = 12.0_PREC    ! Excluded volume distance
   real(PREC), parameter :: EV2 = EV * EV


   ! Coordinates for the first three beads.
   xyz(:, 1) = origin(:)
   xyz(:, 2) = origin(:) + (/BL, 0.0_PREC, 0.0_PREC/)
   xyz(:, 3) = origin(:) + (/BL*(1.0_PREC + cos(PI - BA)), BL*sin(PI - BA), 0.0_PREC/)

   imp = 4

   do while (imp <= nmp)
      print *, imp
      flush(6)
      ! Generate a dihedral angle at random
      dih = -PI + genrand64_real1() * 2 * PI

      ! Derive a new coordinate by NeRF
      call NeRF(xyz(:, imp-3), xyz(:, imp-2), xyz(:, imp-1), BL, BA, dih, trial)

      ! Check if any clash
      flg = .False.
      do jmp = 1, imp-4

         v(:) = trial(:) - xyz(:, jmp)
         d2 = dot_product(v, v)

         if (d2 < EV2) then
             flg = .True.
             exit
         endif
      enddo

      ! If clash, remove the last Nback beads
      if (flg) then
         imp = max(4, imp - 10)
         cycle
      endif

      ! Add the new bead
      xyz(:, imp) = trial(:)

      imp = imp + 1
   enddo

contains

   ! Natural Extension of Reference Frame
   subroutine NeRF(a, b, c, bond, angl, dihd, d)

      real(PREC), intent(in) :: a(3), b(3), c(3)
      real(PREC), intent(in) :: bond, angl, dihd
      real(PREC), intent(out) :: d(3)

      real(PREC) :: t
      real(PREC) :: d2(3), ab(3), bc(3), n(3), m(3,3), n_bc(3)

      t = PI - angl

      d2(1) = bond * cos(t)
      d2(2) = bond * cos(dihd) * sin(t)
      d2(3) = bond * sin(dihd) * sin(t)

      ab(:) = b(:) - a(:)
      bc(:) = c(:) - b(:)

      n_bc(:) = bc(:) / norm2(bc)

      n(:) = cross(ab, n_bc)
      n(:) = n(:) / norm2(n)

      m(:, 1) = n_bc(:)
      m(:, 2) = cross(n ,n_bc)
      m(:, 3) = n(:)
      d(:) = matmul(m, d2) + c(:)

   endsubroutine NeRF

   function cross(a, b)
      real(PREC), intent(in) :: a(3), b(3)
      real(PREC) :: cross(3)

      cross(1) = a(2) * b(3) - a(3) * b(2)
      cross(2) = a(3) * b(1) - a(1) * b(3)
      cross(3) = a(1) * b(2) - a(2) * b(1)
   endfunction cross

endsubroutine gen_random_coil
