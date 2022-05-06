program main

   implicit none

   integer, parameter :: PREC = 8
   real(PREC), parameter :: PI = 4.D0 * atan(1.0D0)

   ! parameters
   real(PREC) :: bp_bond_k
   real(PREC) :: bp_bond_r
   real(PREC) :: bp_angl_k
   real(PREC) :: bp_angl_theta1
   real(PREC) :: bp_angl_theta2
   real(PREC) :: bp_dihd_k
   real(PREC) :: bp_dihd_phi1
   real(PREC) :: bp_dihd_phi2
   real(PREC) :: bp_U0

   ! variables
   real(PREC) :: d, t1_1, t1_2, t2_1, t2_2, p1, p2

   integer :: i
   integer, parameter :: fh = 10 ! file handle
   real(PREC) :: degree

   !!!!!! NHT 22 original
   open(fh, file='bp_energy_NHT22.out', status='unknown')

   call reset_parameters()

   call reset_variables(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   write(fh, '(a)') '# distance'
   do i = 1, 400
      d = 0.05 * i
      write(fh,*) d, bp_energy(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   enddo
   write(fh,*)
   write(fh,*)

   call reset_variables(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   write(fh, '(a)') '# angle'
   do i = 1, 180
      degree = real(i, kind=PREC)
      t1_1 = degree / 180.0_PREC * PI
      write(fh,*) degree, bp_energy(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   enddo
   write(fh,*)
   write(fh,*)

   call reset_variables(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   write(fh, '(a)') '# dihedral'
   do i = -180, 180
      degree = real(i, kind=PREC)
      p1 = degree / 180.0_PREC * PI
      write(fh,*) degree, bp_energy(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   enddo
   write(fh,*)
   write(fh,*)

   !!!!!! para3, refined angle and dihedral, k_theta = 3.2, k_phi = 1.3
   open(fh, file='bp_energy_para3.out', status='unknown')

   call reset_parameters()
   bp_angl_k = 3.2_PREC
   bp_dihd_k = 1.3_PREC
   bp_U0 = -8.40_PREC

   call reset_variables(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   write(fh, '(a)') '# distance'
   do i = 1, 400
      d = 0.05 * i
      write(fh,*) d, bp_energy(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   enddo
   write(fh,*)
   write(fh,*)

   call reset_variables(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   write(fh, '(a)') '# angle'
   do i = 1, 180
      degree = real(i, kind=PREC)
      t1_1 = degree / 180.0_PREC * PI
      write(fh,*) degree, bp_energy(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   enddo
   write(fh,*)
   write(fh,*)

   call reset_variables(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   write(fh, '(a)') '# dihedral'
   do i = -180, 180
      degree = real(i, kind=PREC)
      p1 = degree / 180.0_PREC * PI
      write(fh,*) degree, bp_energy(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   enddo
   write(fh,*)
   write(fh,*)

   close(fh)


   stop

contains

subroutine reset_parameters()
   ! bond
   bp_bond_k =  3.0_PREC     ! A^(-2)
   bp_bond_r = 13.8_PREC     ! A

   ! angle
   bp_angl_k =  1.5_PREC          ! rad^(-2)
   bp_angl_theta1 = 1.8326_PREC   ! rad
   bp_angl_theta2 = 0.9425_PREC   ! rad

   ! dihedral
   bp_dihd_k =  0.5_PREC        ! --
   bp_dihd_phi1 = 1.8326_PREC   ! rad
   bp_dihd_phi2 = 1.1345_PREC   ! rad

   ! Energy scale
   bp_U0 = -5.0_PREC  ! kcal/mol
end subroutine reset_parameters

subroutine reset_variables(d, t1_1, t1_2, t2_1, t2_2, p1, p2)
   real(PREC), intent(out) :: d, t1_1, t1_2, t2_1, t2_2, p1, p2
   d = 13.8_PREC     ! A
   t1_1 = 1.8326_PREC   ! rad
   t1_2 = 1.8326_PREC   ! rad
   t2_1 = 0.9425_PREC   ! rad
   t2_2 = 0.9425_PREC   ! rad
   p1 = PI - 1.8326_PREC   ! rad
   p2 = PI - 1.1345_PREC   ! rad
end subroutine reset_variables

pure function bp_energy(d, t1_1, t1_2, t2_1, t2_2, p1, p2) result (Ebp)

   integer, parameter :: PREC = 8

   real(PREC) :: Ebp  ! Result energy

   real(PREC), intent(in) :: d    ! distance
   real(PREC), intent(in) :: t1_1 ! theta i, j, j-1
   real(PREC), intent(in) :: t1_2 ! theta i-1, i, j
   real(PREC), intent(in) :: t2_1 ! theta i, j, j+1
   real(PREC), intent(in) :: t2_2 ! theta i+1, i, j
   real(PREC), intent(in) :: p1   ! phi j-1, j, i, i-1
   real(PREC), intent(in) :: p2   ! phi j+1, j, i, i+1
   
   real(PREC) :: u

   u = bp_bond_k * (d - bp_bond_r)**2

   u = u + bp_angl_k * (t1_1 - bp_angl_theta1)**2

   u = u + bp_angl_k * (t1_2 - bp_angl_theta1)**2

   u = u + bp_angl_k * (t2_1 - bp_angl_theta2)**2

   u = u + bp_angl_k * (t2_2 - bp_angl_theta2)**2

   u = u + bp_dihd_k * (1.0 + cos(p1 + bp_dihd_phi1))

   u = u + bp_dihd_k * (1.0 + cos(p2 + bp_dihd_phi2))

   Ebp = bp_U0 * exp(-u)

end function bp_energy

end program
