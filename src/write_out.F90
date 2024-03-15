subroutine write_out_header(irep)

   use var_io, only : hdl_out
   use var_replica, only : flg_replica
   use var_potential, only : flg_stage, flg_ele, flg_twz

   integer, intent(in) :: irep

   integer :: icol

   ! .out file header
   if (flg_replica) then
      write(hdl_out(irep), '(a)', advance='no') '#(1)nframe   (2)R (3)T   (4)Ekin       (5)Epot       (6)Ebond     '
                                                !123456789012 1234 123456 1234567890123 1234567890123 1234567890123'
      write(hdl_out(irep), '(a)', advance='no') ' (7)Eangl      (8)Edih       (9)Ebp        (10)Eexv     '
                                                ! 1234567890123 1234567890123 1234567890123 1234567890123'
      icol = 10
   else
      write(hdl_out(irep), '(a)', advance='no') '#(1)nframe   (2)T   (3)Ekin       (4)Epot       (5)Ebond     '
                                                !123456789012 123456 1234567890123 1234567890123 1234567890123'
      write(hdl_out(irep), '(a)', advance='no') ' (6)Eangl      (7)Edih       (8)Ebp        (9)Eexv      '
                                                ! 1234567890123 1234567890123 1234567890123 1234567890123'
      icol = 9
   endif

   if (flg_ele) then
      icol = icol + 1
                                                     ! 1   23     4567890123
      write(hdl_out(irep), '(a,i2,a)', advance='no') ' (', icol, ')Eele     '
   endif

   if (flg_stage) then
      icol = icol + 1
      write(hdl_out(irep), '(a,i2,a)', advance='no') ' (', icol, ')Estage   '
   endif

   if (flg_twz) then
      icol = icol + 1
      write(hdl_out(irep), '(a,i2,a)', advance='no') ' (', icol, ')Etweezers'
   endif

   write(hdl_out(irep), '(a)') ''
endsubroutine write_out_header

subroutine write_out(irep, istep_out, rep_label, tK)

   use const, only : PREC
   use const_idx, only : ENE
   use var_io, only : hdl_out
   use var_replica, only : flg_replica
   use var_potential, only : flg_stage, flg_ele, flg_twz
   use var_state, only : Ekinetic, energies

   integer, intent(in) :: irep
   integer, intent(in) :: istep_out
   integer, intent(in) :: rep_label
   real(PREC), intent(in) :: tK

   integer :: i
   character(len=37) :: out_fmt

   ! Format in .out file
   if (flg_replica) then
      write(out_fmt, '(a23,i2,a12)') '(i12, 1x, i4, 1x, f6.2,', ENE%EXV+2, '(1x, g13.6))'
   else
      write(out_fmt, '(a15,i2,a12)') '(i12, 1x, f6.2,', ENE%EXV+2, '(1x, g13.6))'
   endif

   if (flg_replica) then
      write(hdl_out(irep), out_fmt, advance='no') istep_out, rep_label, tK, Ekinetic(irep), (energies(i, irep), i=0,ENE%EXV)
   else
      write(hdl_out(irep), out_fmt, advance='no') istep_out, tK, Ekinetic(irep), (energies(i, irep), i=0,ENE%EXV)
   endif

   if (flg_ele) then
      write(hdl_out(irep), '(1x, g13.6)', advance='no') energies(ENE%ELE, irep)
   endif

   if (flg_stage) then
      write(hdl_out(irep), '(1x, g13.6)', advance='no') energies(ENE%STAGE, irep)
   endif

   if (flg_twz) then
      write(hdl_out(irep), '(1x, g13.6)', advance='no') energies(ENE%TWZ, irep)
   endif

   write(hdl_out(irep), '(a)') ''

endsubroutine write_out
