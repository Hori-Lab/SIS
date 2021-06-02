module var_io

   logical, save :: flg_out_bp = .False.
   integer, parameter :: KIND_OUT_BP = 2

   integer, parameter :: hdl_dcd = 20
   integer, parameter :: hdl_out = 21
   integer, parameter :: hdl_bp  = 22

   integer, save :: iopen_hdl = 30

end module var_io
