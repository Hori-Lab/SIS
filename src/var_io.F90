module var_io
   
   use const

   logical, save :: flg_out_bp = .False.
   logical, save :: flg_out_bpe = .False.
   integer, parameter :: KIND_OUT_BP  = 2  ! Defines the format of bp output file.
   integer, parameter :: KIND_OUT_BPE = 4  ! Defines the format of bp output file.

   integer, parameter :: hdl_dcd = 20
   integer, parameter :: hdl_out = 21
   integer, parameter :: hdl_bp  = 22
   integer, parameter :: hdl_bpe = 23

   integer, save :: iopen_hdl = 30

   character(CHAR_FILE_PATH), save :: cfile_out
   character(CHAR_FILE_PATH), save :: cfile_bp
   character(len=:), allocatable, save :: cfile_ff, cfile_prefix, cfile_dcd_in, cfile_pdb_ini

end module var_io
