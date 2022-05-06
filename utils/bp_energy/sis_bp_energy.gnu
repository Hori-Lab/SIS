set term postscript color 20 lw 4 solid
set out 'sis_bp_energy.ps'

set title 'Distance dependence'
set xl 'r / A'
set yl 'Ubp / kcal/mol'
set xr [8:20]
p 'bp_energy_NHT22.out' i 0 u 1:2 w l title 'NHT22' \
, 'bp_energy_para3.out' i 0 u 1:2 w l title 'PARA3'

set title 'Angle dependence'
set xl 'theta / degree'
set xr [*:*]
set yl 'Ubp / kcal/mol'
set key bottom right
p 'bp_energy_NHT22.out' i 1 u 1:2 w l title 'NHT22' \
, 'bp_energy_para3.out' i 1 u 1:2 w l title 'PARA3'

set title 'Dihedral dependence'
set xl 'phi/ degree'
set yl 'Ubp / kcal/mol'
set key bottom right
p 'bp_energy_NHT22.out' i 2 u 1:2 w l title 'NHT22' \
, 'bp_energy_para3.out' i 2 u 1:2 w l title 'PARA3'
