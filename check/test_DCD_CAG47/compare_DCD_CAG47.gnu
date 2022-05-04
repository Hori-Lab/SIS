set term postscript color 20 lw 2 solid
set out 'compare_DCD_CAG47.ps'

set title 'U_{bond}'
p './CAG47.energy.out' u 1:2 w lp lw 2 title 'OpenMM'\
, './test_DCD_CAG47.out' u ($1*100):5 w lp title 'sis'\
, './test_DCD_CAG47_HTN.out' u ($1*100):5 w lp title 'sis(HTN)'

set title 'U_{angl}'
p './CAG47.energy.out' u 1:3 w lp lw 2 title 'OpenMM'\
, './test_DCD_CAG47.out' u ($1*100):6 w lp title 'sis'\
, './test_DCD_CAG47_HTN.out' u ($1*100):6 w lp title 'sis(HTN)'

set title 'U_{bp}'
p './CAG47.energy.out' u 1:7 w lp lw 2 title 'OpenMM'\
, './test_DCD_CAG47.out' u ($1*100):7 w lp title 'sis'\
, './test_DCD_CAG47_HTN.out' u ($1*100):7 w lp title 'sis(HTN)'

set title 'U_{exv}'
p './CAG47.energy.out' u 1:4 w lp lw 2 title 'OpenMM'\
, './test_DCD_CAG47.out' u ($1*100):8 w lp title 'sis'\
, './test_DCD_CAG47_HTN.out' u ($1*100):8 w lp title 'sis(HTN)'

set title 'U_{ele}'
p './CAG47.energy.out' u 1:5 w lp lw 2 title 'OpenMM'\
, './test_DCD_CAG47.out' u ($1*100):9 w lp title 'sis'\
, './test_DCD_CAG47_HTN.out' u ($1*100):9 w lp title 'sis(HTN)'

set title 'U_{total}'
p './CAG47.energy.out' u 1:($2+$3+$7+$4+$5) w lp lw 2 title 'OpenMM'\
, './test_DCD_CAG47.out' u ($1*100):4 w lp title 'sis'\
, './test_DCD_CAG47_HTN.out' u ($1*100):4 w lp title 'sis(HTN)'
