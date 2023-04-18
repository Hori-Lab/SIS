set term pdfcairo color solid font "Arial, 12" lw 2
set out 'bond.pdf'

#kr = 15.0
r0 = 5.84

set grid
set xr [4:8]
set yr [0:6]

set yl 'U_{bond} / kcal/mol'
set xl 'r / A'

set key bottom right

p 0.5 * 15.0 * (x - r0)**2 w l lw 2 title 'k = 15.0'\
, 0.5 * 10.0 * (x - r0)**2 w l lw 2 title 'k = 10.0'\
, 0.5 *  5.0 * (x - r0)**2 w l lw 2 title 'k =  5.0'
