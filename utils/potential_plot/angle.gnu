set term pdfcairo color solid font "Arial, 12" lw 2
set out 'angle.pdf'

#kr = 10.0
t0 = 2.643

set grid
set xr [0:pi]
set yr [0:4]

set yl 'U_{angle} / kcal/mol'
set xl 'theta'

set key bottom left

set xtics (0, '1/4 pi' 0.25*pi, '1/2 pi' 0.5*pi, '3/4 pi' 0.75*pi, 'pi' pi)

set title 'Harmonic'
p 0.5 * 10.0 * (x - t0)**2 w l lw 2 title 'k = 10.0'\
, 0.5 *  5.0 * (x - t0)**2 w l lw 2 title 'k =  5.0'\
, 0.5 *  3.0 * (x - t0)**2 w l lw 2 title 'k =  3.0'

set title 'Restricted bending'
p 0.5 * 10.0 * (cos(x) - cos(t0))**2 / sin(x)**4 w l lw 2 title 'k = 10.0'\
, 0.5 *  5.0 * (cos(x) - cos(t0))**2 / sin(x)**4 w l lw 2 title 'k =  5.0'\
, 0.5 *  3.0 * (cos(x) - cos(t0))**2 / sin(x)**4 w l lw 2 title 'k =  3.0'
