set term postscript color 20 lw 4 solid
set out 'stage.ps'

set grid
set xr [5:20]
set yr [-3:15]

set yl 'U_{stage} / kcal/mol'
set xl 'z / A'

set key top right

d = 10.0
eps = 1.2

p 4 * eps * ((d/x)**12 - (d/x)**6) w l lw 2 notitle
