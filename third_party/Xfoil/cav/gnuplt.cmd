set title "Cavitation Diagram"
set ylabel "CL"
set xlabel "Cpmin"
set grid
set xrange [0:-6]
set xtics 0.5
set yrange [-0.5:1.5]
set ytics 0.5
 
plot "cavplt.tmp1" using 6:2 w lines lw 2 t "E817mod1 10%"\
, "cavplt.tmp2" using 6:2 w lines lw 2 t "H105-T10"\
, "cavplt.tmp3" using 6:2 w lines lw 2 t "DB10C"\
 
set term postscript color
set output "cav-plt.ps"
replot
