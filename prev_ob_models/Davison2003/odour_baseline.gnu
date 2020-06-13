# gnuplot file to produce a figure similar to Figure 8
# in Davison A.P., Feng J. and Brown D. (2003) J. Neurophysiol.

# Note that the figure produced will not match exactly the
# published figure due to differences in the sequence of
# random numbers used to set up the network.

# gnuplot is available from http://www.gnuplot.info

set term postscript portrait enhanced mono solid "Helvetica" 8
set output "odour_baseline.eps"

set size 0.49,0.5
set multiplot

set size 0.445,0.115
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set noxtics
set noytics
set nokey

set noborder
set origin 0,0.375
set label 1 '{/Helvetica=14 A}' at screen 0.01,0.475
plot [0:3000][-0.5:35.5] "odour_baseline.ras" using 4:3 with points pointtype 7 pointsize 0.3

set origin 0,0.25
set label 1 '{/Helvetica=14 B}' at screen 0.01,0.35
set border 1 linewidth 0.5
plot [0:3000] "odour_baseline.smhist" with lines linewidth 0.5

set origin 0,0.125
set noborder
set label 1 '{/Helvetica-Bold=14 C}' at screen 0.01,0.225
plot [0:3000] "odour_baseline.gran.ras" using 4:3 with dots

set origin 0,0.0
set border 1 linewidth 0.5
set label 1 '{/Helvetica-Bold=14 D}' at screen 0.01,0.10
set arrow from 200,70 to 700,70 nohead linewidth 2
set label 2 '500 ms' at 250,55
plot [0:3000] "odour_baseline.gran.smhist" with lines linewidth 0.5

set size 0.05,0.115
set origin 0.45,0.375
set nolabel 1
set border 2 lw 0.5
set noarrow
#plot "input1.nbar" using 1:($2+1) with steps linetype 2 linewidth 0.5
set arrow from 5,-1.5 to 25,-1.5 nohead lw 2
set label 2 '{/=8 10 s^{/=6 -1}}' at -1,-5
#plot [0:34][] "odour_baseline.nbar" using 1:($2+1.1) with steps linewidth 0.5

set origin 0.45,0.125
set noarrow
set nolabel
#plot [0:34][] "odour_baseline2.gran.nbar" using 1:($2+1) with steps linewidth 0.2
