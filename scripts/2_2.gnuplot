set term svg enhanced mouse size 1600, 900 font "Helvetica,18"
set key autotitle columnheader
set output 'graph.svg'

set xlabel "size of array"
set ylabel "time"
set key left top
set logscale y

#plot for [col=2:4] 'data' using 1:col with lines title col lw 3
plot 'data' using 1:2 with lines title "serial" lw 3, \
     'data' using 1:3 with lines title "shift" lw 3, \
     'data' using 1:4 with lines title "random" lw 3, \
