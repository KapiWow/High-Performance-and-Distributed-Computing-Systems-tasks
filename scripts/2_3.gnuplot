set term svg enhanced mouse size 1600, 900 font "Helvetica,18"
set key autotitle columnheader
set output 'graph.svg'

set xlabel "size of array"
set ylabel "time"
set key left top

plot for [col=2:9] 'data' using 1:col with lines title col lw 3
