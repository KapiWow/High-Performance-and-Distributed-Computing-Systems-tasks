set term svg enhanced mouse size 1600, 900 font "Helvetica,18"
set output "2_1_1.svg"
set key autotitle columnheader

set arrow from 8192, graph 0 to 8192, graph 1 nohead lw 3
set arrow from 65536, graph 0 to 65536, graph 1 nohead lw 3
set label "L1 cache" at 8192+2000, 3000
set label "L2 cache" at 65536+2000, 3000
set xlabel "size of array"
set ylabel "time"
set key right bottom

plot for [col=2:2] 'data2_1_1' using 1:col smooth bezier \
        title "bezier" lw 3, \
     for [col=2:2] 'data2_1_1' using 1:col with lines title "lines" lw 3
