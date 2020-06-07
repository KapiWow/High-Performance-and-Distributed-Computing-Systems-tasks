set term svg enhanced mouse size 1600, 900 font "Helvetica,18"
set output "2_1_2.svg"
set key autotitle columnheader

set arrow from 517152, graph 0 to 517152, graph 1 nohead lw 3
set arrow from 1017152, graph 0 to 1017152, graph 1 nohead lw 3
set arrow from 2097152, graph 0 to 2097152, graph 1 nohead lw 3
set label "L3 cache" at 2097152+20000, 10000 front
set label "1/2 L3 cache" at 1017152+20000, 10000 front
set label "1/4 L3 cache" at 517152+20000, 10000 front
set xlabel "size of array"
set ylabel "time"
set key right bottom

plot for [col=2:2] 'data2_1_2' using 1:col smooth bezier \
        title "bezier" lw 2, \
     for [col=2:2] 'data2_1_2' using 1:col with lines title "lines" lw 2
