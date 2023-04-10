#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color solid
set output "u0_u.ps"

set datafile separator ','
set grid
#set xrange [0:1]
#set yrange [0:1]
plot "x_u0.csv" using 1:2 with lines title "u0", "x_u.csv" using 1:2 with lines title "u"

set terminal postscript eps enhanced color solid
set output "x_p.ps"

set datafile separator ','
set grid
#set xrange [0:1]
#set yrange [0:1]
plot "x_p.csv" using 1:2 with lines title "p"

