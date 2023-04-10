#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color solid
set output "u0_u.ps"

set datafile separator ','
set grid

set hidden3d

#set pm3d at bs
set contour

set dgrid3d 100,100 qnorm 3
set xlabel "x" font "Helvetica Bold ,18"
set ylabel "y" rotate by 90 font "Helvetica Bold ,18"
set zlabel "z" font "Helvetica Bold ,18"
splot "u0.csv" using 1:2:3 with lines title "u0", "u.csv" using 1:2:3 with lines title "u"

set terminal postscript eps enhanced color solid
set output "x_p.ps"

set datafile separator ','
set grid
set hidden3d
#set pm3d at bs
set contour

set dgrid3d 100,100 qnorm 3
set xlabel "x" font "Helvetica Bold ,18"
set ylabel "y" rotate by 90 font "Helvetica Bold ,18"
set zlabel "z" font "Helvetica Bold ,18"
splot "p.csv" using 1:2:3 with lines title "p"
