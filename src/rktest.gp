#!/usr/bin/gnuplot
# -*- mode:sh -*-

set terminal pdf color enhanced size 4,3
set output "out/rk/rktest.pdf"
set title 'rktest: x(0) = 1 --[ {\\dot x = 1/x} ]--> x(1)'
set log x
set log y
set format x "10^{%L}"
set format y "10^{%L}"
set xlabel 'number of eval(\dot x)'
set ylabel 'abs(error)'
plot \
  "out/rk/rkeuler.txt" u 1:(abs($5)) title 'Euler', \
  "out/rk/rkmid.txt" u 1:(abs($5)) title 'midpoint', \
  "out/rk/rkrk4.txt" u 1:(abs($5)) title 'RK4', \
  "out/rk/rkcv8.txt" u 1:(abs($5)) title 'Cooper Verner 8'
